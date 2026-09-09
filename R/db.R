## Query layer over the SQLite lipid library (see scripts/lbm2_to_sqlite.py).

library(DBI)
library(RSQLite)

#' Open a read-only connection to the lipid library.
lipid_db <- function(path = "DB/lipids.sqlite") {
  if (!file.exists(path)) {
    stop("Lipid database not found at '", path, "'.\n",
         "Build it first:\n",
         "  python3 scripts/lbm2_to_sqlite.py DB/<file>.lbm2 ", path,
         call. = FALSE)
  }
  dbConnect(RSQLite::SQLite(), path, flags = RSQLite::SQLITE_RO)
}

#' Lipid classes with record counts, m/z range and retention-time range.
lipid_classes <- function(con) {
  dbGetQuery(con, "SELECT lipid_class, n, min_mz, max_mz, min_rt, max_rt
                     FROM lipid_class_summary ORDER BY lipid_class")
}

#' Adducts with ion mode and record counts.
lipid_adducts <- function(con) {
  dbGetQuery(con, "SELECT adduct, ion_mode, n FROM adduct_summary
                    WHERE adduct <> '' ORDER BY n DESC")
}

#' Total number of records in the library.
lipid_count <- function(con) {
  sum(dbGetQuery(con, "SELECT n FROM lipid_class_summary")$n)
}

# Helper: "?,?,?" placeholder list for an IN clause.
in_clause <- function(x) paste0("(", paste(rep("?", length(x)), collapse = ","), ")")

#' Search the library by any combination of name, precursor m/z, class,
#' ion mode and adduct.
#'
#' @param name      Substring of the lipid name, e.g. "18:1/16:0".
#' @param mz,tol    Precursor m/z and half-window, in Da.
#' @param class     Zero or more exact lipid classes.
#' @param mode      "Positive" or "Negative" (NULL for either).
#' @param adduct    Zero or more exact adduct strings.
#' @param limit     Maximum rows returned.
#'
#' Name matching picks one of two strategies, which matters a lot:
#'  * With a class or m/z filter present, LIKE is driven off ix_lipid_class_mz
#'    and only scans that slice (~20 ms). Routing through the trigram index
#'    instead would materialise every matching name first -- "18:1" alone hits
#'    156k names -- and costs ~315 ms.
#'  * With the name as the only filter, LIKE degrades to a full scan for rare
#'    substrings (~410 ms), so the trigram index wins (~1 ms).
#' The trigram tokeniser cannot match patterns shorter than 3 characters, so
#' those always fall back to LIKE.
search_lipids <- function(con, name = NULL, mz = NULL, tol = 0.01,
                          class = NULL, mode = NULL, adduct = NULL,
                          limit = 500) {
  where <- character()
  params <- list()

  has_name <- !is.null(name) && nzchar(name)
  selective <- !is.null(mz) || length(class) > 0

  if (has_name && nchar(name) >= 3 && !selective) {
    where <- c(where, "l.name IN (SELECT name FROM name_fts WHERE name MATCH ?)")
    params <- c(params, list(paste0('"', gsub('"', '""', name), '"')))
  } else if (has_name) {
    where <- c(where, "l.name LIKE ?")
    params <- c(params, list(paste0("%", name, "%")))
  }

  if (!is.null(mz)) {
    where <- c(where, "l.precursor_mz BETWEEN ? AND ?")
    params <- c(params, list(mz - tol, mz + tol))
  }
  if (length(class) > 0) {
    where <- c(where, paste("l.lipid_class IN", in_clause(class)))
    params <- c(params, as.list(class))
  }
  if (!is.null(mode) && nzchar(mode)) {
    where <- c(where, "l.ion_mode = ?")
    params <- c(params, list(mode))
  }
  if (length(adduct) > 0) {
    where <- c(where, paste("l.adduct IN", in_clause(adduct)))
    params <- c(params, as.list(adduct))
  }

  sql <- paste0(
    "SELECT l.id, l.name, l.lipid_class, l.precursor_mz, l.adduct,
            l.ion_mode, l.retention_time, l.ccs, l.formula, l.n_peaks
       FROM lipid l",
    if (length(where)) paste0(" WHERE ", paste(where, collapse = " AND ")) else "",
    " ORDER BY l.precursor_mz LIMIT ?")

  dbGetQuery(con, sql, params = c(params, list(limit)))
}

#' Fetch one record's MS/MS spectrum as a data.frame of mz / intensity.
#'
#' Peaks are stored as a BLOB of little-endian float32 pairs; the source
#' values are float32, so this round-trips exactly.
get_spectrum <- function(con, id) {
  r <- dbGetQuery(con, "SELECT id, name, lipid_class, precursor_mz, adduct,
                               ion_mode, retention_time, ccs, formula, exact_mass,
                               inchikey, smiles, n_peaks, peaks
                          FROM lipid WHERE id = ?", params = list(id))
  if (nrow(r) == 0L) return(NULL)
  v <- readBin(r$peaks[[1]], "numeric", n = r$n_peaks * 2L, size = 4L)
  peaks <- data.frame(mz = v[c(TRUE, FALSE)], intensity = v[c(FALSE, TRUE)])
  peaks <- peaks[order(peaks$mz), , drop = FALSE]
  peaks$rel <- 100 * peaks$intensity / max(peaks$intensity, 1)
  list(meta = r[, setdiff(names(r), "peaks")], peaks = peaks)
}
