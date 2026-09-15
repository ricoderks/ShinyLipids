## Whole-library MS/MS search: match a pasted query spectrum against every
## record in DB/lipids.sqlite.
##
## The library holds ~5.3M spectra and ~35M peaks, so the search cannot walk
## the records. Instead every peak in the library is pulled into one flat
## array sorted by m/z; a query peak then resolves to a contiguous slice of
## that array with two binary searches, and a 30-peak query touches tens of
## thousands of peaks rather than tens of millions. Building the index costs
## ~10 s and ~1 GB of RAM once per process; every search after that is well
## under a second.
##
## See R/db.R for the record-level query layer this sits beside.

library(data.table)

# The library's peaks are the only thing worth parallelising here; data.table
# picks a sensible default but a container may report a single core.
setDTthreads(0)

# Stein & Scott (1994) weighting for the weighted dot product: w = mz^a * I^b.
# Their optimum for library matching was a = 3, b = 0.6, which is what NIST
# and MS-DIAL inherited; MassBank-style scoring uses a = 2, b = 0.5 instead,
# so both are exposed rather than baked in.
MZ_POWER  <- 3
INT_POWER <- 0.6

# ------------------------------------------------------------ parsing -------

#' Parse a pasted peak list into a spectrum.
#'
#' Accepts whatever a user is likely to paste: comma, tab, semicolon or
#' whitespace between the two columns, with or without a header, and with the
#' surrounding metadata of an MSP block still attached. Any line whose first
#' field is not a number is skipped, which disposes of headers ("m/z\tint"),
#' MSP fields ("Num Peaks: 12") and blank lines in one rule.
#'
#' A third column, if present, is ignored -- MSP peak lists often carry an
#' annotation string there.
#'
#' @return list(peaks = data.frame(mz, intensity, rel), n_skipped) where `rel`
#'   is intensity as a percentage of the base peak, matching get_spectrum().
#'   `peaks` has zero rows if nothing parsed.
parse_peak_text <- function(txt) {
  empty <- list(peaks = data.frame(mz = numeric(0), intensity = numeric(0),
                                   rel = numeric(0)),
                n_skipped = 0L)
  if (is.null(txt) || !nzchar(trimws(txt))) return(empty)

  lines <- strsplit(txt, "\r\n|\r|\n")[[1]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (!length(lines)) return(empty)

  # Split on a comma / semicolon / tab / run of spaces, each optionally padded.
  fields <- strsplit(lines, "[ \t]*[,;\t][ \t]*|[ \t]+")

  mz  <- suppressWarnings(as.numeric(vapply(fields, function(f)
    if (length(f) >= 1L) f[1] else NA_character_, character(1))))
  int <- suppressWarnings(as.numeric(vapply(fields, function(f)
    if (length(f) >= 2L) f[2] else NA_character_, character(1))))

  ok <- !is.na(mz) & !is.na(int) & is.finite(mz) & is.finite(int) &
        mz > 0 & int > 0
  n_skipped <- sum(!ok)
  if (!any(ok)) return(list(peaks = empty$peaks, n_skipped = n_skipped))

  p <- data.frame(mz = mz[ok], intensity = int[ok])
  p <- p[order(p$mz), , drop = FALSE]
  p$rel <- 100 * p$intensity / max(p$intensity)
  rownames(p) <- NULL
  list(peaks = p, n_skipped = n_skipped)
}

#' Wrap a parsed peak list so it can be handed to the spectrum plotting code.
#'
#' The plot helpers expect the shape get_spectrum() returns. A pasted spectrum
#' has no precursor and no library metadata, so precursor_mz is whatever the
#' user typed (possibly NA) and the plot skips the precursor marker for it.
as_query_spectrum <- function(peaks, precursor_mz = NA_real_,
                              name = "Query spectrum") {
  list(meta = data.frame(name = name, precursor_mz = as.numeric(precursor_mz),
                         stringsAsFactors = FALSE),
       peaks = peaks)
}

# ------------------------------------------------------- match tolerance ----

#' Half-window in Da for a given m/z, from a Da and a ppm setting.
#'
#' Both are honoured at once and the wider of the two wins, so a single pair of
#' inputs covers "0.01 Da", "10 ppm" and "10 ppm but never tighter than
#' 0.005 Da". Setting either to 0 (or NA) disables that half.
mz_tolerance <- function(mz, tol_da = 0.01, tol_ppm = 0) {
  da  <- if (isTRUE(is.finite(tol_da)))  tol_da  else 0
  ppm <- if (isTRUE(is.finite(tol_ppm))) tol_ppm else 0
  pmax(da, mz * ppm / 1e6)
}

# --------------------------------------------------------- peak index -------

# The index is expensive to build and immutable once built, so it is cached
# for the lifetime of the R process and shared by every Shiny session.
.index_cache <- new.env(parent = emptyenv())

#' Build the flat, m/z-sorted index of every peak in the library.
#'
#' Norms for the plain dot product are precomputed here; the weighted norms
#' depend on the exponents the user picks, so those are computed on demand and
#' memoised by peak_norms_weighted().
#'
#' @param progress Optional function(fraction, message) for a progress bar.
#' @return An environment (cheap to pass around) with:
#'   n_rec, id, precursor_mz, adduct, ion_mode   -- one entry per record
#'   norm_plain                                  -- one entry per record
#'   mz, rel, rec                                -- one entry per peak, m/z order
build_peak_index <- function(con, progress = NULL) {
  say <- function(f, msg) if (is.function(progress)) progress(f, msg)

  say(0.05, "Reading spectra from the library…")
  r <- DBI::dbGetQuery(con, "SELECT id, precursor_mz, adduct, ion_mode, n_peaks,
                                    peaks FROM lipid ORDER BY id")

  # One readBin over every blob concatenated is ~30x faster than decoding each
  # record separately; the blobs are all little-endian float32 pairs, so they
  # concatenate without a separator.
  say(0.35, "Decoding peak lists…")
  raw <- unlist(r$peaks, use.names = FALSE)
  v <- readBin(raw, "numeric", n = length(raw) %/% 4L, size = 4L)
  rm(raw); r$peaks <- NULL; invisible(gc(FALSE))

  mz  <- v[c(TRUE, FALSE)]
  int <- v[c(FALSE, TRUE)]
  rm(v); invisible(gc(FALSE))

  n_rec <- nrow(r)
  rec <- rep.int(seq_len(n_rec), r$n_peaks)

  # Relative to each record's base peak, so library intensities are comparable
  # with the query's, which parse_peak_text() normalises the same way.
  say(0.60, "Normalising intensities…")
  D <- data.table(rec = rec, int = int)
  base <- D[, .(m = max(int)), by = rec]$m
  rel <- int / base[rec]
  rm(D, int, base); invisible(gc(FALSE))

  norm_plain <- sqrt(data.table(g = rec, x = rel * rel)[, sum(x), keyby = g]$V1)

  say(0.80, "Sorting peaks by m/z…")
  ord <- order(mz, method = "radix")

  ix <- new.env(parent = emptyenv())
  ix$n_rec        <- n_rec
  ix$id           <- as.integer(r$id)
  ix$precursor_mz <- r$precursor_mz
  ix$adduct       <- r$adduct
  ix$ion_mode     <- r$ion_mode
  ix$norm_plain   <- norm_plain
  ix$mz  <- mz[ord]
  ix$rel <- rel[ord]
  ix$rec <- rec[ord]
  ix$n_peaks_total <- length(ord)
  ix$wnorm <- new.env(parent = emptyenv())   # memoised weighted norms
  rm(mz, rel, rec, ord); invisible(gc(FALSE))

  # Warm the norms for the default exponents: they are ~1.5 s over 35M peaks,
  # and paying that here keeps the first search as fast as every later one.
  say(0.92, "Precomputing weighted norms\u2026")
  peak_norms_weighted(ix, MZ_POWER, INT_POWER)

  say(1, "Index ready")
  ix
}

#' The cached peak index, built on first use.
peak_index <- function(con, progress = NULL) {
  if (is.null(.index_cache$ix)) .index_cache$ix <- build_peak_index(con, progress)
  .index_cache$ix
}

#' Whether the index has already been built (so the UI can warn about the wait).
peak_index_ready <- function() !is.null(.index_cache$ix)

#' Per-record norms of the weighted intensities, memoised per exponent pair.
#'
#' Recomputing these is ~1.5 s over 35M peaks, which is the whole cost of a
#' first search with new exponents; every later search with the same pair is
#' free. Users change the exponents rarely, if ever, so the cache effectively
#' never misses.
peak_norms_weighted <- function(ix, mz_power, int_power) {
  key <- sprintf("%.6g_%.6g", mz_power, int_power)
  if (is.null(ix$wnorm[[key]])) {
    w <- ix$mz^mz_power * ix$rel^int_power
    ix$wnorm[[key]] <-
      sqrt(data.table(g = ix$rec, x = w * w)[, sum(x), keyby = g]$V1)
  }
  ix$wnorm[[key]]
}

# ------------------------------------------------------------- search -------

#' Score a query spectrum against every record the filters allow.
#'
#' Three scores are reported, all as cosine similarities on 0-1 so they sit on
#' the same scale as the mirror plot's score:
#'
#'  * dot          -- cosine on relative intensities. Unmatched peaks on either
#'                    side contribute nothing to the numerator but still count
#'                    in both norms, so extra fragments are penalised.
#'  * weighted dot -- the same, on w = mz^a * I^b (Stein & Scott). Upweighting
#'                    high-m/z, high-intensity fragments suppresses matches
#'                    that rest on a crowd of small, unspecific low-mass peaks.
#'  * reverse dot  -- cosine on relative intensities with query peaks that have
#'                    no library counterpart dropped from the query norm. This
#'                    is the classical reverse match: it ignores co-isolated
#'                    contamination in the query, at the cost of rewarding
#'                    library records with very few peaks.
#'
#' Peaks are paired greedily, strongest intensity product first, so a peak on
#' either side is used at most once -- the same rule as match_peaks(), and the
#' one matchms calls "cosine greedy".
#'
#' @param q          data.frame(mz, rel) -- the query peaks (rel as % of base).
#' @param mode       "Positive"/"Negative", or "" / NULL for either.
#' @param adduct     Zero or more exact adduct strings; empty means any.
#' @param precursor  Optional query precursor m/z; records outside tolerance
#'                   of it are dropped.
#' @param min_match  Minimum number of matched peaks for a record to be kept.
#' @param rank_by    "wdot", "dot" or "rdot" -- which score orders the result.
#' @param top_n      Rows returned.
#' @return data.frame(id, precursor_mz, adduct, ion_mode, n_match, dot, wdot,
#'   rdot), ordered by `rank_by` descending. Zero rows if nothing matched.
search_spectrum <- function(ix, q, tol_da = 0.01, tol_ppm = 0,
                            mode = NULL, adduct = NULL, precursor = NULL,
                            min_match = 1L, mz_power = MZ_POWER,
                            int_power = INT_POWER, rank_by = "wdot",
                            top_n = 200L) {
  none <- data.frame(id = integer(0), precursor_mz = numeric(0),
                     adduct = character(0), ion_mode = character(0),
                     n_match = integer(0), dot = numeric(0),
                     wdot = numeric(0), rdot = numeric(0))
  if (nrow(q) == 0L) return(with_counts(none, NA, NA))

  q_mz  <- q$mz
  q_rel <- q$rel / 100                       # base peak = 1, as in the index
  q_wp  <- q_rel                             # plain weights
  q_ww  <- q_mz^mz_power * q_rel^int_power   # weighted weights
  nq_p  <- sqrt(sum(q_wp * q_wp))
  nq_w  <- sqrt(sum(q_ww * q_ww))
  if (nq_p == 0 || nq_w == 0) return(with_counts(none, NA, NA))

  # --- record-level filters -------------------------------------------------
  keep <- keep_mask(ix, mode, adduct, precursor, tol_da, tol_ppm)
  n_filtered <- if (is.null(keep)) ix$n_rec else sum(keep)
  if (n_filtered == 0) return(with_counts(none, 0, 0))

  # --- candidate peaks: two binary searches per query peak ------------------
  tol <- mz_tolerance(q_mz, tol_da, tol_ppm)
  lo <- findInterval(q_mz - tol, ix$mz)
  hi <- findInterval(q_mz + tol, ix$mz)
  cnt <- hi - lo
  if (sum(cnt) == 0) return(with_counts(none, n_filtered, 0))

  pidx <- sequence(cnt, from = lo + 1L)      # positions in the sorted arrays
  qidx <- rep.int(seq_along(q_mz), cnt)      # which query peak each came from

  H <- data.table(rec = ix$rec[pidx], qi = qidx, li = pidx,
                  lrel = ix$rel[pidx])
  if (!is.null(keep)) {
    H <- H[keep[H$rec]]
    if (nrow(H) == 0L) return(with_counts(none, n_filtered, 0))
  }

  # --- greedy one-to-one pairing -------------------------------------------
  # Strongest intensity product first, then drop any later pair that reuses a
  # peak on either side. Sorting by -prod inside each record and taking the
  # first occurrence per (record, query peak) and then per (record, library
  # peak) is exactly the greedy rule, without a loop.
  H[, prod := q_wp[qi] * lrel]
  setorder(H, rec, -prod)
  H <- H[!duplicated(H, by = c("rec", "qi"))]
  H <- H[!duplicated(H, by = c("rec", "li"))]

  # --- per-record sums ------------------------------------------------------
  H[, `:=`(num_w = q_ww[qi] * (ix$mz[li]^mz_power * lrel^int_power),
           mq_p  = q_wp[qi] * q_wp[qi])]
  S <- H[, .(n_match = .N, num_p = sum(prod), num_w = sum(num_w),
             mq_p = sum(mq_p)), keyby = rec]
  n_hit <- nrow(S)
  if (min_match > 1L) S <- S[n_match >= min_match]
  if (nrow(S) == 0L) return(with_counts(none, n_filtered, n_hit))

  nl_p <- ix$norm_plain[S$rec]
  nl_w <- peak_norms_weighted(ix, mz_power, int_power)[S$rec]

  S[, `:=`(dot  = num_p / (nq_p * nl_p),
           wdot = num_w / (nq_w * nl_w),
           rdot = num_p / (sqrt(mq_p) * nl_p))]

  key <- switch(rank_by, dot = S$dot, rdot = S$rdot, S$wdot)
  # Ties on the score are common -- isomeric records share a spectrum exactly
  # -- so break them on the match count to keep the order stable.
  o <- order(-key, -S$n_match, method = "radix")
  o <- utils::head(o, top_n)

  with_counts(data.frame(
    id           = ix$id[S$rec[o]],
    precursor_mz = ix$precursor_mz[S$rec[o]],
    adduct       = ix$adduct[S$rec[o]],
    ion_mode     = ix$ion_mode[S$rec[o]],
    n_match      = as.integer(S$n_match[o]),
    dot          = S$dot[o],
    wdot         = S$wdot[o],
    rdot         = S$rdot[o],
    stringsAsFactors = FALSE
  ), n_filtered, n_hit)
}

#' Which library records the polarity / adduct / precursor filters allow.
#'
#' @return A logical vector over records, or NULL when nothing is filtered --
#'   the caller then skips the subset entirely rather than paying for a
#'   5.3M-element all-TRUE mask.
keep_mask <- function(ix, mode = NULL, adduct = NULL, precursor = NULL,
                      tol_da = 0.01, tol_ppm = 0) {
  keep <- NULL
  and <- function(x) if (is.null(keep)) x else keep & x

  if (!is.null(mode) && length(mode) == 1L && nzchar(mode)) {
    keep <- and(ix$ion_mode == mode)
  }
  if (length(adduct) > 0) keep <- and(ix$adduct %in% adduct)
  if (!is.null(precursor) && length(precursor) == 1L && is.finite(precursor)) {
    keep <- and(abs(ix$precursor_mz - precursor) <=
                  mz_tolerance(precursor, tol_da, tol_ppm))
  }
  keep
}

#' Tag a result frame with how wide the search actually was.
#'
#' `n_filtered` is how many library records the polarity/adduct/precursor
#' filters allowed, `n_hit` how many of those shared at least one peak with the
#' query before `min_match` was applied. Attributes rather than columns so the
#' frame still goes straight into a table or a CSV.
with_counts <- function(df, n_filtered, n_hit) {
  attr(df, "n_filtered") <- as.numeric(n_filtered)
  attr(df, "n_hit") <- as.numeric(n_hit)
  df
}

#' The same three scores for a single query/library pair.
#'
#' Used by the mirror plot so the numbers under the plot are produced by the
#' same rules as the hit table rather than recomputed a second way.
#'
#' @return list(match = list(a, b) index vectors into q / lib, dot, wdot, rdot)
score_pair <- function(q, lib, tol_da = 0.01, tol_ppm = 0,
                       mz_power = MZ_POWER, int_power = INT_POWER) {
  none <- list(match = list(a = integer(0), b = integer(0)),
               dot = NA_real_, wdot = NA_real_, rdot = NA_real_)
  if (nrow(q) == 0L || nrow(lib) == 0L) return(none)

  qa <- q$rel / 100
  lb <- lib$rel / 100
  q_ww <- q$mz^mz_power * qa^int_power
  l_ww <- lib$mz^mz_power * lb^int_power

  tol <- mz_tolerance(q$mz, tol_da, tol_ppm)
  pairs <- do.call(rbind, lapply(seq_along(q$mz), function(i) {
    j <- which(abs(lib$mz - q$mz[i]) <= tol[i])
    if (!length(j)) NULL else data.frame(a = i, b = j, prod = qa[i] * lb[j])
  }))
  if (is.null(pairs) || nrow(pairs) == 0L) return(none)

  pairs <- pairs[order(-pairs$prod), , drop = FALSE]
  pairs <- pairs[!duplicated(pairs$a), , drop = FALSE]
  pairs <- pairs[!duplicated(pairs$b), , drop = FALSE]

  nq_p <- sqrt(sum(qa^2)); nl_p <- sqrt(sum(lb^2))
  nq_w <- sqrt(sum(q_ww^2)); nl_w <- sqrt(sum(l_ww^2))
  num_p <- sum(qa[pairs$a] * lb[pairs$b])
  num_w <- sum(q_ww[pairs$a] * l_ww[pairs$b])

  list(match = list(a = pairs$a, b = pairs$b),
       dot  = num_p / (nq_p * nl_p),
       wdot = num_w / (nq_w * nl_w),
       rdot = num_p / (sqrt(sum(qa[pairs$a]^2)) * nl_p))
}
