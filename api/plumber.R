## Read-only REST API over the MS-DIAL lipid library.
##
## Runs as its own process alongside the Shiny app, both reading
## DB/lipids.sqlite and sharing the query layer in R/db.R. It is deliberately
## not layered behind Shiny: read-only SQLite is just a shared file, so the two
## are peers and either can restart without affecting the other.
##
## Start with:  Rscript api/run.R          (see that file for host/port options)
## Docs at:     http://<host>:<port>/__docs__/

library(plumber)
library(DBI)

# plumber sets the working directory to this file's own folder while sourcing
# it, so every path here is anchored to the project root instead. run.R exports
# it; the ".." fallback covers running this file directly.
ROOT <- Sys.getenv("LIPID_PROJECT_ROOT", unset = "..")

source(file.path(ROOT, "R", "db.R"))
source(file.path(ROOT, "R", "spectrum.R"))

# Data lives outside the code tree so the container image carries none of it:
# the database is 2.4 GB and is mounted read-only at runtime.
DB_PATH <- Sys.getenv("LIPID_DB", unset = file.path(ROOT, "DB", "lipids.sqlite"))

con <- lipid_db(DB_PATH)

CLASSES <- lipid_classes(con)
ADDUCTS <- lipid_adducts(con)

# Hard ceiling on any single response. The library holds 5.3M records, so an
# uncapped limit is a denial-of-service on yourself.
MAX_LIMIT <- 1000L

# Records whose name and stored formula contradict each other (see
# scripts/export_mismatches.R). They are reported with flagged = true rather
# than hidden: silently serving them is how bad data spreads, and silently
# dropping them would make counts inexplicable to a caller.
FLAGGED <- local({
  f <- Sys.getenv("LIPID_FLAGGED",
                  unset = file.path(ROOT, "name_formula_mismatches.csv"))
  if (!file.exists(f)) {
    warning("no ", f, "; every record will report flagged = false", call. = FALSE)
    return(integer(0))
  }
  sort(unique(read.csv(f, stringsAsFactors = FALSE)$id))
})

# Split "PC,PE" or a repeated ?class=PC&class=PE into a character vector.
as_multi <- function(x) {
  if (is.null(x)) return(NULL)
  x <- unlist(strsplit(as.character(x), ",", fixed = TRUE))
  x <- trimws(x)
  x <- x[nzchar(x)]
  if (length(x)) x else NULL
}

# Drop NULL entries: jsonlite renders a NULL list element as "{}", so echoing
# unset filters back would litter the response with empty objects.
compact <- function(x) x[!vapply(x, is.null, logical(1))]

as_num <- function(x) {
  if (is.null(x) || identical(x, "")) return(NULL)
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) NULL else v
}

decorate <- function(df) {
  df$flagged <- df$id %in% FLAGGED
  df
}

#* @apiTitle ShinyLipids API
#* @apiDescription Read-only access to an MS-DIAL in-silico lipid MS/MS
#*   library: 5,275,482 records across 203 lipid classes, searchable by name,
#*   precursor m/z, class, ion mode and adduct, with MS/MS peak lists.
#*
#*   Records whose shorthand name contradicts their stored molecular formula
#*   are returned with `flagged: true` rather than removed.
#* @apiVersion 1.0.0

#* Allow browser clients. The API is read-only, so a permissive origin is safe.
#* `forward` is not a parameter -- plumber would try to fill it from the
#* request -- it is called as plumber::forward().
#* @filter cors
function(req, res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  res$setHeader("Access-Control-Allow-Methods", "GET, OPTIONS")
  if (identical(req$REQUEST_METHOD, "OPTIONS")) {
    res$status <- 200L
    return(list())
  }
  plumber::forward()
}

#* Service health and library size
#* @get /health
#* @serializer json list(digits = 8, na = "null", auto_unbox = TRUE)
function() {
  list(
    status  = "ok",
    records = sum(CLASSES$n),
    classes = nrow(CLASSES),
    adducts = nrow(ADDUCTS),
    flagged_records = length(FLAGGED)
  )
}

#* Lipid classes, with record counts and m/z range
#* @get /classes
#* @serializer json list(digits = 8, na = "null", auto_unbox = TRUE)
function() CLASSES

#* Adducts, with ion mode and record counts
#* @get /adducts
#* @serializer json list(digits = 8, na = "null", auto_unbox = TRUE)
function() ADDUCTS

#* Search the library
#*
#* Combine any of the filters. `class` and `adduct` accept a comma-separated
#* list. Results are ordered by precursor m/z.
#*
#* @param name:str Substring of the lipid name, e.g. "18:0_20:4"
#* @param mz:number Precursor m/z
#* @param tol:number Tolerance around mz (default 0.01)
#* @param tol_unit:str "da" (default) or "ppm"
#* @param class:str Lipid class, comma-separated for several
#* @param mode:str Ion mode: "Positive" or "Negative"
#* @param adduct:str Adduct, comma-separated for several
#* @param limit:int Maximum records to return (default 100, max 1000)
#* @get /lipids
#* @serializer json list(digits = 8, na = "null", auto_unbox = TRUE)
function(name = "", mz = "", tol = 0.01, tol_unit = "da", class = "",
         mode = "", adduct = "", limit = 100L, res) {

  mzv <- as_num(mz)
  tolv <- as_num(tol)
  if (is.null(tolv) || tolv < 0) tolv <- 0.01
  if (identical(tolower(tol_unit), "ppm") && !is.null(mzv)) {
    tolv <- mzv * tolv / 1e6
  }

  lim <- suppressWarnings(as.integer(limit))
  if (is.na(lim) || lim < 1L) lim <- 100L
  if (lim > MAX_LIMIT) lim <- MAX_LIMIT

  modev <- if (nzchar(mode)) mode else NULL
  if (!is.null(modev) && !modev %in% c("Positive", "Negative")) {
    res$status <- 400L
    return(list(error = "mode must be 'Positive' or 'Negative'"))
  }

  # Fetch one extra row to report truncation without paying for a COUNT(*).
  hits <- search_lipids(con,
    name = if (nzchar(name)) name else NULL,
    mz = mzv, tol = tolv,
    class = as_multi(class), mode = modev, adduct = as_multi(adduct),
    limit = lim + 1L)

  truncated <- nrow(hits) > lim
  hits <- utils::head(hits, lim)

  list(
    query = compact(list(
      name = if (nzchar(name)) name else NULL, mz = mzv, tol_da = tolv,
      class = as_multi(class), mode = modev, adduct = as_multi(adduct),
      limit = lim)),
    returned = nrow(hits),
    truncated = truncated,
    results = decorate(hits)
  )
}

#* One record by id
#* @param id:int Record id
#* @get /lipids/<id:int>
#* @serializer json list(digits = 8, na = "null", auto_unbox = TRUE)
function(id, res) {
  s <- get_spectrum(con, id)
  if (is.null(s)) {
    res$status <- 404L
    return(list(error = "no such record", id = id))
  }
  m <- decorate(s$meta)
  m$n_peaks <- nrow(s$peaks)
  # as.list so a single record serialises as an object, not a 1-element array.
  as.list(m)
}

#* MS/MS peak list for one record
#*
#* Peaks are returned in m/z order. `intensity` is as stored; `rel` is the
#* percentage of the base peak.
#* @param id:int Record id
#* @get /lipids/<id:int>/spectrum
#* @serializer json list(digits = 8, na = "null", auto_unbox = TRUE)
function(id, res) {
  s <- get_spectrum(con, id)
  if (is.null(s)) {
    res$status <- 404L
    return(list(error = "no such record", id = id))
  }
  list(
    id = id,
    name = s$meta$name,
    precursor_mz = s$meta$precursor_mz,
    adduct = s$meta$adduct,
    ion_mode = s$meta$ion_mode,
    flagged = id %in% FLAGGED,
    n_peaks = nrow(s$peaks),
    peaks = s$peaks[, c("mz", "intensity", "rel")]
  )
}
