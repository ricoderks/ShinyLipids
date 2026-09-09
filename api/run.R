## Launch the lipid library REST API.
##
##   Rscript api/run.R                    # 127.0.0.1:7413, local only
##   LIPID_API_HOST=0.0.0.0 Rscript api/run.R
##   LIPID_API_PORT=8000 Rscript api/run.R
##
## Must be run from the project root: api/plumber.R sources R/db.R and opens
## DB/lipids.sqlite by relative path.
##
## Binding to 127.0.0.1 is the default on purpose. The API has no
## authentication and no rate limiting, so exposing it on 0.0.0.0 puts an
## unauthenticated reader of a 2.4 GB database on the network -- do that only
## behind a reverse proxy that adds whatever access control you need.

if (!file.exists("api/plumber.R")) {
  stop("run this from the project root: Rscript api/run.R", call. = FALSE)
}

# plumber chdirs into api/ while sourcing plumber.R; hand it the real root.
# Respect an existing value so the container can set it in the image.
if (!nzchar(Sys.getenv("LIPID_PROJECT_ROOT"))) {
  Sys.setenv(LIPID_PROJECT_ROOT = normalizePath(getwd()))
}

host <- Sys.getenv("LIPID_API_HOST", "127.0.0.1")
port <- as.integer(Sys.getenv("LIPID_API_PORT", "7413"))

message("lipid API on http://", host, ":", port,
        "   docs at http://", host, ":", port, "/__docs__/")

plumber::pr("api/plumber.R") |>
  plumber::pr_set_docs("swagger") |>
  plumber::pr_run(host = host, port = port)
