# Read-only REST API over the MS-DIAL lipid library.
#
# The image carries code only. The database (2.4 GB) and the flagged-record
# export are mounted at /data at runtime, so the image stays small, rebuilds
# are fast, and the same image serves any version of the library.
#
#   docker build -t shinylipids-api .
#   docker run --rm -p 7413:7413 \
#     -v "$PWD/DB/lipids.sqlite:/data/lipids.sqlite:ro" \
#     -v "$PWD/name_formula_mismatches.csv:/data/name_formula_mismatches.csv:ro" \
#     shinylipids-api
#
# Docs then at http://localhost:7413/__docs__/

FROM rocker/r-ver:4.6.1

# libsodium is for the `sodium` package, a hard dependency of plumber (it is
# used for encrypted session cookies, which this read-only API never sets).
# curl is used by the HEALTHCHECK below.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libsodium-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        curl \
    && rm -rf /var/lib/apt/lists/*

# rocker/r-ver pins a Posit Package Manager snapshot matching the R version, so
# these resolve to binaries and to the same versions on every rebuild.
RUN install2.r --error --skipinstalled --ncpus -1 \
        plumber \
        DBI \
        RSQLite \
    && rm -rf /tmp/downloaded_packages

WORKDIR /app

# Only the code the API actually loads. R/db.R and R/spectrum.R are shared with
# the Shiny app; nothing Shiny-specific is needed here.
COPY R/db.R R/spectrum.R ./R/
COPY api/ ./api/

ENV LIPID_PROJECT_ROOT=/app \
    LIPID_DB=/data/lipids.sqlite \
    LIPID_FLAGGED=/data/name_formula_mismatches.csv \
    LIPID_API_HOST=0.0.0.0 \
    LIPID_API_PORT=7413

# Binding to 0.0.0.0 is required inside a container -- it is the container's
# own interface, not the host's. Publish the port deliberately (-p 127.0.0.1:
# 7413:7413 to keep it local), and put a proxy in front before exposing it:
# the API has no authentication and no rate limiting.
EXPOSE 7413

# Run unprivileged. The API only ever reads, and /data should be mounted :ro.
RUN useradd --uid 10001 --create-home --shell /usr/sbin/nologin apiuser
USER 10001

HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
    CMD curl -fsS http://127.0.0.1:7413/health || exit 1

CMD ["Rscript", "api/run.R"]
