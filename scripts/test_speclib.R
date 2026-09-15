## Checks for the whole-library spectral search (R/speclib.R).
##
## Usage:  Rscript scripts/test_speclib.R   (from the project root)
##
## Runs against the real DB/lipids.sqlite, so it builds the full peak index
## first (~10 s, ~1 GB). Record 3620661 is PC 18:0_20:4 [M+HCOO]-, a nine-peak
## negative-mode spectrum that is a good stand-in for a real query.

suppressMessages({
  source("R/db.R")
  source("R/speclib.R")
})

REF_ID <- 3620661L

fails <- 0L
ok <- function(label, cond) {
  pass <- isTRUE(cond)
  if (!pass) fails <<- fails + 1L
  cat(sprintf("%-50s %s\n", label, if (pass) "PASS" else "** FAIL **"))
}

# ----------------------------------------------------------- parsing -------

p <- parse_peak_text("283.2643,100\n303.2330,55.3")
ok("comma-separated pairs parse", nrow(p$peaks) == 2L)
p <- parse_peak_text("m/z\tintensity\n283.2643\t1000\n303.2330\t553")
ok("header line is skipped", nrow(p$peaks) == 2L && p$n_skipped == 1L)
ok("intensities normalise to the base peak", abs(p$peaks$rel[1] - 100) < 1e-9)
p <- parse_peak_text("NAME: PE 18:0_20:4\nNum Peaks: 2\n283.2643 999 \"FA 18:0\"\n303.2330 550")
ok("MSP metadata and a third column are ignored", nrow(p$peaks) == 2L)
ok("semicolons work", nrow(parse_peak_text("283.2643; 10\n303.2330 ;5")$peaks) == 2L)
ok("blank input gives no peaks", nrow(parse_peak_text("  ")$peaks) == 0L)
ok("non-positive intensities are dropped",
   nrow(parse_peak_text("300,0\n301,-5\n302,10")$peaks) == 1L)

ok("Da tolerance alone", mz_tolerance(800, 0.01, 0) == 0.01)
ok("ppm tolerance alone", abs(mz_tolerance(800, 0, 10) - 0.008) < 1e-12)
ok("the wider of the two wins",
   identical(mz_tolerance(c(100, 800), 0.005, 10), c(0.005, 0.008)))

# ------------------------------------------------------------- MSP ---------

msp <- paste(
  "NAME: PC 34:1", "PRECURSORMZ: 760.5851", "PRECURSORTYPE: [M+H]+",
  "IONMODE: Positive", "RETENTIONTIME: 9.91", "Num Peaks: 3",
  "184.0733\t999", "478.3292\t120", "760.5851\t300", "",
  "NAME: PE 36:2", "PRECURSORMZ: 742.5387", "PRECURSORTYPE: [M-H]-",
  "IONMODE: Negative", "Num Peaks: 2",
  "281.2486\t999", "283.2643\t850", "",
  sep = "\n")
m <- parse_msp(msp)
ok("two records parse from one file", nrow(m$meta) == 2L)
ok("headers are read", m$meta$precursor_mz[1] == 760.5851 &&
     m$meta$adduct[2] == "[M-H]-" && m$meta$ion_mode[2] == "Negative")
ok("peaks land with the right record",
   nrow(m$peaks[[1]]) == 3L && nrow(m$peaks[[2]]) == 2L)
ok("peaks normalise per record", abs(max(m$peaks[[2]]$rel) - 100) < 1e-9)
ok("retention time is optional", is.na(m$meta$retention_time[2]))

alias <- parse_msp(paste("Name: A", "Precursor_type: [M+Na]+",
                         "PRECURSOR_MZ: 100.5", "Polarity: POSITIVE",
                         "Num Peaks: 1", "50.1 100", sep = "\n"))
ok("field-name dialects are accepted",
   alias$meta$adduct == "[M+Na]+" && alias$meta$precursor_mz == 100.5 &&
     alias$meta$ion_mode == "Positive")

noblank <- parse_msp(paste("NAME: A", "Num Peaks: 2",
                           "50.1\t100\t\"frag: x\"", "60.2\t50",
                           "NAME: B", "Num Peaks: 1", "70.0\t10", sep = "\n"))
ok("records split on NAME without blank lines", nrow(noblank$meta) == 2L)
ok("an annotation in a third peak column is ignored",
   nrow(noblank$peaks[[1]]) == 2L)

noname <- parse_msp("50.1 100\n60.2 50\n\n70.0 10\n80.0 5\n")
ok("a file with no NAME splits on blank lines", nrow(noname$meta) == 2L)
ok("unnamed records get a placeholder name", noname$meta$name[1] == "spectrum 1")

ok("an unreadable polarity stays empty",
   parse_msp("NAME: A\nIONMODE: N/A\nNum Peaks: 1\n50.1 100")$meta$ion_mode == "")
drop <- parse_msp("NAME: header only\nPRECURSORMZ: 100\n\nNAME: real\nNum Peaks: 1\n50.1 100")
ok("a record with no peaks is dropped",
   nrow(drop$meta) == 1L && drop$n_dropped == 1L)
ok("an empty file yields no records", nrow(parse_msp("")$meta) == 0L)
ok("CRLF line endings work",
   nrow(parse_msp("NAME: A\r\nNum Peaks: 1\r\n50.1\t100\r\n")$meta) == 1L)

# The sample the app offers for download has to be one the app can read back.
if (file.exists("example.msp")) {
  ex <- parse_msp(readLines("example.msp", warn = FALSE))
  ok("the shipped example.msp parses",
     nrow(ex$meta) > 0 && all(ex$meta$n_peaks > 0))
  ok("the shipped example.msp declares a precursor and polarity per record",
     all(is.finite(ex$meta$precursor_mz)) && all(nzchar(ex$meta$ion_mode)))
} else {
  cat("example.msp not present -- skipping the sample-file checks\n")
}

# ------------------------------------------------------------ search -------

con <- lipid_db()
cat("\nbuilding the peak index (one-off, ~10 s)...\n")
ix <- peak_index(con)
cat(sprintf("indexed %s peaks over %s records\n\n",
            format(ix$n_peaks_total, big.mark = ","),
            format(ix$n_rec, big.mark = ",")))

q <- get_spectrum(con, REF_ID)$peaks

r <- search_spectrum(ix, q, top_n = 5)
ok("an exact self-match ranks first", r$id[1] == REF_ID)
ok("an exact self-match scores 1 on all three",
   all(abs(c(r$dot[1], r$wdot[1], r$rdot[1]) - 1) < 1e-9))
ok("reverse dot is never below dot", all(r$rdot >= r$dot - 1e-12))
ok("every score lies in [0, 1]",
   all(c(r$dot, r$wdot, r$rdot) >= -1e-12 &
       c(r$dot, r$wdot, r$rdot) <= 1 + 1e-9))

# A realistic experimental spectrum: every peak 6 ppm high, intensities noisy.
set.seed(42)
qp <- data.frame(mz = q$mz * (1 + 6e-6),
                 rel = pmin(100, pmax(1, q$rel * runif(nrow(q), 0.9, 1.1))))
r_pp <- search_spectrum(ix, qp, tol_da = 0, tol_ppm = 10, top_n = 5)
r_da <- search_spectrum(ix, qp, tol_da = 0.001, tol_ppm = 0, top_n = 5)
ok("a 10 ppm window catches a 6 ppm drift",
   r_pp$id[1] == REF_ID && r_pp$n_match[1] == nrow(q))
ok("a 0.001 Da window does not", r_da$n_match[1] < nrow(q))

ok("the polarity filter holds",
   all(search_spectrum(ix, q, mode = "Positive", top_n = 5)$ion_mode == "Positive"))
ok("the adduct filter holds",
   all(search_spectrum(ix, q, adduct = "[M-H]-", top_n = 5)$adduct == "[M-H]-"))
r_pre <- search_spectrum(ix, q, precursor = 854.5917, tol_da = 0.01, top_n = 500)
ok("the precursor filter holds",
   all(abs(r_pre$precursor_mz - 854.5917) <= 0.01))
ok("the filters narrow the candidate count",
   attr(r_pre, "n_filtered") < attr(r, "n_filtered"))
ok("min_match is respected",
   all(search_spectrum(ix, q, min_match = 5, top_n = 500)$n_match >= 5))

ok("an empty query returns nothing", nrow(search_spectrum(ix, q[0, ])) == 0L)
ok("an m/z beyond the library returns nothing",
   nrow(search_spectrum(ix, data.frame(mz = 12345.6, rel = 100))) == 0L)
ok("contradictory filters return nothing",
   nrow(search_spectrum(ix, q, mode = "Positive", adduct = "[M-H]-")) == 0L)

for (k in c("dot", "wdot", "rdot")) {
  rr <- search_spectrum(ix, q, rank_by = k, top_n = 50)
  ok(sprintf("rank_by = '%s' orders the result", k), !is.unsorted(rev(rr[[k]])))
}

r_mb <- search_spectrum(ix, q, mz_power = 2, int_power = 0.5, top_n = 5)
ok("MassBank exponents still self-match exactly",
   r_mb$id[1] == REF_ID && abs(r_mb$wdot[1] - 1) < 1e-9)
ok("weighted norms are memoised per exponent pair", length(ls(ix$wnorm)) == 2L)

# The mirror plot rescores the selected pair on its own; if the two ever
# disagree the table and the number under the plot would contradict.
# A whole MSP file, each record searched under its own header. Every spectrum
# here is a library record put back through the search with a realistic ppm
# drift and intensity noise, so each must find its own source at rank 1.
set.seed(11)
batch_ids <- dbGetQuery(con, "SELECT id FROM lipid
                               WHERE n_peaks BETWEEN 6 AND 14
                                 AND adduct IN ('[M+HCOO]-','[M+H]+','[M+NH4]+','[M-H]-')
                               ORDER BY id LIMIT 12 OFFSET 40000")$id
msp_lines <- unlist(lapply(batch_ids, function(id) {
  sp <- get_spectrum(con, id); mm <- sp$meta
  pk <- sp$peaks
  pk$mz <- pk$mz * (1 + runif(nrow(pk), 2e-6, 6e-6))
  pk$rel <- pmax(0.5, pk$rel * runif(nrow(pk), 0.85, 1.15))
  c(sprintf("NAME: %s", mm$name),
    sprintf("PRECURSORMZ: %.4f", mm$precursor_mz * (1 + 4e-6)),
    sprintf("PRECURSORTYPE: %s", mm$adduct),
    sprintf("IONMODE: %s", mm$ion_mode),
    sprintf("Num Peaks: %d", nrow(pk)),
    sprintf("%.4f\t%.1f", pk$mz, pk$rel), "")
}))
batch <- parse_msp(msp_lines)
ok("every record survives the MSP round trip",
   nrow(batch$meta) == length(batch_ids))

t0 <- Sys.time()
batch_hits <- lapply(seq_len(nrow(batch$meta)), function(i) {
  mm <- batch$meta[i, ]
  search_spectrum(ix, batch$peaks[[i]], tol_da = 0, tol_ppm = 10,
                  mode = mm$ion_mode, adduct = mm$adduct,
                  precursor = mm$precursor_mz, min_match = 2L, top_n = 25)
})
per <- 1000 * as.numeric(Sys.time() - t0, units = "secs") / nrow(batch$meta)
ok("each MSP record finds its own source record first",
   all(vapply(seq_along(batch_ids), function(i)
     nrow(batch_hits[[i]]) > 0 && batch_hits[[i]]$id[1] == batch_ids[i], TRUE)))
ok(sprintf("a batch search stays interactive (%.0f ms per spectrum)", per),
   per < 200)

r20 <- search_spectrum(ix, q, top_n = 20)
delta <- vapply(seq_len(nrow(r20)), function(i) {
  sp <- score_pair(q, get_spectrum(con, r20$id[i])$peaks)
  max(abs(c(sp$dot - r20$dot[i], sp$wdot - r20$wdot[i], sp$rdot - r20$rdot[i],
            length(sp$match$a) - r20$n_match[i])))
}, numeric(1))
ok("score_pair() agrees with the bulk scorer over 20 hits", max(delta) < 1e-9)

cat(sprintf("\n%s\n", if (fails == 0L) "all checks passed"
                      else sprintf("%d check(s) failed", fails)))
quit(status = if (fails == 0L) 0L else 1L)
