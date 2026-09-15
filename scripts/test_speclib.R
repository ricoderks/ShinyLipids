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
