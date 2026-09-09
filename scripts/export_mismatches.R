## Export library records whose name and stored formula disagree.
##
## The lipid shorthand in `name` implies a chain composition, and therefore a
## neutral mass. Where that mass cannot be reconciled with the stored
## `exact_mass`/`formula`, the two disagree -- the entry contradicts itself.
##
## See R/lipid_names.R for the parser and the reconciliation method. Records in
## classes with no chain shorthand at all (CoQ, Vitamin_E, Unknown, ...) are
## excluded: there is nothing to contradict, so they are not mismatches.
##
## Usage:  Rscript scripts/export_mismatches.R [outfile.csv]

source("R/db.R")
source("R/lipid_names.R")

out <- commandArgs(trailingOnly = TRUE)
out <- if (length(out)) out[1] else "name_formula_mismatches.csv"

con <- lipid_db()
message("reading library ...")
d <- dbGetQuery(con, "SELECT id, name, lipid_class, adduct, ion_mode,
                             precursor_mz, formula, exact_mass FROM lipid")

message("parsing ", format(nrow(d), big.mark = ","), " names ...")
p <- parse_lipid_names(d$name)
p$lipid_class <- d$lipid_class
p$exact_mass  <- d$exact_mass
r <- reconcile_masses(p)

bad <- !r$ok & r$n_chains > 0

# Deuterated standards are NOT mismatches: the formula spells the label out
# ("C45H71[2H]7O2") and is perfectly consistent with the name -- it is the
# reconciliation model that has no deuterium term. Excluding them keeps
# correct entries out of an error report.
is_labelled <- grepl("[2H]", d$formula, fixed = TRUE)
bad <- bad & !is_labelled

message(format(sum(bad), big.mark = ","),
        " records disagree (excluded ", format(sum(!r$ok & is_labelled), big.mark = ","),
        " deuterated standards)")

# Express each mass gap in whole atomic units, so the discrepancy is readable
# rather than a bare number of daltons. Search the smallest combination of
# oxygen / CH2 / H2 that accounts for the residual.
explain <- local({
  grid <- expand.grid(O = -3:3, CH2 = -3:3, H2 = -3:3)
  grid$mass <- grid$O * M_O + grid$CH2 * M_CH2 + grid$H2 * M_H2
  grid$cost <- abs(grid$O) + abs(grid$CH2) + abs(grid$H2)
  grid <- grid[order(grid$cost), ]
  function(res) {
    vapply(res, function(x) {
      hit <- which(abs(grid$mass - x) < 0.01)
      if (!length(hit)) return("unclassified")
      g <- grid[hit[1], ]
      if (g$cost == 0) return("none")
      part <- c(
        if (g$O)   sprintf("%+d O", g$O),
        if (g$CH2) sprintf("%+d CH2", g$CH2),
        if (g$H2)  sprintf("%+d H2", g$H2)
      )
      paste("formula vs name:", paste(part, collapse = " "))
    }, "")
  }
})

e <- data.frame(
  id              = d$id[bad],
  name            = d$name[bad],
  lipid_class     = d$lipid_class[bad],
  adduct          = d$adduct[bad],
  ion_mode        = d$ion_mode[bad],
  precursor_mz    = round(d$precursor_mz[bad], 4),
  formula         = d$formula[bad],
  exact_mass      = round(d$exact_mass[bad], 5),
  parsed_chains   = r$chains[bad],
  parsed_n_chains = r$n_chains[bad],
  parsed_total_c  = r$total_c[bad],
  parsed_total_db = r$total_db[bad],
  parsed_total_o  = r$total_o[bad],
  deuterated      = r$labelled[bad],
  mass_gap_da     = round(r$residual[bad], 5),
  discrepancy     = explain(r$residual[bad]),
  stringsAsFactors = FALSE
)
e <- e[order(e$lipid_class, e$name, e$adduct), ]

# Only a gap that resolves to whole atoms is reported as a contradiction. An
# "unclassified" gap means the model could not account for it -- e.g. the
# bile-acid conjugate notation ";G;"/";T;", which adds glycine (57.021) or
# taurine (107.004) that the parser does not model. Those are gaps in this
# script, not proven library errors, so they go to a separate file rather than
# into an error report.
unexplained <- e[e$discrepancy == "unclassified", ]
e <- e[e$discrepancy != "unclassified", ]

write.csv(e, out, row.names = FALSE)
message("wrote ", out, " (", format(nrow(e), big.mark = ","), " rows)")

if (nrow(unexplained)) {
  out2 <- sub("\\.csv$", "_unexplained.csv", out)
  write.csv(unexplained, out2, row.names = FALSE)
  message("wrote ", out2, " (", format(nrow(unexplained), big.mark = ","),
          " rows whose mass gap this script cannot attribute)")
}

cat("\n=== by class ===\n")
print(as.data.frame(table(class = e$lipid_class), stringsAsFactors = FALSE),
      row.names = FALSE)

cat("\n=== by discrepancy ===\n")
print(as.data.frame(table(discrepancy = e$discrepancy), stringsAsFactors = FALSE),
      row.names = FALSE)

cat("\n=== distinct names affected ===\n")
u <- unique(e[, c("name", "formula", "discrepancy")])
cat(format(nrow(u), big.mark = ","), "distinct name/formula pairs\n")

DBI::dbDisconnect(con)
