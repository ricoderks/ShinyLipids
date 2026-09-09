## PROTOTYPE -- not wired into the app yet.
##
## Parses MS-DIAL lipid shorthand into chain-level columns, so that "any lipid
## containing 20:4" can be a structural query instead of a substring match.
## Substring search is wrong in both directions: 34% of names containing the
## text "20:4" have no 20:4 acyl chain (they are ether O-20:4 or oxidised
## 20:4;O), and it misses the same species written at sum-composition level
## (PE 38:4 vs PE 18:0_20:4).
##
## Validate with reconcile_masses() before trusting any of this: a parse that
## silently drops a chain is worse than no parse at all.

# Monoisotopic masses used to reconstruct a mass from a chain composition.
M_CH2 <- 14.0156500642   # each additional carbon, as CH2
M_H2  <-  2.0156500642   # each double bond removes H2
M_O   <- 15.9949146221   # each ";O" oxygen

# One chain token: an optional ether/alkyl link, carbons:double-bonds, and an
# optional ";O" oxygen count. Deuterium tags like "(d7)" carry no colon and so
# are skipped by construction.
CHAIN_RE <- "(?:[OP]-)?[0-9]+:[0-9]+(?:;O[0-9]*)?"

#' Parse a vector of lipid names into chain-level features.
#'
#' @return data.frame with one row per input name:
#'   n_chains, total_c, total_db, total_o, n_ether, sn_resolved,
#'   sum_composition, labelled, chains (the tokens, "|"-joined)
parse_lipid_names <- function(names) {
  # Internal standards carry an alias: "PC 34:1(d7)|PC 16:0_18:1(d7)". The
  # right-hand side is the more informative molecular-species form.
  nm <- sub("^.*\\|", "", names)
  body <- sub("^[^ ]+ ", "", nm)          # drop the class prefix

  toks <- regmatches(body, gregexpr(CHAIN_RE, body, perl = TRUE))
  lens <- lengths(toks)
  flat <- unlist(toks, use.names = FALSE)
  idx  <- rep.int(seq_along(toks), lens)

  ether <- grepl("^[OP]-", flat)
  core  <- sub("^[OP]-", "", flat)
  cs    <- as.integer(sub(":.*$", "", core))
  rest  <- sub("^[0-9]+:", "", core)
  dbs   <- as.integer(sub(";.*$", "", rest))

  has_ox <- grepl(";O", core, fixed = TRUE)
  ox <- integer(length(core))
  if (any(has_ox)) {
    digits <- sub("^.*;O", "", core[has_ox])
    n <- suppressWarnings(as.integer(digits))
    n[is.na(n)] <- 1L                     # bare ";O" means one oxygen
    ox[has_ox] <- n
  }

  # A plasmenyl ("P-") vinyl-ether chain carries one more degree of
  # unsaturation than the plasmanyl ("O-") alkyl form: MS-DIAL's own
  # convention makes PE P-18:0 equivalent to PE O-18:1. Counting it keeps
  # total_db chemically meaningful and reconciles the EtherPE class.
  plasmenyl <- grepl("^P-", flat)
  dbs <- dbs + as.integer(plasmenyl)

  # rowsum() over the flattened tokens is far faster than per-name apply().
  agg <- function(v) {
    out <- numeric(length(toks))
    if (!length(v)) return(out)
    r <- rowsum(v, idx, reorder = TRUE)
    out[as.integer(rownames(r))] <- r[, 1]
    out
  }

  # Branch methyls are written outside the chain token -- "18:0(methyl)" is a
  # 19-carbon branched chain -- so they must be added back explicitly.
  n_methyl <- lengths(gregexpr("(methyl)", body, fixed = TRUE))
  n_methyl[!grepl("(methyl)", body, fixed = TRUE)] <- 0L

  # Positional hydroxyls, written as ";2OH,3OH" or "(2OH)", each contribute an
  # oxygen that the ";O<n>" token does not cover.
  has_oh <- grepl("[0-9]+OH", body)
  n_oh <- integer(length(body))
  n_oh[has_oh] <- lengths(gregexpr("[0-9]+OH", body[has_oh], perl = TRUE))

  # The "(FA 12:0)" shorthand names an esterified species: the same molecule
  # that the "/" form writes out in full (HexCer 24:2;O3(FA 12:0) and
  # HexCer 12:2;O2/24:1;O2 are both C42H77NO10). The ester is a condensation,
  # so the two notations encode different double-bond and oxygen counts for
  # one compound -- reconciliation has to compare like with like.
  fa_notation <- grepl("(FA ", body, fixed = TRUE)

  data.frame(
    name            = names,
    n_chains        = lens,
    n_methyl        = n_methyl,
    n_plasmenyl     = agg(as.numeric(plasmenyl)),
    total_c         = agg(cs) + n_methyl,
    total_db        = agg(dbs),
    total_o         = agg(ox) + n_oh,
    fa_notation     = fa_notation,
    n_ether         = agg(as.numeric(ether)),
    sn_resolved     = grepl("/", body, fixed = TRUE),
    sum_composition = lens == 1L & !grepl("[_/]", body),
    labelled        = grepl("\\(d[0-9]+\\)|_d[0-9]+", names),
    chains          = vapply(toks, paste, "", collapse = "|"),
    stringsAsFactors = FALSE
  )
}

#' Check a parse against the library's own exact masses.
#'
#' Within one lipid class the neutral mass is an exact linear function of the
#' chain composition: base + C*CH2 - DB*H2 + O*oxygen. The base is whatever the
#' head group weighs, so it need not be known in advance -- it is read off the
#' data as the modal residual. A name whose deviation from that base is not
#' ~0 was parsed wrongly (a dropped chain shifts it by whole CH2 units).
#'
#' This is the point of the exercise: it catches silent misparses without
#' needing per-class head-group chemistry.
reconcile_masses <- function(df, tol = 0.005) {
  resid <- df$exact_mass - (df$total_c * M_CH2 - df$total_db * M_H2 +
                            df$total_o * M_O)
  # Modal residual per group, to 3 dp -- robust to a minority of bad parses.
  #
  # The group is class *and* chain count, not class alone. A class can hold
  # structurally different forms whose head groups differ: in Cer_EOS the
  # 3-chain esterified species is a condensation (one H2O lighter) of the
  # 2-chain free form, so pooling them puts a constant offset on one subgroup
  # and makes a correct parse look wrong.
  key  <- paste(df$lipid_class, df$n_chains, df$fa_notation)
  base <- tapply(round(resid, 3), key, function(x) {
    t <- table(x); as.numeric(names(t)[which.max(t)])
  })
  df$residual <- resid - base[key]
  df$ok <- abs(df$residual) < tol
  df
}
