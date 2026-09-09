## Structural interpretation of MS/MS peaks.
##
## Answers "what is this fragment?" for a peak in a library spectrum, so the
## hover text can say "[M+H-RCOOH]+ - neutral loss of NH3 and 16:0 as free
## acid" instead of only repeating the m/z the user is already pointing at.
##
## The approach is to predict rather than to look up: the library stores no
## fragment annotations, so candidate ions are computed from the lipid's own
## name (which carries the chain composition) plus its adduct, and then matched
## back against the observed peaks. A peak with no candidate within tolerance
## is simply left unannotated -- silence is the honest answer for the classes
## whose chemistry is not encoded here, and it is far better than guessing.

# ------------------------------------------------------------- masses -------
# Monoisotopic, CODATA/IUPAC values. The electron mass matters at the 4th
# decimal, which is exactly the precision the peak labels show.
E_H <- 1.00782503207
E_C <- 12
E_N <- 14.0030740048
E_O <- 15.9949146196
E_P <- 30.97376163
E_D <- 2.01410177811                     # deuterium, 2H
E_E <- 0.00054857991

# What one deuterium substitution adds. The labelled standards in the library
# are labelled by substitution, not by addition, so every rule that touches a
# labelled group shifts by a multiple of this rather than by a whole hydrogen.
D_SHIFT <- E_D - E_H                     # 1.006277

PROTON   <- E_H - E_E                    # 1.007276
M_H2     <- 2 * E_H
M_H2O    <- 2 * E_H + E_O                # 18.010565
M_NH3    <- E_N + 3 * E_H                # 17.026549
M_CH2    <- E_C + 2 * E_H                # 14.015650
M_CH3CAT <- E_C + 3 * E_H - E_E          # 15.022928, lost as a methyl cation
M_CD3CAT <- E_C + 3 * E_D - E_E          # 18.041757, the same from a d9 choline
M_CH2O   <- E_C + 2 * E_H + E_O          # 30.010565, formaldehyde
M_CH4O   <- E_C + 4 * E_H + E_O          # 32.026215, methanol
M_C2H2   <- 2 * E_C + 2 * E_H            # 26.015650, the C1-C2 of a sphingoid base
M_C3H4O  <- 3 * E_C + 4 * E_H + E_O      # 56.026215, C1-C3 of a phyto base
M_C2H5N  <- 2 * E_C + 5 * E_H + E_N      # 43.042199, aziridine
M_NME3   <- 3 * E_C + 9 * E_H + E_N      # 59.073499, trimethylamine
M_HCOOH  <- 2 * E_H + E_C + 2 * E_O      # 46.005479
M_CH3COOH<- 2 * E_C + 4 * E_H + 2 * E_O  # 60.021129
M_H2CO3  <- 2 * E_H + E_C + 3 * E_O      # 62.000394
M_NA     <- 22.9897692809

M_GLYCEROL <- 3 * E_C + 8 * E_H + 3 * E_O   #  92.047344
M_GPA      <- 3 * E_C + 9 * E_H + 6 * E_O + E_P  # 172.013674, glycerophosphate
M_H3PO4    <- 3 * E_H + E_P + 4 * E_O       #  97.976896
M_HEXOSE   <- 6 * E_C + 12 * E_H + 6 * E_O  # 180.063388, free hexose
M_HEX_RES  <- M_HEXOSE - M_H2O              # 162.052824, hexose residue
M_PCHOLINE <- 5 * E_C + 14 * E_H + E_N + 4 * E_O + E_P  # 183.066045
M_PETN     <- 2 * E_C + 8 * E_H + E_N + 4 * E_O + E_P   # 141.019094, phosphoethanolamine

# ------------------------------------------------------------- adducts ------
# n = molecules per ion, z = signed charge, d = mass added to n*M.
# m/z = (n * M + d) / |z|
ADDUCT_DEF <- list(
  "[M+H]+"      = c(n = 1, z =  1, d =  PROTON),
  "[M+NH4]+"    = c(n = 1, z =  1, d =  M_NH3 + PROTON),
  "[M+Na]+"     = c(n = 1, z =  1, d =  M_NA - E_E),
  "[M+H-H2O]+"  = c(n = 1, z =  1, d =  PROTON - M_H2O),
  "[M+2H]2+"    = c(n = 1, z =  2, d =  2 * PROTON),
  "[M+2NH4]2+"  = c(n = 1, z =  2, d =  2 * (M_NH3 + PROTON)),
  "[2M+H]+"     = c(n = 2, z =  1, d =  PROTON),
  "[M-H]-"      = c(n = 1, z = -1, d = -PROTON),
  "[M-2H]2-"    = c(n = 1, z = -2, d = -2 * PROTON),
  "[2M-H]-"     = c(n = 2, z = -1, d = -PROTON),
  "[M+HCOO]-"   = c(n = 1, z = -1, d =  M_HCOOH - PROTON),
  "[M+CH3COO]-" = c(n = 1, z = -1, d =  M_CH3COOH - PROTON),
  "[M+HCO3]-"   = c(n = 1, z = -1, d =  M_H2CO3 - PROTON)
)

# What the adduct has to shed to reach the protonated / deprotonated molecule.
# Used only for wording, so that a TG peak reads "loss of NH3" rather than the
# uninformative "[M+H]+".
ADDUCT_SHED <- c(
  "[M+NH4]+"    = "NH3",
  "[M+2NH4]2+"  = "NH3",
  "[M+HCOO]-"   = "HCOOH",
  "[M+CH3COO]-" = "CH3COOH",
  "[M+HCO3]-"   = "H2CO3"
)

#' Neutral monoisotopic mass behind a precursor ion.
#'
#' Prefers the library's own exact_mass and falls back to inverting the adduct,
#' so records predating that column still annotate.
neutral_mass <- function(precursor_mz, adduct, exact_mass = NA_real_) {
  if (!is.null(exact_mass) && length(exact_mass) == 1L &&
      !is.na(exact_mass) && exact_mass > 0) return(as.numeric(exact_mass))
  a <- ADDUCT_DEF[[adduct %or% ""]]
  if (is.null(a)) return(NA_real_)
  (precursor_mz * abs(a[["z"]]) - a[["d"]]) / a[["n"]]
}

# Not `%||%`: base R defines that one, and this is deliberately stricter --
# it treats NA and a zero-length value as missing too, which is what an absent
# adduct or exact mass looks like coming out of the database.
`%or%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x[1])) y else x

# ------------------------------------------------------- name parsing -------

# One chain token: optional ether/vinyl-ether link, carbons:double-bonds, and
# an optional ";O<n>" oxygen count. Deuterium tags like "(d7)" carry no colon
# and are skipped by construction.
FRAG_CHAIN_RE <- "(?:[OP]-)?[0-9]+:[0-9]+(?:;O[0-9]*)?"

# Classes built on a sphingoid base rather than glycerol. The first chain of
# such a name is the long-chain base, not an acyl chain, and it fragments
# completely differently -- annotating it as a fatty acid would be wrong.
SPHINGO_RE <- "Cer|^SM|^ASM|^SL$|^SPB|^MIPC$"

#' Split a lipid shorthand name into its chains.
#'
#' @return data.frame(label, c, db, ox, kind, on_hg) with kind one of "acyl",
#'   "ether" or "lcb", and on_hg marking a chain esterified to the head group
#'   rather than to the backbone. A zero-row frame when nothing parses.
parse_chains <- function(name, lipid_class = "") {
  if (is.null(name) || is.na(name) || !nzchar(name)) return(empty_chains())
  # Internal standards carry an alias, "PC 34:1(d7)|PC 16:0_18:1(d7)"; the
  # right-hand side is the more informative molecular-species form.
  nm <- sub("^.*\\|", "", name)

  m <- gregexpr(FRAG_CHAIN_RE, nm, perl = TRUE)
  toks <- regmatches(nm, m)[[1]]
  if (!length(toks)) return(empty_chains())
  starts <- as.integer(m[[1]])
  ends   <- starts + attr(m[[1]], "match.length") - 1L

  core <- sub("^[OP]-", "", toks)
  cs   <- as.integer(sub(":.*$", "", core))
  rest <- sub("^[0-9]+:", "", core)
  dbs  <- as.integer(sub(";.*$", "", rest))

  ox <- integer(length(core))
  has_ox <- grepl(";O", core, fixed = TRUE)
  if (any(has_ox)) {
    n <- suppressWarnings(as.integer(sub("^.*;O", "", core[has_ox])))
    n[is.na(n)] <- 1L                       # a bare ";O" means one oxygen
    ox[has_ox] <- n
  }

  # MS-DIAL writes a plasmenyl (vinyl-ether) chain as "P-18:0" for what is
  # compositionally O-18:1, so the extra unsaturation has to be added back.
  dbs <- dbs + as.integer(grepl("^P-", toks))

  # A positional hydroxyl is written outside the chain token -- the alpha-
  # hydroxy acyl of Cer_AS and Cer_AP is "11:0(2OH)", not "11:0;O". Dropping it
  # put every fragment of those classes 16 Da out, so each "<n>OH" is charged
  # to the chain token it follows.
  gaps <- substring(nm, ends + 1L, c(utils::tail(starts, -1L), nchar(nm) + 1L) - 1L)
  ox <- ox + vapply(gregexpr("[0-9]+OH", gaps, perl = TRUE),
                    function(g) sum(g > 0), 0L)

  kind <- ifelse(grepl("^[OP]-", toks), "ether", "acyl")

  # AHexCer, ASHexCer and ADGGA write the fatty acid esterified to the head
  # group in parentheses -- "AHexCer (O-12:0)12:2;O2/12:1". Despite the "O-"
  # that is an ester, not an ether, and it fragments like any other acyl
  # chain. The parentheses are the only signal, and they are used for nothing
  # else: bare "O-"/"P-" tokens are always genuine ether links.
  on_hg <- starts > 1L & substring(nm, starts - 1L, starts - 1L) == "("
  kind[on_hg] <- "acyl"

  # The long-chain base is the token written immediately before the "/" that
  # separates base from N-acyl chain. Without a "/" the name is at sum-
  # composition level and no chain can be assigned a role.
  if (grepl(SPHINGO_RE, lipid_class) && grepl("/", nm, fixed = TRUE)) {
    slash <- gregexpr("/", nm, fixed = TRUE)[[1]]
    lcb <- which((ends + 1L) %in% slash)
    if (length(lcb)) kind[lcb[1]] <- "lcb"
  }

  data.frame(label = toks, c = cs, db = dbs, ox = ox, kind = kind,
             on_hg = on_hg, stringsAsFactors = FALSE)
}

empty_chains <- function() {
  data.frame(label = character(0), c = integer(0), db = integer(0),
             ox = integer(0), kind = character(0), on_hg = logical(0),
             stringsAsFactors = FALSE)
}

# Free fatty acid RCOOH: C(c) H(2c-2db) O(2+ox).
mass_fa <- function(c, db, ox) c * E_C + (2 * c - 2 * db) * E_H + (2 + ox) * E_O
# The alkyl/alkenyl alcohol behind an ether chain: C(c) H(2c+2-2db) O(1+ox).
mass_alkanol <- function(c, db, ox) c * E_C + (2*c + 2 - 2*db) * E_H + (1 + ox) * E_O
# Sphingoid base: C(c) H(2c+3-2db) N O(ox). d18:1;O2 -> 299.2824.
mass_lcb <- function(c, db, ox) c * E_C + (2*c + 3 - 2*db) * E_H + E_N + ox * E_O

# --------------------------------------------------------- class families ---

#' Collapse 200-odd library classes onto the handful of fragmentation
#' chemistries that actually differ. Anything unmatched returns "" and gets
#' only the class-independent rules.
lipid_family <- function(cls) {
  if (is.null(cls) || is.na(cls)) return("")
  rules <- list(
    c("^(L?PC|EtherL?PC|OxPC|EtherOxPC|L?PC_d5)$",        "PC"),
    c("^SM_d9$",                                           "SM_D9"),
    c("^(SM|ASM)$",                                        "SM"),
    c("^(L?PE|EtherL?PE|OxPE|EtherOxPE|L?PE_d5)$",        "PE"),
    c("^(MMPE|DMPE)$",                                     "MPE"),
    c("^(L?PS|EtherL?PS|OxPS|L?PS_d5)$",                  "PS"),
    c("^(L?PG|EtherL?PG|OxPG|L?PG_d5|BMP|HBMP)$",         "PG"),
    c("^(L?PI|EtherL?PI|OxPI|L?PI_d5)$",                  "PI"),
    c("^(L?PA|EtherL?PA|OxPA)$",                           "PA"),
    c("^(CL|MLCL|DLCL)$",                                  "CL"),
    c("^(M|D|T)G$|^Ether(M|D|T)G$|^OxTG$|^TG_EST$|^TG_d5$|^EtherTG$", "GL"),
    c("^A?SHexCer$",                                       "SHEXCER"),
    c("Hex.*Cer|^Hex[0-9]Cer$",                            "HEXCER"),
    c("Cer",                                               "CER"),
    c("^(WE|CE|SE|ST)$",                                   "WE")
  )
  for (r in rules) if (grepl(r[1], cls, perl = TRUE)) return(r[2])
  ""
}

# Neutral head group lost from the protonated molecule, per family. These are
# the losses that expose the glycerol backbone and so open up the acyl-loss
# series in positive mode.
HG_LOSS_POS <- list(
  PC  = c(mass = M_PCHOLINE,           name = "phosphocholine"),
  SM  = c(mass = M_PCHOLINE,           name = "phosphocholine"),
  # SM_d9 carries its nine deuteriums on the three N-methyls of the choline,
  # so every head-group rule shifts with it while the ceramide half does not.
  SM_D9 = c(mass = M_PCHOLINE + 9 * D_SHIFT, name = "d9-phosphocholine"),
  PE  = c(mass = M_PETN,               name = "phosphoethanolamine"),
  MPE = c(mass = M_PETN + M_CH2,       name = "N-methyl-phosphoethanolamine"),
  PS  = c(mass = M_PETN + 43.98983,    name = "phosphoserine"),
  PG  = c(mass = M_GPA,                name = "phosphoglycerol"),
  PI  = c(mass = 6 * E_C + 13 * E_H + 9 * E_O + E_P,  # C6H13O9P, 260.029722
          name = "inositol phosphate"),
  PA  = c(mass = M_H3PO4,              name = "phosphoric acid"),
  HEXCER = c(mass = M_HEX_RES,         name = "hexose"),
  SHEXCER = c(mass = 2 * E_H + 31.97207069 + 4 * E_O,  # H2SO4, 97.967379
              name = "sulfuric acid")
)

# Diagnostic head-group ions. mz is absolute, not a loss.
HG_IONS_POS <- list(
  PC = list(c(184.073321, "phosphocholine, C5H15NO4P+")),
  SM = list(c(184.073321, "phosphocholine, C5H15NO4P+")),
  SM_D9 = list(c(184.073321 + 9 * D_SHIFT, "d9-phosphocholine, C5H6D9NO4P+"))
)
HG_IONS_NEG <- list(
  PC = list(c(168.042572, "demethylated phosphocholine, C4H11NO4P-"),
            c(224.068786, "demethylated glycerophosphocholine - H2O, C7H15NO5P-")),
  SM = list(c(168.042572, "demethylated phosphocholine, C4H11NO4P-"),
            c( 78.958411, "PO3-")),
  SM_D9 = list(c(168.042572 + 6 * D_SHIFT,
                 "demethylated d9-phosphocholine, C4H5D6NO4P-"),
               c( 78.958411, "PO3-")),
  PE = list(c(140.011817, "phosphoethanolamine - H, C2H7NO4P-"),
            c(196.038082, "glycerophosphoethanolamine - H2O, C5H11NO5P-")),
  PS = list(c( 78.958411, "PO3-"),
            c( 96.969001, "H2PO4-"),
            c(152.995833, "glycerophosphate - H2O, C3H6O5P-")),
  PG = list(c(152.995833, "glycerophosphate - H2O, C3H6O5P-"),
            c(171.006398, "glycerophosphate - H, C3H8O6P-")),
  PA = list(c( 78.958411, "PO3-"),
            c( 96.969001, "H2PO4-"),
            c(152.995833, "glycerophosphate - H2O, C3H6O5P-")),
  PI = list(c(241.011877, "inositol phosphate - H2O, C6H10O8P-"),
            c(259.022442, "inositol phosphate - H, C6H12O9P-"),
            c(297.038092, "glycerophosphoinositol - H2O, C9H14O9P-")),
  CL = list(c(152.995833, "glycerophosphate - H2O, C3H6O5P-"),
            c( 78.958411, "PO3-")),
  HEXCER  = list(c(179.055648, "hexose - H, C6H11O6-"),
                 c(161.045083, "hexose - H2O - H, C6H9O5-")),
  SHEXCER = list(c( 96.960103, "HSO4-"),
                 c( 80.964637, "SO3- radical anion"),
                 c(241.001672, "sulfated hexose - H2O, C6H9O8S-"))
)

# The amine of a phospholipid head group leaves on its own, before the
# phosphate follows: aziridine from an ethanolamine, trimethylamine from a
# choline. Keyed on the head group, so the choline-bearing families stay in
# step with each other. bmPC is deliberately absent: its head group is a CH2
# heavier than phosphocholine, and it is not currently recognised as a family.
HG_AMINE_LOSS <- list(
  PE = c(mass = M_C2H5N, formula = "C2H5N",
         name = "the ethanolamine leaving as aziridine"),
  PC = c(mass = M_NME3,  formula = "C3H9N",
         name = "trimethylamine from the choline"),
  SM = c(mass = M_NME3,  formula = "C3H9N",
         name = "trimethylamine from the choline"),
  SM_D9 = c(mass = M_NME3 + 9 * D_SHIFT, formula = "C3D9N",
            name = "trimethylamine from the labelled choline")
)

# Families whose head group is a choline, and so demethylate and shed
# trimethylamine in negative mode.
CHOLINE_FAMS <- c("PC", "SM", "SM_D9")

# Which backbone ions each sphingolipid subclass actually produces in negative
# mode. Not chemical universals -- this is the library's own in-silico rule
# set, read off the data. Every entry was 100% present or 100% absent across
# thousands of records, measured only on species where the competing
# explanations do not land on the same m/z. That last part matters: an
# [LCB-H-H2O]- and an (N-acyl amide + C1-C2) coincide exactly when the acyl
# has two fewer carbons and one fewer double bond than the base, which is the
# very common 18:1;O2/16:0 case -- and reading the rules off those species is
# how the two ions below came to be labelled as base fragments when they are
# not.
CER_NEG_EMITS <- list(
  # N-acyl amide plus C1-C2 of the base; keeps the C1 hydroxyl in the O2 form.
  amide_c2h2  = "^Cer_N|^HexCer_",
  amide_c2h2o = "^Cer_N",
  amide_nh3   = "^Cer_[AN]D?S(_d[0-9]+)?$",
  # A phytosphingoid base has a third hydroxyl at C4, and cleaving there leaves
  # the amide holding C1-C3 instead of C1-C2. Only the three "P" subclasses do
  # it; Cer_ABP, despite the P, does not, so the pattern is anchored to the
  # exact suffix rather than to "phyto" as a concept.
  amide_c3h4o = "_[AN]P$",
  # Plain water losses off the deprotonated molecule. Too irregular across the
  # subclasses to express as a pattern -- Cer_AS and Cer_ABP do it while Cer_AP
  # does not -- so the members are spelled out. Only Cer_NP goes on to lose a
  # second water.
  m_h2o       = "^(Cer_ABP|Cer_ADS|Cer_AS|Cer_BS|Cer_NP|CerP)$",
  m_2h2o      = "^Cer_NP$",
  # The base stripped of both its C1 and its nitrogen: C1 leaves as methanol
  # and the C2 amine as ammonia, leaving the C3/C4 diol on a chain one carbon
  # shorter (C17H31O2- for an 18:0;O3 base). Cer_NP alone gives it -- the other
  # phytosphingosine classes share the base and do not.
  lcb_ch3oh_nh3 = "^Cer_NP$",
  # The genuine base water-loss ion, which only these two subclasses give.
  lcb_h2o     = "^Cer_AD?S$"
)
cer_emits <- function(rule, cls) {
  !is.null(cls) && !is.na(cls) &&
    grepl(CER_NEG_EMITS[[rule]], as.character(cls), perl = TRUE)
}

# A labelled sphingoid base moves every ion built from it. Cer_NS_d7 carries
# its seven deuteriums on the base rather than on the N-acyl: its base ions all
# sit +7 D up while its amide ions do not move at all.
LCB_LABEL_D <- c(Cer_NS_d7 = 7)
lcb_shift <- function(cls) {
  n <- unname(LCB_LABEL_D[as.character(cls)])
  if (is.na(n)) 0 else n * D_SHIFT
}

# Head-group ions that need the sodium, and so only arise from a sodiated
# precursor. Every PC and LPC record in the library carries two of these three:
# the diacyl species gives C2H5O4PNa+, the lyso species gives the choline
# cation, and both give C5H12N+. Listing all three for the one family costs
# nothing -- the absent one simply never matches a peak -- and saves splitting
# LPC out into a family of its own purely to express which pair it shows.
HG_IONS_POS_NA <- list(
  PC = list(c( 86.096426, "C5H12N+, the choline cation less water"),
            c(104.106990, "choline cation, C5H14NO+"),
            c(146.981766,
              "C2H5O4PNa+, phosphocholine less trimethylamine, sodiated"))
)

# Head-group neutral losses in negative mode, which open the acyl-loss series
# the same way the positive-mode losses do.
HG_LOSS_NEG <- list(
  PS      = c(mass = 87.032028,  name = "serine residue, C3H5NO2"),
  HEXCER  = c(mass = M_HEX_RES,  name = "the hexose"),
  SHEXCER = c(mass = M_HEX_RES + 31.97207069 + 3 * E_O,   # + SO3
                     name = "the sulfated hexose")
)

# -------------------------------------------------------- candidate build ---

new_bag <- function() new.env(parent = emptyenv(), size = 256L)

add_cand <- function(bag, mz, ion, note, prio) {
  if (!is.finite(mz) || mz <= 0) return(invisible(NULL))
  k <- as.character(length(ls(bag, all.names = TRUE)) + 1L)
  assign(k, list(mz = mz, ion = ion, note = note, prio = prio), envir = bag)
  invisible(NULL)
}

# "[M+H]+" -> "[M+H", so a further loss can be appended inside the bracket.
ion_stem <- function(ion) sub("\\].$", "", ion)

bag_frame <- function(bag) {
  xs <- mget(ls(bag, all.names = TRUE), envir = bag)
  if (!length(xs)) {
    return(data.frame(mz = numeric(0), ion = character(0), note = character(0),
                      prio = integer(0), stringsAsFactors = FALSE))
  }
  d <- data.frame(
    mz   = vapply(xs, `[[`, 0,  "mz"),
    ion  = vapply(xs, `[[`, "", "ion"),
    note = vapply(xs, `[[`, "", "note"),
    prio = vapply(xs, `[[`, 0L, "prio"),
    row.names = NULL, stringsAsFactors = FALSE
  )
  # A symmetric molecule reaches the same ion by several routes -- a
  # cardiolipin with two identical chain pairs generates its phosphatidate
  # anion twice. Keep the best-ranked wording for each distinct ion.
  d <- d[order(d$prio, d$mz), , drop = FALSE]
  d[!duplicated(paste(round(d$mz, 4), d$ion)), , drop = FALSE]
}

# "16:0_18:1" for the chains left behind by a loss -- the part of a TG or PC
# annotation the user actually reads.
remainder <- function(ch, drop) {
  keep <- ch$label[-drop]
  if (!length(keep)) return("")
  paste(keep, collapse = "_")
}

#' All fragment ions predicted for one library record.
#'
#' @param meta One row of the `lipid` table (name, lipid_class, precursor_mz,
#'   adduct, ion_mode, exact_mass).
#' @return data.frame(mz, ion, note, prio); prio orders competing explanations
#'   for the same m/z, lowest first.
fragment_candidates <- function(meta) {
  bag <- new_bag()
  adduct <- as.character(meta$adduct %or% "")
  M <- neutral_mass(meta$precursor_mz, adduct,
                    if (!is.null(meta$exact_mass)) meta$exact_mass else NA_real_)
  fam <- lipid_family(meta$lipid_class)
  ch  <- parse_chains(meta$name, meta$lipid_class)

  add_cand(bag, as.numeric(meta$precursor_mz),
           if (nzchar(adduct)) paste("precursor", adduct) else "precursor ion",
           "the selected precursor, unfragmented", 0L)

  pos <- identical(as.character(meta$ion_mode), "Positive") ||
         grepl("\\+$", adduct)
  if (!is.finite(M)) return(bag_frame(bag))

  acyl  <- which(ch$kind == "acyl")
  ether <- which(ch$kind == "ether")
  lcb   <- which(ch$kind == "lcb")
  shed  <- unname(ADDUCT_SHED[adduct])

  # A "base" is an ion the acyl-loss series is enumerated from: the protonated
  # or deprotonated molecule, and whatever the head group chemistry produces
  # before any chain leaves. `pre` is the wording carried into every note built
  # on top of it, so the NH3 of an ammoniated TG is never silently dropped.
  bases <- list()
  # A base is registered for the loss series either way, but it is only worth
  # listing as its own explanation when it is not simply the precursor again:
  # for a [M-H]- record "[M-H]-, the intact ionised molecule" says nothing the
  # precursor line has not already said.
  base <- function(mz, ion, pre, prio) {
    bases[[length(bases) + 1L]] <<- list(mz = mz, ion = ion, pre = pre)
    if (nzchar(pre)) add_cand(bag, mz, ion, paste("neutral loss of", pre), prio)
  }

  if (pos) {
    # A sodiated lipid keeps its sodium all the way through: there is no
    # [M+H]+ for anything to leave from, so every loss runs off the precursor
    # itself and the ion labels have to say so.
    sodiated <- identical(adduct, "[M+Na]+")
    P    <- if (sodiated) as.numeric(meta$precursor_mz) else M + PROTON
    stem <- if (sodiated) "[M+Na" else "[M+H"

    if (sodiated) {
      base(P, "[M+Na]+", "", 2L)
    } else {
      base(P, "[M+H]+", if (!is.na(shed)) shed else "",
           if (!is.na(shed)) 1L else 2L)
    }
    pre_h <- if (!is.na(shed)) paste0(shed, " + H2O") else "H2O"
    base(P - M_H2O, paste0(stem, "-H2O]+"), pre_h, 3L)
    add_cand(bag, P - 2 * M_H2O, paste0(stem, "-2H2O]+"),
             paste("neutral loss of", sub("H2O$", "2 H2O", pre_h)), 4L)
    if (sum(ch$ox) >= 3L) {
      add_cand(bag, P - 3 * M_H2O, paste0(stem, "-3H2O]+"),
               paste("neutral loss of", sub("H2O$", "3 H2O", pre_h)), 5L)
    }

    hg <- HG_LOSS_POS[[fam]]
    if (!is.null(hg)) {
      hg_m <- as.numeric(hg[["mass"]])
      hg_n <- hg[["name"]]
      tag  <- if (fam == "HEXCER") "Hex" else "HG"
      pre_hg <- paste(c(if (!is.na(shed)) shed, hg_n), collapse = " + ")
      base(P - hg_m, sprintf("%s-%s]+", stem, tag), pre_hg, 2L)

      if (sodiated) {
        # The head group takes the sodium with it as often as it leaves as the
        # free acid: for PE that is the 163.0011 loss, which is the base peak
        # of every sodiated record in the library.
        base(P - (hg_m + M_NA - E_H), sprintf("%s-%s(Na)]+", stem, tag),
             sprintf("%s as its sodium salt, taking the sodium with it; the same ion as [M+H-%s]+",
                     hg_n, tag), 1L)
        # The same head group can keep the sodium and the charge instead.
        add_cand(bag, hg_m + M_NA - E_E, "head-group ion",
                 sprintf("sodiated %s", hg_n), 1L)
        # A sphingomyelin takes a water off the sphingoid backbone along with
        # the sodiated head group -- C5H15NO5PNa, 223.0586 -- which is the one
        # loss every sodiated SM record shows on top of the two above. The
        # product is the ceramide-like cation, minus water. Glycerophospho-
        # lipids do not do this, so it stays keyed to the SM family.
        if (fam %in% c("SM", "SM_D9")) {
          base(P - (hg_m + M_NA - E_H + M_H2O),
               sprintf("%s-%s(Na)-H2O]+", stem, tag),
               sprintf("%s as its sodium salt, taking a water off the backbone",
                       hg_n), 1L)
        }
      }
    }

    am <- HG_AMINE_LOSS[[fam]]
    if (sodiated && !is.null(am)) {
      base(P - as.numeric(am[["mass"]]),
           sprintf("%s-%s]+", stem, am[["formula"]]),
           sprintf("%s, %s", am[["formula"]], am[["name"]]), 1L)
    }

    for (h in HG_IONS_POS[[fam]] %or% list()) {
      add_cand(bag, as.numeric(h[1]), "head-group ion", h[2], 1L)
    }
    if (sodiated) for (h in HG_IONS_POS_NA[[fam]] %or% list()) {
      add_cand(bag, as.numeric(h[1]), "head-group ion", h[2], 1L)
    }
  } else {
    D <- M - PROTON
    base(D, "[M-H]-", if (!is.na(shed)) shed else "", if (!is.na(shed)) 1L else 2L)
    if (fam %in% CHOLINE_FAMS) {
      # The nine deuteriums of an SM_d9 sit on the three N-methyls, so the
      # methyl that leaves is a CD3 and the amine that follows is C3D9N.
      d9 <- identical(fam, "SM_D9")
      me <- if (d9) "CD3" else "CH3"
      pre <- paste(c(if (!is.na(shed)) shed, me), collapse = " + ")
      base(M - (if (d9) M_CD3CAT else M_CH3CAT), sprintf("[M-%s]-", me),
           paste0(pre, " (demethylation of the choline)"), 1L)
      pre <- paste(c(if (!is.na(shed)) shed, sprintf("N(%s)3", me)),
                   collapse = " + ")
      base(D - as.numeric(HG_AMINE_LOSS[[fam]][["mass"]]),
           sprintf("[M-H-N(%s)3]-", me),
           paste0(pre, " (trimethylamine from the choline)"), 2L)
    }
    hg <- HG_LOSS_NEG[[fam]]
    if (!is.null(hg)) {
      pre <- paste(c(if (!is.na(shed)) shed, hg[["name"]]), collapse = " + ")
      base(D - as.numeric(hg[["mass"]]), "[M-H-HG]-", pre, 2L)
    }
    for (h in HG_IONS_NEG[[fam]] %or% list()) {
      add_cand(bag, as.numeric(h[1]), "head-group ion", h[2], 1L)
    }
  }

  # --- chain-derived ions ---------------------------------------------------
  fa <- mass_fa(ch$c, ch$db, ch$ox)
  al <- mass_alkanol(ch$c, ch$db, ch$ox)

  for (i in acyl) {
    if (pos) {
      add_cand(bag, fa[i] - M_H2O + PROTON, sprintf("acylium %s", ch$label[i]),
               sprintf("RCO+ of %s", ch$label[i]), 2L)
      add_cand(bag, fa[i] + PROTON, sprintf("[FA %s +H]+", ch$label[i]),
               sprintf("protonated free acid of %s", ch$label[i]), 3L)
    } else {
      add_cand(bag, fa[i] - PROTON, sprintf("[FA %s -H]-", ch$label[i]),
               sprintf("carboxylate anion of %s", ch$label[i]), 1L)
      if (ch$ox[i] > 0L) {
        add_cand(bag, fa[i] - PROTON - M_H2O - (E_C + E_O),
                 sprintf("[FA %s -H-H2O-CO]-", ch$label[i]),
                 sprintf("hydroxy-acyl anion of %s, minus water and CO", ch$label[i]), 4L)
      }
    }
  }

  # Losing a chain either as the free acid or as the ketene is the workhorse
  # rule: it is what turns a TG into its DG-type ions and a PC into its LPC
  # ions. Both are enumerated from every base, so "loss of NH3 then RCOOH" and
  # "loss of the head group then RCOOH" are both covered.
  sgn <- if (pos) "+" else "-"
  for (b in bases) {
    for (i in acyl) {
      rem <- remainder(ch, i)
      tail <- if (nzchar(rem)) sprintf("; leaves %s", rem) else ""
      pre <- if (nzchar(b$pre)) paste0(b$pre, " + ") else ""
      add_cand(bag, b$mz - fa[i],
               sprintf("%s-RCOOH]%s", ion_stem(b$ion), sgn),
               sprintf("neutral loss of %s%s as the free acid (RCOOH)%s",
                       pre, ch$label[i], tail), 3L)
      add_cand(bag, b$mz - (fa[i] - M_H2O),
               sprintf("%s-ketene]%s", ion_stem(b$ion), sgn),
               sprintf("neutral loss of %s%s as the ketene (RCOOH - H2O)%s",
                       pre, ch$label[i], tail), 4L)
    }
    for (i in ether) {
      pre <- if (nzchar(b$pre)) paste0(b$pre, " + ") else ""
      add_cand(bag, b$mz - al[i],
               sprintf("%s-ROH]%s", ion_stem(b$ion), sgn),
               sprintf("neutral loss of %s%s as the alkyl alcohol",
                       pre, ch$label[i]), 4L)
    }
    # Two chains at once: how a TG reaches its monoacyl ions and a CL its
    # lyso-phosphatidate ions.
    if (length(acyl) >= 3L) {
      for (k in seq_along(acyl)) for (l in seq_along(acyl)) {
        if (l <= k) next
        i <- acyl[k]; j <- acyl[l]
        rem <- remainder(ch, c(i, j))
        tail <- if (nzchar(rem)) sprintf("; leaves %s", rem) else ""
        add_cand(bag, b$mz - fa[i] - fa[j],
                 sprintf("%s-2 RCOOH]%s", ion_stem(b$ion), sgn),
                 sprintf("neutral loss of %s and %s as free acids%s",
                         ch$label[i], ch$label[j], tail), 5L)
        for (q in c(i, j)) {
          o <- if (q == i) j else i
          add_cand(bag, b$mz - fa[q] - (fa[o] - M_H2O),
                   sprintf("%s-RCOOH-ketene]%s", ion_stem(b$ion), sgn),
                   sprintf("neutral loss of %s as the free acid and %s as the ketene%s",
                           ch$label[q], ch$label[o], tail), 5L)
        }
      }
    }
  }

  # Monoacylglycerol ions -- one chain still on the glycerol, both waters gone.
  # For a diacyl lipid this coincides with the single-chain loss, but for the
  # three-chain lipids it is a distinct ion and often the only evidence of the
  # third chain (HBMP puts three of them in every positive spectrum).
  if (pos && fam %in% c("GL", "PA", "PC", "PE", "PS", "PG", "PI", "MPE", "CL")) {
    for (i in acyl) {
      add_cand(bag, fa[i] - M_H2O + M_GLYCEROL - M_H2O + PROTON,
               sprintf("[MG %s +H-H2O]+", ch$label[i]),
               sprintf("monoacylglycerol ion of %s, minus water", ch$label[i]), 3L)
    }
  }

  # --- sphingoid base -------------------------------------------------------
  # The long-chain base carries the charge in most ceramide spectra, and its
  # water losses are the ions people actually use to call the base: 264.2686
  # is d18:1 in every positive ceramide spectrum in the library.
  d_lcb <- lcb_shift(meta$lipid_class)
  if (length(lcb) == 1L) {
    i <- lcb
    L <- mass_lcb(ch$c[i], ch$db[i], ch$ox[i]) + d_lcb
    lbl <- ch$label[i]
    if (pos) {
      add_cand(bag, L + PROTON, sprintf("[LCB %s +H]+", lbl),
               sprintf("protonated sphingoid base %s", lbl), 2L)
      add_cand(bag, L + PROTON - M_H2O, sprintf("[LCB %s +H-H2O]+", lbl),
               sprintf("sphingoid base %s, minus one water", lbl), 1L)
      add_cand(bag, L + PROTON - 2 * M_H2O, sprintf("[LCB %s +H-2H2O]+", lbl),
               sprintf("sphingoid base %s, minus two waters", lbl), 1L)
      add_cand(bag, L + PROTON - M_H2O - M_CH2O,
               sprintf("[LCB %s +H-H2O-CH2O]+", lbl),
               sprintf("sphingoid base %s, minus water and formaldehyde", lbl), 3L)
      add_cand(bag, L + PROTON - M_H2O - M_C2H5N,
               sprintf("[LCB %s +H-H2O-C2H5N]+", lbl),
               sprintf("sphingoid base %s, minus water and C2H5N", lbl), 3L)
      if (ch$ox[lcb] >= 3L) {
        add_cand(bag, L + PROTON - 3 * M_H2O, sprintf("[LCB %s +H-3H2O]+", lbl),
                 sprintf("phytosphingoid base %s, minus three waters", lbl), 1L)
      }
      add_cand(bag, L + PROTON - 2 * M_H2O - M_CH2O,
               sprintf("[LCB %s +H-2H2O-CH2O]+", lbl),
               sprintf("sphingoid base %s, minus two waters and formaldehyde", lbl), 4L)
      for (j in acyl) {
        amide <- mass_fa(ch$c[j], ch$db[j], ch$ox[j]) + E_N + E_H - E_O
        add_cand(bag, amide + PROTON, sprintf("[N-acyl %s amide +H]+", ch$label[j]),
                 sprintf("N-acyl chain %s as its protonated amide", ch$label[j]), 2L)
      }
    } else {
      add_cand(bag, L - PROTON, sprintf("[LCB %s -H]-", lbl),
               sprintf("sphingoid base %s", lbl), 2L)
      if (cer_emits("lcb_h2o", meta$lipid_class)) {
        add_cand(bag, L - PROTON - M_H2O, sprintf("[LCB %s -H-H2O]-", lbl),
                 sprintf("sphingoid base %s, minus one water", lbl), 1L)
      }
      add_cand(bag, L - PROTON - M_H2O - M_NH3,
               sprintf("[LCB %s -H-H2O-NH3]-", lbl),
               sprintf("sphingoid base %s, minus water and ammonia", lbl), 3L)
      add_cand(bag, L - PROTON - M_H2O - M_C2H5N,
               sprintf("[LCB %s -H-H2O-C2H5N]-", lbl),
               sprintf("sphingoid base %s, minus water and C2H5N", lbl), 4L)
      if (cer_emits("lcb_ch3oh_nh3", meta$lipid_class)) {
        add_cand(bag, L - PROTON - M_CH4O - M_NH3,
                 sprintf("[LCB %s -H-CH3OH-NH3]-", lbl),
                 sprintf("phytosphingoid base %s, minus methanol from C1 and the amine as ammonia",
                         lbl), 3L)
      }
      if (ch$ox[lcb] >= 3L) {
        add_cand(bag, L - PROTON - M_H2O - M_C2H5N - M_CH2O,
                 sprintf("[LCB %s -H-H2O-C2H5N-CH2O]-", lbl),
                 sprintf("phytosphingoid base %s, minus water, C2H5N and formaldehyde",
                         lbl), 3L)
      }
    }
    # The N-acyl chain leaves as an amide rather than a fatty acid, so it needs
    # its own ion: for Cer 18:1;O2/16:0 that is m/z 254.2489, not 255.2330.
    # It also shows up still carrying C1-C2 of the sphingoid base, cleaved at
    # the C2-C3 bond -- which is what m/z 280.2646 and 296.2595 really are in
    # an N-acyl ceramide, not the base fragments they coincide with there.
    if (!pos) {
      cls <- meta$lipid_class
      for (j in acyl) {
        amide <- fa[j] + E_N + E_H - E_O
        lb <- ch$label[j]
        add_cand(bag, amide - PROTON, sprintf("[N-acyl %s amide -H]-", lb),
                 sprintf("N-acyl chain %s as its amide anion", lb), 2L)
        if (cer_emits("amide_nh3", cls)) {
          add_cand(bag, amide - PROTON - M_NH3,
                   sprintf("[N-acyl %s amide -H-NH3]-", lb),
                   sprintf("N-acyl chain %s, the amide anion less ammonia", lb), 3L)
        }
        if (cer_emits("amide_c2h2", cls)) {
          add_cand(bag, amide - PROTON + M_C2H2,
                   sprintf("[N-acyl %s + C1-C2]-", lb),
                   sprintf("N-acyl chain %s still carrying C1-C2 of the sphingoid base",
                           lb), 2L)
        }
        if (cer_emits("amide_c3h4o", cls)) {
          add_cand(bag, amide - PROTON + M_C3H4O,
                   sprintf("[N-acyl %s + C1-C3]-", lb),
                   sprintf("N-acyl chain %s carrying C1-C3 of the phytosphingoid base and a hydroxyl",
                           lb), 2L)
        }
        if (cer_emits("amide_c2h2o", cls)) {
          add_cand(bag, amide - PROTON + M_C2H2 + E_O,
                   sprintf("[N-acyl %s + C1-C2+OH]-", lb),
                   sprintf("N-acyl chain %s with C1-C2 of the sphingoid base and the C1 hydroxyl",
                           lb), 2L)
        }
      }
    }
  }
  # --- ceramide core --------------------------------------------------------
  # Every sphingolipid carrying a head group -- a hexose, a sulfated hexose, an
  # acylated hexose, a phosphocholine -- sheds it to leave the bare ceramide,
  # and that ion with its water losses is how the backbone gets read off the
  # spectrum. Building the core from the base and the N-acyl chains, instead of
  # subtracting a per-class head-group mass, covers AHexCer, ASHexCer and the
  # Cer_E** families without having to name any of them.
  core_acyl <- which(ch$kind == "acyl" & !ch$on_hg)
  if (length(lcb) == 1L && length(core_acyl)) {
    L <- mass_lcb(ch$c[lcb], ch$db[lcb], ch$ox[lcb]) + d_lcb
    core <- L + sum(fa[core_acyl] - M_H2O)
    core_lbl <- paste(c(ch$label[lcb], ch$label[core_acyl]), collapse = "/")
    if (abs(core - M) > 0.5) {          # otherwise it is only the precursor
      if (pos) {
        add_cand(bag, core + PROTON, "[Cer+H]+",
                 sprintf("ceramide core %s, head group lost", core_lbl), 3L)
        add_cand(bag, core + PROTON - M_H2O, "[Cer+H-H2O]+",
                 sprintf("ceramide core %s, head group and water lost", core_lbl), 3L)
        add_cand(bag, core + PROTON - 2 * M_H2O, "[Cer+H-2H2O]+",
                 sprintf("ceramide core %s, head group and two waters lost", core_lbl), 4L)
      } else {
        add_cand(bag, core - PROTON, "[Cer-H]-",
                 sprintf("ceramide core %s, head group lost", core_lbl), 3L)
        add_cand(bag, core - PROTON - M_H2O, "[Cer-H-H2O]-",
                 sprintf("ceramide core %s, head group and water lost", core_lbl), 3L)
      }
    }
  }

  # The fatty acid esterified to the hexose of an AHexCer departs as an
  # oxocarbenium that still carries the sugar.
  if (pos) for (i in which(ch$on_hg)) {
    add_cand(bag, fa[i] - M_H2O + M_HEX_RES + PROTON,
             sprintf("[Hex+%s]+", ch$label[i]),
             sprintf("hexose still carrying its esterified %s", ch$label[i]), 2L)
  }

  # Ceramide backbone losses, which are what distinguishes a ceramide
  # [M-H]- spectrum from a plain fatty-acyl one.
  if (!pos && fam %in% c("CER", "HEXCER", "SHEXCER") && length(lcb) == 1L) {
    for (b in bases) {
      pre <- if (nzchar(b$pre)) paste0(b$pre, " + ") else ""
      stem <- ion_stem(b$ion)
      add_cand(bag, b$mz - M_CH2O, paste0(stem, "-CH2O]-"),
               paste0("loss of ", pre, "formaldehyde from C1"), 3L)
      add_cand(bag, b$mz - M_CH4O, paste0(stem, "-CH3OH]-"),
               paste0("loss of ", pre, "methanol from C1"), 4L)
      add_cand(bag, b$mz - M_H2O - M_CH2O, paste0(stem, "-H2O-CH2O]-"),
               paste0("loss of ", pre, "water and formaldehyde"), 4L)
      if (cer_emits("m_h2o", meta$lipid_class)) {
        add_cand(bag, b$mz - M_H2O, paste0(stem, "-H2O]-"),
                 paste0("loss of ", pre, "water"), 3L)
      }
      if (cer_emits("m_2h2o", meta$lipid_class)) {
        add_cand(bag, b$mz - 2 * M_H2O, paste0(stem, "-2H2O]-"),
                 paste0("loss of ", pre, "two waters"), 4L)
      }
    }
  }

  # The counterpart of the [LCB-H-H2O-C2H5N]- ion: the same bond breaks, but the
  # charge stays with the head group and the N-acyl chain. Both halves show up
  # in the hexosyl and sulfatide records, and seeing them as a pair is what
  # confirms the base assignment.
  if (!pos && length(lcb) == 1L) {
    L <- mass_lcb(ch$c[lcb], ch$db[lcb], ch$ox[lcb]) + d_lcb
    lost <- L - (2 * E_C + 7 * E_H + E_N + E_O)      # LCB less C2H7NO
    D <- M - PROTON
    add_cand(bag, D - lost, "[M-H-LCB residue]-",
             sprintf("loss of the sphingoid base %s less C2H7NO; head group and N-acyl retained",
                     ch$label[lcb]), 3L)
    add_cand(bag, D - lost - M_H2O, "[M-H-LCB residue-H2O]-",
             sprintf("loss of the sphingoid base %s less C2H7NO, and water",
                     ch$label[lcb]), 4L)
  }

  # --- cardiolipin phosphatidate ions ---------------------------------------
  # A cardiolipin fragments into its two phosphatidate halves; those anions,
  # and their own acyl losses, are the whole basis of assigning CL chains.
  if (!pos && fam == "CL" && length(acyl) >= 2L) {
    res <- fa - M_H2O
    for (k in seq_along(acyl)) {
      i <- acyl[k]
      add_cand(bag, M_GPA + res[i] - PROTON, sprintf("[LPA %s -H]-", ch$label[i]),
               sprintf("lyso-phosphatidate carrying %s", ch$label[i]), 4L)
      for (l in seq_along(acyl)) {
        if (l <= k) next
        j <- acyl[l]
        pa <- M_GPA + res[i] + res[j] - PROTON
        # Adjacent chains in the name are the two on one glycerol; any other
        # pair would need a rearrangement, so rank it lower. The pair is named
        # in a fixed order so the two routes to a symmetric cardiolipin's
        # phosphatidate collapse to one entry.
        canonical <- (l == k + 1L) && (k %% 2L == 1L)
        pair <- paste(sort(c(ch$label[i], ch$label[j])), collapse = "_")
        add_cand(bag, pa, sprintf("[PA %s -H]-", pair),
                 sprintf("phosphatidate half carrying %s and %s",
                         ch$label[i], ch$label[j]), if (canonical) 2L else 5L)
        for (q in c(i, j)) {
          other <- if (q == i) j else i
          add_cand(bag, pa - fa[q],
                   sprintf("[PA %s -H-RCOOH]-", pair),
                   sprintf("phosphatidate half minus %s; leaves %s",
                           ch$label[q], ch$label[other]),
                   if (canonical) 3L else 6L)
        }
      }
    }
  }

  bag_frame(bag)
}

# ------------------------------------------------------------- matching -----

#' Explain each peak of a spectrum, as HTML ready for a plotly tooltip.
#'
#' @param tol Half-window in Da. Deliberately tighter than the app's spectral
#'   match tolerance: a candidate list has many near-coincident entries, and a
#'   loose window would hand the user a confident-looking wrong assignment.
#' @param max_n How many competing explanations to show for one peak.
#' @param width Wrap column. A plotly tooltip does not wrap for itself, so an
#'   unwrapped "neutral loss of HCOOH + CH3 ..." note would stretch the hover
#'   box across the whole card.
#' @return character vector, one entry per row of `peaks`; "" where no
#'   candidate matched.
annotate_peaks <- function(meta, peaks, tol = 0.005, max_n = 2L, width = 46L) {
  out <- character(nrow(peaks))
  if (nrow(peaks) == 0L) return(out)
  cand <- tryCatch(fragment_candidates(meta), error = function(e) NULL)
  if (is.null(cand) || nrow(cand) == 0L) return(out)

  cand <- cand[order(cand$mz), , drop = FALSE]
  lo <- findInterval(peaks$mz - tol, cand$mz)          # last row strictly below
  hi <- findInterval(peaks$mz + tol, cand$mz)          # last row within
  for (i in seq_len(nrow(peaks))) {
    if (hi[i] <= lo[i]) next
    j <- (lo[i] + 1L):hi[i]
    j <- j[order(cand$prio[j], abs(cand$mz[j] - peaks$mz[i]))]
    j <- j[!duplicated(cand$ion[j])]
    # A second explanation earns its place only if it is as well ranked as the
    # first. Otherwise every acyl loss would drag along the water-plus-ketene
    # route to the same ion, which is the same chemistry stated worse.
    j <- utils::head(j[cand$prio[j] == cand$prio[j[1]]], max_n)
    out[i] <- paste(sprintf("<b>%s</b><br>%s", cand$ion[j],
                            wrap_html(cand$note[j], width)),
                    collapse = "<br><br>")
  }
  out
}

wrap_html <- function(x, width) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "<br>"),
         "", USE.NAMES = FALSE)
}
