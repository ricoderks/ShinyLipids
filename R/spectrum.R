# Rotated peak-label geometry. Plotly's yshift is in pixels while the axis
# range is in percent, so converting needs an assumed plot height (~230px for
# the default 0-145% range): 1% is about 1.59px.
LANE_PX   <- 46   # vertical pitch between stacked label lanes, px
LANE_PCT  <- 29   # the same pitch expressed in % of base peak
LABEL_PCT <- 25   # length of a rotated "806.592" label, in % of base peak

SPEC_SOURCE <- "spectrum"   # plotly event source id for the spectrum plot

# Mirror-plot palette. Matched fragments keep their spectrum's colour so the
# top/bottom identity survives; unmatched ones drop back to grey, which makes
# the shared fragments the thing the eye lands on.
TOP_COLOR   <- "#2c3e50"
REF_COLOR   <- "#8e44ad"
UNMATCHED   <- "#b3bcc4"

#' Lay out rotated m/z labels for the peaks visible in an m/z window.
#'
#' Labels are rotated vertically: horizontal ones need ~60 Da of clear space at
#' a typical plot width, which would suppress exactly the peaks that matter
#' most -- the two fatty acyl anions of a diacyl lipid sit ~20 Da apart
#' (283.264 / 303.233 for PE 18:0_20:4).
#'
#' "Max labels" is a promise: take the `label_n` most intense peaks *in the
#' window* and label all of them. Labels that would still collide are stacked
#' into a higher lane rather than dropped -- dropping made the control look
#' broken, since near-coincident fragments are common (PE 18:0_20:4 has peaks
#' 1.96 Da apart, roughly one pixel; PC 16:0_18:0 [M+HCOO]- two 10.02 Da apart).
#'
#' Spacing is a fraction of the *visible* span, not the full spectrum, so
#' zooming in genuinely resolves crowded regions instead of keeping the layout
#' that was computed for the whole spectrum.
#'
#' @param dir  1 draws labels upward from each peak, -1 downward for the
#'   mirrored reference spectrum.
#' @param color Label colour.
#' @return list(annotations, y_top) -- y_top is a magnitude in both directions.
spectrum_labels <- function(pk, xlim, label_n, enabled = TRUE, dir = 1,
                            color = TOP_COLOR) {
  none <- list(annotations = list(), y_top = 145)
  if (!enabled || nrow(pk) == 0) return(none)

  vis <- which(pk$mz >= xlim[1] & pk$mz <= xlim[2])
  if (!length(vis)) return(none)

  sel <- vis[utils::head(order(-pk$rel[vis]), label_n)]
  sel <- sel[order(pk$mz[sel])]                        # place left to right
  min_gap <- max(xlim[2] - xlim[1], 1e-9) * 0.02       # ~1 rotated label wide

  lane_last <- numeric(0)   # rightmost m/z already labelled, per lane
  lab_rel <- numeric(0); lab_lane <- integer(0); ann <- list()
  for (i in sel) {
    lane <- which(pk$mz[i] - lane_last > min_gap)[1]
    if (is.na(lane)) {
      lane_last <- c(lane_last, -Inf)
      lane <- length(lane_last)
    }
    lane_last[lane] <- pk$mz[i]
    lab_rel  <- c(lab_rel, pk$rel[i])
    lab_lane <- c(lab_lane, lane)
    ann[[length(ann) + 1L]] <- list(
      x = pk$mz[i], y = dir * pk$rel[i], text = sprintf("%.3f", pk$mz[i]),
      showarrow = FALSE, textangle = -90,
      yanchor = if (dir > 0) "bottom" else "top",
      yshift = dir * (5 + (lane - 1) * LANE_PX),
      font = list(size = 10, color = color)
    )
  }

  # Headroom is driven by the tallest peak that actually carries a label in
  # each lane, not by the lane count: a second lane sitting over two 5% peaks
  # needs no extra room. Reserving per lane regardless squashed the spectrum
  # into the bottom of the card for the ~7% of spectra needing three lanes.
  list(annotations = ann,
       y_top = max(145, max(lab_rel + (lab_lane - 1) * LANE_PCT + LABEL_PCT) + 5))
}

#' Pair up peaks between two spectra that agree in m/z within `tol` Da.
#'
#' Greedy nearest-neighbour, walking the more intense peaks first so a strong
#' fragment claims its partner before a weak neighbour can take it. Each peak
#' is used at most once, so two peaks 1.96 Da apart cannot both claim the same
#' partner at a 0.01 Da tolerance.
#'
#' @return list(a, b) of matched index vectors, aligned pairwise.
match_peaks <- function(pk_a, pk_b, tol = 0.01) {
  if (nrow(pk_a) == 0 || nrow(pk_b) == 0) {
    return(list(a = integer(0), b = integer(0)))
  }
  taken <- rep(FALSE, nrow(pk_b))
  ia <- integer(0); ib <- integer(0)
  for (i in order(-pk_a$rel)) {
    d <- abs(pk_b$mz - pk_a$mz[i])
    d[taken] <- Inf
    j <- which.min(d)
    if (length(j) == 1L && is.finite(d[j]) && d[j] <= tol) {
      taken[j] <- TRUE
      ia <- c(ia, i); ib <- c(ib, j)
    }
  }
  list(a = ia, b = ib)
}

#' Cosine similarity between two spectra over their matched peaks.
#'
#' The usual normalised dot product on relative intensities: unmatched peaks
#' contribute nothing to the numerator but still count in each spectrum's norm,
#' so a spectrum with many unshared fragments is penalised. Returns NA when
#' either spectrum is empty.
spectral_similarity <- function(pk_a, pk_b, m) {
  if (nrow(pk_a) == 0 || nrow(pk_b) == 0) return(NA_real_)
  denom <- sqrt(sum(pk_a$rel^2)) * sqrt(sum(pk_b$rel^2))
  if (denom == 0) return(NA_real_)
  sum(pk_a$rel[m$a] * pk_b$rel[m$b]) / denom
}

#' Full annotation set and y range for one spectrum or a mirrored pair.
#'
#' Kept in one place because the initial render and the zoom handler must agree
#' exactly -- plotly's relayout replaces *all* annotations, so the spectrum name
#' captions have to be rebuilt alongside the peak labels on every zoom.
#'
#' @return list(annotations, y_top, y_bot)
spectrum_layout <- function(top, ref, xlim, label_n, enabled = TRUE) {
  lt <- spectrum_labels(top$peaks, xlim, label_n, enabled, dir = 1,
                        color = if (is.null(ref)) TOP_COLOR else TOP_COLOR)
  ann <- lt$annotations
  y_top <- lt$y_top
  y_bot <- 0

  if (!is.null(ref)) {
    lr <- spectrum_labels(ref$peaks, xlim, label_n, enabled, dir = -1,
                          color = REF_COLOR)
    ann <- c(ann, lr$annotations)
    y_bot <- lr$y_top

    caption <- function(text, color, y, anchor) list(
      x = 0.01, y = y, xref = "paper", yref = "paper", text = text,
      showarrow = FALSE, xanchor = "left", yanchor = anchor,
      font = list(size = 11, color = color)
    )
    ann <- c(ann, list(
      caption(top$meta$name, TOP_COLOR, 0.99, "top"),
      caption(ref$meta$name, REF_COLOR, 0.01, "bottom")
    ))
  }
  list(annotations = ann, y_top = y_top, y_bot = y_bot)
}

#' Render a spectrum's peaks as tab-separated text for the clipboard.
#'
#' Tab-separated with a header row so it pastes straight into Excel as three
#' columns. Peaks are in m/z order, matching the plot. m/z carries 4 decimals:
#' the source values are float32, so more digits would be invented precision.
peak_list_tsv <- function(s) {
  if (is.null(s) || nrow(s$peaks) == 0) return("")
  p <- s$peaks[order(s$peaks$mz), , drop = FALSE]
  paste0(
    "m/z\tintensity\trel_intensity_pct\n",
    paste(sprintf("%.4f\t%.0f\t%.2f", p$mz, p$intensity, p$rel), collapse = "\n"),
    "\n"
  )
}
