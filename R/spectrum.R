# Rotated peak-label geometry. Plotly's yshift is in pixels while the axis
# range is in percent, so converting needs an assumed plot height (~230px for
# the default 0-145% range): 1% is about 1.59px.
LANE_PX   <- 46   # vertical pitch between stacked label lanes, px
LANE_PCT  <- 29   # the same pitch expressed in % of base peak
LABEL_PCT <- 25   # length of a rotated "806.592" label, in % of base peak

# plotly event source ids. The two tabs each own one so a zoom on the library
# browser's spectrum does not re-lay the search tab's mirror plot, and vice
# versa -- event_data() is keyed on this and nothing else.
SPEC_SOURCE  <- "spectrum"
QSPEC_SOURCE <- "query-spectrum"

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

#' m/z window that comfortably contains one spectrum, or a mirrored pair.
#'
#' The window must cover both spectra so the mirror halves stay aligned. A
#' pasted query spectrum has no precursor, so NA precursors drop out rather
#' than poisoning the range.
spectrum_xlim <- function(top, ref = NULL) {
  mzs <- c(top$peaks$mz, top$meta$precursor_mz)
  if (!is.null(ref)) mzs <- c(mzs, ref$peaks$mz, ref$meta$precursor_mz)
  mzs <- mzs[is.finite(mzs)]
  if (!length(mzs)) return(c(0, 1000))
  rng <- range(mzs)
  pad <- max(diff(rng) * 0.06, 5)
  c(rng[1] - pad, rng[2] + pad)
}

#' A plotly placeholder for "nothing selected yet".
spectrum_placeholder <- function(text) {
  plot_ly(type = "scatter", mode = "markers") |>
    layout(
      xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
      annotations = list(
        text = text, showarrow = FALSE, xref = "paper", yref = "paper",
        x = 0.5, y = 0.5, font = list(size = 14, color = "#7b8a8b")
      )
    ) |> config(displayModeBar = FALSE)
}

#' Draw one spectrum, or two mirrored against each other.
#'
#' Shared by the library browser and the spectrum search so the two views
#' cannot drift apart; `source` keeps their zoom events separate.
#'
#' @param top,ref   Spectra in the shape get_spectrum() returns. `ref` NULL
#'   draws a single, upward spectrum.
#' @param match_a,match_b Indices of the matched peaks in `top` / `ref`.
#'   NULL means "treat every peak of `top` as matched", which is what a lone
#'   spectrum wants.
#' @param ann_top,ann_ref Per-peak explanation strings, or NULL for none.
#' @return A plotly object.
spectrum_plot <- function(top, ref = NULL, match_a = NULL, match_b = integer(0),
                          ann_top = NULL, ann_ref = NULL,
                          label_peaks = TRUE, label_n = 8,
                          show_precursor = TRUE, source = SPEC_SOURCE,
                          xlim = NULL) {
  if (is.null(xlim)) xlim <- spectrum_xlim(top, ref)
  if (is.null(match_a)) match_a <- seq_len(nrow(top$peaks))
  if (is.null(ann_top)) ann_top <- rep("", nrow(top$peaks))
  if (!is.null(ref) && is.null(ann_ref)) ann_ref <- rep("", nrow(ref$peaks))

  # Labels are laid out before the traces because the axis range, and so the
  # precursor line's height, depends on how many lanes they need.
  lay <- spectrum_layout(top, ref, xlim, label_n, enabled = isTRUE(label_peaks))
  y_top <- lay$y_top
  y_bot <- lay$y_bot

  # One spectrum, drawn upward (dir = 1) or mirrored downward (dir = -1).
  # Matched peaks keep the spectrum's colour, unmatched fall back to grey.
  add_spectrum <- function(p, sp, dir, color, matched_idx, ann) {
    pk <- sp$peaks
    is_match <- seq_len(nrow(pk)) %in% matched_idx
    for (grp in list(list(i = !is_match, col = UNMATCHED),
                     list(i = is_match,  col = color))) {
      if (!any(grp$i)) next
      q <- pk[grp$i, , drop = FALSE]
      qa <- ann[grp$i]
      p <- p |>
        add_segments(
          x = q$mz, xend = q$mz, y = 0, yend = dir * q$rel,
          line = list(color = grp$col, width = 1.5),
          hoverinfo = "none", showlegend = FALSE
        ) |>
        add_markers(
          x = q$mz, y = dir * q$rel,
          marker = list(color = grp$col, size = 5, opacity = 0.01),
          hoverinfo = "text", showlegend = FALSE,
          text = paste0(
            sprintf("%s<br>m/z %.4f<br>%.1f%% base peak<br>abs %.0f",
                    sp$meta$name, q$mz, q$rel, q$intensity),
            ifelse(nzchar(qa), paste0("<br><br>", qa), ""))
        )
    }
    # A pasted query spectrum need not have a precursor at all.
    if (isTRUE(show_precursor) && is.finite(sp$meta$precursor_mz)) {
      p <- p |> add_segments(
        x = sp$meta$precursor_mz, xend = sp$meta$precursor_mz, y = 0,
        yend = dir * (if (dir > 0) y_top else y_bot) * 0.95,
        line = list(color = "#e74c3c", width = 1, dash = "dot"),
        hoverinfo = "text", showlegend = FALSE,
        text = sprintf("precursor m/z %.4f", sp$meta$precursor_mz)
      )
    }
    p
  }

  p <- plot_ly(source = source) |>
    add_spectrum(top, 1, TOP_COLOR, match_a, ann_top)
  if (!is.null(ref)) {
    p <- p |> add_spectrum(ref, -1, REF_COLOR, match_b, ann_ref) |>
      add_segments(x = xlim[1], xend = xlim[2], y = 0, yend = 0,
                   line = list(color = "#95a5a6", width = 1),
                   hoverinfo = "none", showlegend = FALSE)
  }

  # In mirror mode the axis runs negative, but intensity is a magnitude on
  # both halves, so the tick labels stay positive.
  ticks <- if (is.null(ref)) c(0, 25, 50, 75, 100) else c(-100, -50, 0, 50, 100)

  p |> layout(
    xaxis = list(title = "m/z", range = xlim, zeroline = FALSE),
    # Headroom for the rotated peak labels above a 100% base peak.
    # ticktext overrides ticksuffix, so the "%" has to be baked in here.
    yaxis = list(title = "relative intensity (%)",
                 range = c(-y_bot, y_top), zeroline = FALSE,
                 tickvals = ticks, ticktext = paste0(abs(ticks), "%")),
    annotations = lay$annotations, hovermode = "closest",
    # Fragment explanations run to several lines; centred text makes them
    # much harder to read than the two-line default tooltip was.
    hoverlabel = list(align = "left"),
    margin = list(t = 20, r = 10)
  ) |> config(displaylogo = FALSE,
              modeBarButtonsToRemove = list("select2d", "lasso2d"))
}

#' Re-lay the peak labels for the window a zoom/pan event left behind.
#'
#' Labels are laid out server-side at render, so a zoom would otherwise keep
#' the spacing computed for the whole spectrum and never resolve a crowded
#' region.
#'
#' @return A list for plotly's relayout, or NULL when the event is not a range
#'   change -- our own annotation push and a resize both arrive here too, and
#'   acting on them would loop.
spectrum_relayout <- function(ed, top, ref = NULL, label_n = 8,
                              label_peaks = TRUE) {
  if (is.null(ed)) return(NULL)
  full <- spectrum_xlim(top, ref)

  # A drag zoom or pan sends the two bounds as separate keys; a programmatic
  # relayout sends them as one array; double-click to reset sends autorange.
  xlim <- if (!is.null(ed[["xaxis.range[0]"]]) && !is.null(ed[["xaxis.range[1]"]])) {
    sort(c(as.numeric(ed[["xaxis.range[0]"]]), as.numeric(ed[["xaxis.range[1]"]])))
  } else if (length(ed[["xaxis.range"]]) == 2L) {
    sort(as.numeric(ed[["xaxis.range"]]))
  } else if (isTRUE(ed[["xaxis.autorange"]])) {
    full
  } else {
    return(NULL)
  }

  lay <- spectrum_layout(top, ref, xlim, label_n, enabled = isTRUE(label_peaks))
  upd <- list(annotations = lay$annotations)
  # Only reclaim vertical headroom if the user has not set their own y range
  # (box zoom changes both axes); overriding it would fight them.
  if (!any(grepl("^yaxis\\.range", names(ed)))) {
    upd[["yaxis.range"]] <- c(-lay$y_bot, lay$y_top)
  }
  upd
}
