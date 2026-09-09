## ShinyLipids -- browse and visualise MS/MS spectra from an MS-DIAL lipid library.
##
## The library lives in DB/lipids.sqlite, built from the .lbm2 file by
## scripts/lbm2_to_sqlite.py. See R/db.R for the query layer and
## R/spectrum.R for the peak-label layout.

library(shiny)
library(bslib)
library(plotly)
library(DT)

source("R/db.R")
source("R/spectrum.R")
source("R/fragments.R")

con <- lipid_db()
onStop(function() DBI::dbDisconnect(con))

CLASSES   <- lipid_classes(con)
ADDUCTS   <- lipid_adducts(con)
N_RECORDS <- lipid_count(con)

MAX_ROWS <- 500

APP_TITLE <- "ShinyLipids — MS/MS spectrum browser - v0.3.2"
LIBRARY_FILE <- "Msp20251120132005_NCDK_conventional_converted_dev.lbm2"


# ---------------------------------------------------------------- UI --------

ui <- page_sidebar(
  # A tag (rather than a bare string) as the title means bslib drops it into the
  # navbar's container-fluid untouched; that container is a flex row with
  # space-between, so the About link lands on the right-hand side by itself.
  title = tagList(
    h1(APP_TITLE, class = "bslib-page-title navbar-brand"),
    tags$ul(
      class = "navbar-nav ms-auto",
      tags$li(
        class = "nav-item",
        actionLink("about", "About", class = "nav-link",
                   icon = bsicons::bs_icon("info-circle"))
      )
    )
  ),
  window_title = APP_TITLE,
  theme = bs_theme(version = 5, preset = "flatly"),
  class = "bslib-page-dashboard",
  fillable = TRUE,

  # Keep the results table inside its card: fixed layout plus ellipsis on the
  # name column, rather than letting DT overflow into a horizontal scrollbar.
  tags$head(tags$style(HTML("
    #hits table { table-layout: fixed; width: 100% !important; }
    #hits td, #hits th {
      overflow: hidden; text-overflow: ellipsis; white-space: nowrap;
    }
    #hits tbody tr { cursor: pointer; }
    .dataTables_scrollBody { border-bottom: none !important; }
    #copy_peaks:disabled { opacity: .35; }
  "))),

  # The peak list is pushed to the client whenever the selection changes and
  # copied synchronously on click. Doing the round trip on click instead would
  # risk the browser rejecting the write: navigator.clipboard needs transient
  # user activation, which a server round trip can outlive.
  tags$head(tags$script(HTML("
    document.addEventListener('DOMContentLoaded', function () {
      var peakText = '';

      Shiny.addCustomMessageHandler('spectrumPeaks', function (msg) {
        peakText = msg.text || '';
        var b = document.getElementById('copy_peaks');
        if (b) b.disabled = !peakText;
      });

      // navigator.clipboard only exists in a secure context, so an app served
      // over plain http from an internal host needs the textarea fallback.
      function legacyCopy(text) {
        var ta = document.createElement('textarea');
        ta.value = text;
        ta.setAttribute('readonly', '');
        ta.style.position = 'fixed';
        ta.style.opacity = '0';
        document.body.appendChild(ta);
        ta.select();
        var ok = false;
        try { ok = document.execCommand('copy'); } catch (e) { ok = false; }
        document.body.removeChild(ta);
        return ok;
      }

      function flash(btn, msg) {
        if (btn.dataset.busy) return;
        btn.dataset.busy = '1';
        var original = btn.innerHTML;
        btn.innerHTML = '<span class=\"small\">' + msg + '</span>';
        setTimeout(function () {
          btn.innerHTML = original;
          delete btn.dataset.busy;
        }, 1400);
      }

      document.addEventListener('click', function (e) {
        var btn = e.target.closest ? e.target.closest('#copy_peaks') : null;
        if (!btn || !peakText) return;
        if (navigator.clipboard && window.isSecureContext) {
          navigator.clipboard.writeText(peakText).then(
            function () { flash(btn, 'Copied'); },
            function () { flash(btn, legacyCopy(peakText) ? 'Copied' : 'Copy failed'); }
          );
        } else {
          flash(btn, legacyCopy(peakText) ? 'Copied' : 'Copy failed');
        }
      });
    });
  "))),

  sidebar = sidebar(
    width = 320,
    title = "Search",

    textInput("name", "Lipid name contains", placeholder = "e.g. PC 34:1"),

    layout_columns(
      col_widths = c(7, 5),
      numericInput("mz", "Precursor m/z", value = NA, min = 0, step = 0.001),
      numericInput("tol", "Tolerance", value = 0.01, min = 0, step = 0.001)
    ),
    radioButtons("tol_unit", NULL, choices = c("Da" = "da", "ppm" = "ppm"),
                 selected = "da", inline = TRUE),

    selectizeInput(
      "class", "Lipid class", choices = NULL, multiple = TRUE,
      options = list(placeholder = "All classes", plugins = list("remove_button"))
    ),

    radioButtons("mode", "Ion mode",
                 choices = c("Any" = "", "Positive" = "Positive",
                             "Negative" = "Negative"),
                 selected = "", inline = TRUE),

    selectizeInput(
      "adduct", "Adduct", choices = NULL, multiple = TRUE,
      options = list(placeholder = "All adducts", plugins = list("remove_button"))
    ),

    actionButton("reset", "Reset filters", class = "btn-outline-secondary btn-sm",
                 icon = icon("rotate-left"))
  ),

  layout_column_wrap(
    width = 1/3, fill = FALSE, heights_equal = "row",
    value_box(
      title = "Library", value = format(N_RECORDS, big.mark = ","),
      showcase = bsicons::bs_icon("database"), theme = "primary",
      paste(nrow(CLASSES), "lipid classes")
    ),
    value_box(
      title = "Matches", value = textOutput("n_matches", inline = TRUE),
      showcase = bsicons::bs_icon("funnel"), theme = "secondary",
      textOutput("match_note", inline = TRUE)
    ),
    value_box(
      title = "Selected",
      # Lipid names run long (e.g. "CL 12:0_18:0_28:0_22:5"); shrink the value
      # so a long name does not stretch the whole KPI row.
      value = tags$span(style = "font-size: 1.2rem; line-height: 1.3;",
                        textOutput("sel_name", inline = TRUE)),
      showcase = bsicons::bs_icon("graph-up"), theme = "success",
      textOutput("sel_sub", inline = TRUE)
    )
  ),

  layout_columns(
    col_widths = breakpoints(sm = 12, lg = c(6, 6)),

    card(
      full_screen = TRUE, min_height = 320,
      card_header(
        "Results",
        tooltip(
          bsicons::bs_icon("info-circle", title = "About the result list"),
          paste("Click a row to plot its MS/MS spectrum. At most",
                MAX_ROWS, "rows are returned; narrow the filters to see more.")
        ),
        class = "d-flex justify-content-between align-items-center"
      ),
      DTOutput("hits")
    ),

    card(
      full_screen = TRUE, min_height = 320,
      card_header(
        textOutput("spec_title", inline = TRUE),
        tags$div(
          class = "d-flex align-items-center gap-3",
          uiOutput("compare_btn", inline = TRUE),
          tooltip(
            tags$button(
              id = "copy_peaks", type = "button", disabled = NA,
              class = "btn btn-sm btn-link p-0 border-0 text-decoration-none",
              `aria-label` = "Copy peak list",
              bsicons::bs_icon("clipboard")
            ),
            "Copy the peak list as tab-separated text"
          ),
          popover(
            bsicons::bs_icon("gear", title = "Plot options"),
            title = "Plot options",
            checkboxInput("label_peaks", "Label peaks", TRUE),
            sliderInput("label_n", "Max labels", min = 1, max = 20, value = 8),
            checkboxInput("show_precursor", "Mark precursor m/z", TRUE),
            checkboxInput("explain_peaks", "Explain fragments on hover", TRUE),
            tags$p(
              class = "small text-muted mt-n2 mb-3",
              bsicons::bs_icon("exclamation-triangle-fill",
                               class = "text-warning me-1"),
              "AI-generated. Assignments are predicted from the lipid's name, ",
              "formula and adduct by matching computed fragment masses against ",
              "the peaks. They are not curated reference data and only a part ",
              "of them has been checked against the literature."
            ),
            numericInput("match_tol", "Match tolerance (Da)",
                         value = 0.01, min = 0.0001, max = 1, step = 0.005)
          )
        ),
        class = "d-flex justify-content-between align-items-center"
      ),
      plotlyOutput("spectrum"),
      card_footer(uiOutput("meta"))
    )
  )
)

# ------------------------------------------------------------ server --------

server <- function(input, output, session) {

  updateSelectizeInput(
    session, "class", server = TRUE,
    choices = setNames(CLASSES$lipid_class,
                       sprintf("%s  (%s)", CLASSES$lipid_class,
                               format(CLASSES$n, big.mark = ",", trim = TRUE)))
  )
  updateSelectizeInput(
    session, "adduct", server = TRUE,
    choices = setNames(ADDUCTS$adduct,
                       sprintf("%s  (%s)", ADDUCTS$adduct,
                               format(ADDUCTS$n, big.mark = ",", trim = TRUE)))
  )

  observeEvent(input$about, {
    showModal(modalDialog(
      title = "About ShinyLipids",
      easyClose = TRUE,
      footer = modalButton("Close"),

      tags$p(
        class = "mb-3",
        "A browser for the in-silico MS/MS spectra of the MS-DIAL lipid ",
        "library. Search by name, precursor m/z, class or adduct, inspect the ",
        "annotated fragment spectrum, and compare two records head-to-head."
      ),

      tags$h6("Spectral library"),
      tags$p(
        class = "small mb-1",
        sprintf("%s records across %d lipid classes and %d adducts, converted ",
                format(N_RECORDS, big.mark = ","), nrow(CLASSES), nrow(ADDUCTS)),
        "from the MS-DIAL ", tags$code(".lbm2"), " library file ",
        tags$code(LIBRARY_FILE), "."
      ),

      tags$h6(class = "mt-3", "Please cite"),
      tags$p(
        class = "small mb-1",
        "Tsugawa H, Ikeda K, Takahashi M, ",
        tags$em("et al."), " ",
        "A lipidome atlas in MS-DIAL 4. ",
        tags$em("Nature Biotechnology"), " 38, 1159\u20131163 (2020). ",
        tags$a(href = "https://doi.org/10.1038/s41587-020-0531-2",
               target = "_blank", rel = "noopener",
               "doi:10.1038/s41587-020-0531-2")
      ),
      tags$p(
        class = "small mb-0",
        "Tsugawa H, Cajka T, Kind T, ",
        tags$em("et al."), " ",
        "MS-DIAL: data-independent MS/MS deconvolution for comprehensive ",
        "metabolome analysis. ",
        tags$em("Nature Methods"), " 12, 523\u2013526 (2015). ",
        tags$a(href = "https://doi.org/10.1038/nmeth.3393",
               target = "_blank", rel = "noopener",
               "doi:10.1038/nmeth.3393")
      ),

      tags$p(
        class = "small text-muted mt-3 mb-0 pt-2 border-top",
        bsicons::bs_icon("exclamation-triangle-fill"),
        tags$span(class = "ms-1",
                  "Fragment annotations shown in the spectrum are ",
                  "AI-generated predictions, not curated library data.")
      )
    ))
  })

  observeEvent(input$reset, {
    updateTextInput(session, "name", value = "")
    updateNumericInput(session, "mz", value = NA)
    updateSelectizeInput(session, "class", selected = character(0))
    updateSelectizeInput(session, "adduct", selected = character(0))
    updateRadioButtons(session, "mode", selected = "")
  })

  # Debounced so typing a name does not fire a query per keystroke.
  criteria <- reactive({
    mz <- if (isTruthy(input$mz) && !is.na(input$mz)) input$mz else NULL
    # Clearing the tolerance box yields NA, which would silently turn the
    # BETWEEN into a no-match instead of an error; fall back to the default.
    raw_tol <- if (isTruthy(input$tol) && !is.na(input$tol)) input$tol else 0.01
    tol <- if (is.null(mz)) 0.01
           else if (identical(input$tol_unit, "ppm")) mz * raw_tol / 1e6
           else raw_tol

    list(name = input$name, mz = mz, tol = tol,
         class = input$class, mode = input$mode, adduct = input$adduct)
  }) |> debounce(350)

  # Fetch one extra row to detect truncation without paying for a COUNT(*).
  hits <- reactive({
    cr <- criteria()
    search_lipids(con, name = cr$name, mz = cr$mz, tol = cr$tol,
                  class = cr$class, mode = cr$mode, adduct = cr$adduct,
                  limit = MAX_ROWS + 1)
  })

  shown <- reactive(utils::head(hits(), MAX_ROWS))

  output$n_matches <- renderText({
    n <- nrow(shown())
    if (nrow(hits()) > MAX_ROWS) paste0(format(n, big.mark = ","), "+")
    else format(n, big.mark = ",")
  })
  output$match_note <- renderText({
    if (nrow(hits()) > MAX_ROWS) "showing first 500 — narrow filters"
    else if (nrow(shown()) == 0) "no lipids match"
    else "all matches shown"
  })

  output$hits <- renderDT({
    d <- shown()
    # Ion mode and peak count are omitted: the adduct string already carries
    # the polarity, and both appear in the spectrum card's footer on select.
    datatable(
      data.frame(
        Name    = d$name,
        Class   = d$lipid_class,
        `m/z`   = sprintf("%.4f", d$precursor_mz),
        Adduct  = d$adduct,
        check.names = FALSE
      ),
      selection  = "single",
      rownames   = FALSE,
      fillContainer = TRUE,
      options = list(
        paging = FALSE, scrollY = "100%", scrollCollapse = TRUE,
        scrollX = FALSE, dom = "ti", autoWidth = FALSE,
        columnDefs = list(
          list(className = "dt-right", targets = 2),
          list(width = "40%", targets = 0),
          list(width = "18%", targets = 1),
          list(width = "21%", targets = 2),
          list(width = "21%", targets = 3)
        )
      ),
      # Full name on hover, since long names are truncated to fit.
      callback = JS(
        "table.on('draw.dt', function() {",
        "  table.cells().every(function() {",
        "    var n = this.node(); n.title = n.textContent;",
        "  });",
        "});"
      )
    )
  }, server = TRUE)

  selected <- reactive({
    i <- input$hits_rows_selected
    d <- shown()
    if (is.null(i) || length(i) == 0 || nrow(d) == 0) return(NULL)
    get_spectrum(con, d$id[i])
  })

  output$sel_name <- renderText({
    s <- selected(); if (is.null(s)) "—" else s$meta$name
  })
  output$sel_sub <- renderText({
    s <- selected()
    if (is.null(s)) "click a result row"
    else sprintf("%.4f  %s", s$meta$precursor_mz, s$meta$adduct)
  })
  output$spec_title <- renderText({
    s <- selected()
    if (is.null(s)) "MS/MS spectrum" else paste("MS/MS —", s$meta$name)
  })

  # --- mirror comparison ----------------------------------------------------
  # The pinned spectrum is held by value, not by id: it must survive the user
  # searching again and the row disappearing from the result list, which is the
  # whole point of pinning it.
  reference <- reactiveVal(NULL)

  observeEvent(input$pin_ref, {
    s <- selected()
    if (!is.null(s)) reference(s)
  })
  observeEvent(input$clear_ref, reference(NULL))

  output$compare_btn <- renderUI({
    if (is.null(reference())) {
      tooltip(
        actionButton(
          "pin_ref", NULL, icon = bsicons::bs_icon("layers-half"),
          class = "btn btn-sm btn-link p-0 border-0 text-decoration-none",
          `aria-label` = "Pin as comparison spectrum",
          disabled = if (is.null(selected())) NA else NULL
        ),
        "Pin this spectrum, then pick another to mirror it against"
      )
    } else {
      tooltip(
        actionButton(
          "clear_ref", NULL, icon = bsicons::bs_icon("x-circle-fill"),
          class = "btn btn-sm btn-link p-0 border-0 text-decoration-none",
          style = paste0("color:", REF_COLOR, ";"),
          `aria-label` = "Clear comparison spectrum"
        ),
        paste0("Comparing against ", reference()$meta$name, " — click to clear")
      )
    }
  })

  # Matched peaks and score, recomputed whenever either spectrum or the
  # tolerance changes. NULL when there is nothing to compare.
  comparison <- reactive({
    s <- selected(); r <- reference()
    if (is.null(s) || is.null(r)) return(NULL)
    tol <- if (isTruthy(input$match_tol) && !is.na(input$match_tol)) {
      input$match_tol
    } else 0.01
    m <- match_peaks(s$peaks, r$peaks, tol)
    list(match = m, tol = tol,
         score = spectral_similarity(s$peaks, r$peaks, m))
  })

  # Keep the client's copy buffer in step with the selection; an empty string
  # disables the button.
  observe({
    session$sendCustomMessage("spectrumPeaks",
                              list(text = peak_list_tsv(selected())))
  })

  # Peak explanations, kept out of the plot builder so that the toggle and a
  # change of reference each redo only their own half.
  explain_for <- function(sp) {
    if (is.null(sp)) return(character(0))
    if (!isTRUE(input$explain_peaks)) return(rep("", nrow(sp$peaks)))
    annotate_peaks(sp$meta, sp$peaks)
  }
  explained     <- reactive(explain_for(selected()))
  explained_ref <- reactive(explain_for(reference()))

  output$spectrum <- renderPlotly({
    s <- selected()

    if (is.null(s)) {
      return(
        plot_ly(type = "scatter", mode = "markers") |>
          layout(
            xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
            annotations = list(
              text = "Select a lipid from the results to view its spectrum",
              showarrow = FALSE, xref = "paper", yref = "paper",
              x = 0.5, y = 0.5, font = list(size = 14, color = "#7b8a8b")
            )
          ) |> config(displayModeBar = FALSE)
      )
    }

    r  <- reference()
    cmp <- comparison()

    # The window must cover both spectra so the mirror halves stay aligned.
    mzs <- c(s$peaks$mz, s$meta$precursor_mz)
    if (!is.null(r)) mzs <- c(mzs, r$peaks$mz, r$meta$precursor_mz)
    rng <- range(mzs)
    pad <- max(diff(rng) * 0.06, 5)
    xlim <- c(rng[1] - pad, rng[2] + pad)

    # Labels are laid out before the traces because the axis range, and so the
    # precursor line's height, depends on how many lanes they need.
    lay <- spectrum_layout(s, r, xlim, input$label_n,
                           enabled = isTRUE(input$label_peaks))
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
      if (isTRUE(input$show_precursor)) {
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

    p <- plot_ly(source = SPEC_SOURCE) |>
      add_spectrum(s, 1, TOP_COLOR,
                   if (is.null(cmp)) seq_len(nrow(s$peaks)) else cmp$match$a,
                   explained())
    if (!is.null(r)) {
      p <- p |> add_spectrum(r, -1, REF_COLOR, cmp$match$b, explained_ref()) |>
        add_segments(x = xlim[1], xend = xlim[2], y = 0, yend = 0,
                     line = list(color = "#95a5a6", width = 1),
                     hoverinfo = "none", showlegend = FALSE)
    }

    # In mirror mode the axis runs negative, but intensity is a magnitude on
    # both halves, so the tick labels stay positive.
    ticks <- if (is.null(r)) c(0, 25, 50, 75, 100) else
      c(-100, -50, 0, 50, 100)

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
  })

  # Labels are laid out server-side at render, so a zoom would otherwise keep
  # the spacing computed for the whole spectrum and never resolve a crowded
  # region. Re-lay them out for the visible window on every zoom/pan, pushing
  # the result through a proxy so the user's zoom survives the update.
  observeEvent(event_data("plotly_relayout", source = SPEC_SOURCE), {
    ed <- event_data("plotly_relayout", source = SPEC_SOURCE)
    s <- selected()
    if (is.null(s) || is.null(ed)) return()

    r <- reference()
    mzs <- c(s$peaks$mz, s$meta$precursor_mz)
    if (!is.null(r)) mzs <- c(mzs, r$peaks$mz, r$meta$precursor_mz)
    full <- range(mzs)
    pad <- max(diff(full) * 0.06, 5)
    full <- c(full[1] - pad, full[2] + pad)

    # Zoom/pan sends the new bounds; double-click to reset sends autorange.
    # Anything else (our own annotation push, a resize) is not a range change
    # and must be ignored, or the proxy update below would loop.
    xlim <- if (!is.null(ed[["xaxis.range[0]"]]) && !is.null(ed[["xaxis.range[1]"]])) {
      sort(c(as.numeric(ed[["xaxis.range[0]"]]), as.numeric(ed[["xaxis.range[1]"]])))
    } else if (isTRUE(ed[["xaxis.autorange"]])) {
      full
    } else {
      return()
    }

    lay <- spectrum_layout(s, r, xlim, input$label_n,
                           enabled = isTRUE(input$label_peaks))

    upd <- list(annotations = lay$annotations)
    # Only reclaim vertical headroom if the user has not set their own y range
    # (box zoom changes both axes); overriding it would fight them.
    if (!any(grepl("^yaxis\\.range", names(ed)))) {
      upd[["yaxis.range"]] <- c(-lay$y_bot, lay$y_top)
    }
    plotlyProxyInvoke(plotlyProxy("spectrum", session), "relayout", upd)
  })

  output$meta <- renderUI({
    s <- selected()
    if (is.null(s)) return(NULL)
    m <- s$meta
    field <- function(label, value) {
      if (is.null(value) || is.na(value) || !nzchar(as.character(value))) return(NULL)
      tags$div(class = "me-4 d-inline-block",
               tags$small(class = "text-muted", label), " ",
               tags$span(as.character(value)))
    }
    cmp <- comparison()
    r <- reference()

    tags$div(
      tags$div(
        class = "small",
        field("Class", m$lipid_class),
        field("Formula", m$formula),
        field("Adduct", m$adduct),
        field("Mode", m$ion_mode),
        field("RT", if (!is.na(m$retention_time))
                      sprintf("%.2f min", m$retention_time)),
        field("CCS", if (!is.na(m$ccs) && m$ccs > 0) sprintf("%.1f Å²", m$ccs)),
        field("Peaks", m$n_peaks),
        field("InChIKey", m$inchikey)
      ),
      # The annotations are predictions, and they sit one hover away from
      # looking like curated library content. Say so where they are read,
      # not only behind the options popover.
      if (isTRUE(input$explain_peaks)) tags$div(
        class = "small mt-1 pt-1 border-top",
        style = "color:#a5680a;",
        bsicons::bs_icon("exclamation-triangle-fill"),
        tags$span(
          class = "ms-1",
          "Fragment annotations are AI-generated predictions, not curated ",
          "library data \u2014 not all of them have been verified."
        )
      ),
      if (!is.null(cmp)) tags$div(
        class = "small mt-1 pt-1 border-top",
        tags$span(style = paste0("color:", REF_COLOR, ";"),
                  bsicons::bs_icon("layers-half"), " vs ", r$meta$name),
        tags$span(class = "ms-1 text-muted",
                  sprintf("(%.4f %s)", r$meta$precursor_mz, r$meta$adduct)),
        tags$span(class = "ms-4",
                  tags$small(class = "text-muted", "matched"), " ",
                  sprintf("%d of %d / %d peaks",
                          length(cmp$match$a), nrow(s$peaks), nrow(r$peaks))),
        tags$span(class = "ms-4",
                  tags$small(class = "text-muted", "cosine"), " ",
                  if (is.na(cmp$score)) "—" else sprintf("%.3f", cmp$score)),
        tags$span(class = "ms-4 text-muted",
                  sprintf("± %g Da", cmp$tol))
      )
    )
  })
}

shinyApp(ui, server)
