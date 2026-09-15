## ShinyLipids -- browse and visualise MS/MS spectra from an MS-DIAL lipid library.
##
## Two tabs: a library browser (search by name / m/z / class / adduct) and a
## spectrum search that matches MS/MS spectra -- one pasted peak list, or a
## whole MSP file of them -- against every record.
##
## The library lives in DB/lipids.sqlite, built from the .lbm2 file by
## scripts/lbm2_to_sqlite.py. See R/db.R for the record query layer,
## R/speclib.R for the whole-library spectral search, and R/spectrum.R for the
## plot and peak-label layout.

library(shiny)
library(bslib)
library(plotly)
library(DT)

source("R/db.R")
source("R/spectrum.R")
source("R/fragments.R")
source("R/speclib.R")

con <- lipid_db()
onStop(function() DBI::dbDisconnect(con))

CLASSES   <- lipid_classes(con)
ADDUCTS   <- lipid_adducts(con)
N_RECORDS <- lipid_count(con)
N_PEAKS   <- lipid_peak_count(con)

MAX_ROWS <- 500

# Most intense peaks kept from a query spectrum. Real centroided MS/MS rarely
# reaches this; a profile-mode record would otherwise expand into tens of
# millions of candidate pairs and stall the session.
MAX_QUERY_PEAKS <- 500

# Spectra taken from one uploaded MSP file. At ~10 ms each this is ~20 s of
# searching, about as long as a single-threaded Shiny process should be held
# by one request.
MAX_QUERY_SPECTRA <- 2000

# An MSP file for a full run reaches a few MB; Shiny's own default is 5 MB.
options(shiny.maxRequestSize = 64 * 1024^2)

# A sample upload, served from the app directory rather than www/: it is an
# example input, not an asset the page itself loads.
EXAMPLE_MSP <- "example.msp"

APP_TITLE <- "ShinyLipids — MS/MS spectrum browser - v0.4.1"
LIBRARY_FILE <- "Msp20251120132005_NCDK_conventional_converted_dev.lbm2"


# ---------------------------------------------------------------- UI --------

ui <- page_navbar(
  title = APP_TITLE,
  window_title = APP_TITLE,
  theme = bs_theme(version = 5, preset = "flatly"),
  # page_navbar takes no page-level `class`, so the "bslib-page-dashboard"
  # styling page_sidebar used to carry is gone; it only softened the navbar
  # border, which the tab strip now occupies anyway.
  fillable = TRUE,

  # Both tabs share the table styling and the clipboard helper, so they ride
  # along in the navbar header rather than being repeated per panel.
  header = tagList(
    # Keep the results table inside its card: fixed layout plus ellipsis on the
    # name column, rather than letting DT overflow into a horizontal scrollbar.
    tags$head(tags$style(HTML("
      #hits table, #q_hits table {
        table-layout: fixed; width: 100% !important;
      }
      #hits td, #hits th, #q_hits td, #q_hits th {
        overflow: hidden; text-overflow: ellipsis; white-space: nowrap;
      }
      #hits tbody tr, #q_hits tbody tr { cursor: pointer; }
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
  ),

  nav_panel(
    title = "Library search",
    icon = bsicons::bs_icon("card-list"),
    layout_sidebar(
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
  ),

  # ----------------------------------------------- spectrum search tab ------
  nav_panel(
    title = "Spectrum search",
    icon = bsicons::bs_icon("graph-up-arrow"),
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        title = "Query spectra",

        # Re-rendered rather than updated: Shiny has no way to clear a
        # fileInput's chosen file, so "Clear" replaces the control itself.
        uiOutput("q_file_ui"),

        # What the file has to look like, said where the file is chosen. The
        # sample is only offered when it is actually on disk, so a deployment
        # without it gets no dead link.
        tags$div(
          class = "small text-muted mb-3",
          "Per record: ", tags$code("PRECURSORMZ"), ", ",
          tags$code("IONMODE"), ", optionally ", tags$code("PRECURSORTYPE"),
          ", then the peaks as tab-separated m/z and intensity.",
          if (file.exists(EXAMPLE_MSP)) tagList(
            " ", downloadLink("q_example", "Download an example"), "."
          )
        ),

        checkboxInput(
          "q_file_meta",
          "Use each record's precursor, polarity and adduct", TRUE),
        tags$div(
          class = "small text-muted mb-3",
          "Off, or where a record does not declare one, the settings below ",
          "apply to every spectrum."
        ),

        tags$hr(class = "my-2"),
        tags$label(class = "control-label mb-1", "…or paste one spectrum"),
        textAreaInput(
          "q_text", NULL, rows = 6, resize = "vertical",
          placeholder = paste(
            "283.2643, 100", "303.2330, 55.3", "766.5392, 12",
            "", "Comma, tab, semicolon or spaces all work.",
            sep = "\n")
        ),
        numericInput("q_prec", "Precursor m/z (pasted spectrum)", value = NA,
                     min = 0, step = 0.0001),

        tags$hr(class = "my-2"),

        # Both tolerances apply at once and the wider wins, so "10 ppm but
        # never tighter than 0.005 Da" is a single setting rather than a mode
        # switch. Either box at 0 turns that half off.
        tags$label(class = "control-label mb-1", "Mass tolerance"),
        layout_columns(
          col_widths = c(6, 6), gap = "0.5rem",
          numericInput("q_tol_da", "± Da", value = 0.01, min = 0,
                       step = 0.001),
          numericInput("q_tol_ppm", "± ppm", value = 0, min = 0, step = 1)
        ),
        tags$div(
          class = "small text-muted mb-3",
          "The wider of the two is used at each m/z; 0 disables one."
        ),

        radioButtons("q_mode", "Polarity",
                     choices = c("Any" = "", "Positive" = "Positive",
                                 "Negative" = "Negative"),
                     selected = "", inline = TRUE),

        selectizeInput(
          "q_adduct", "Ion (adduct)", choices = NULL, multiple = TRUE,
          options = list(placeholder = "All adducts",
                         plugins = list("remove_button"))
        ),

        layout_columns(
          col_widths = c(6, 6), gap = "0.5rem",
          numericInput("q_min_match", "Min. matches", value = 2, min = 1,
                       step = 1),
          selectInput("q_rank", "Rank by",
                      choices = c("Weighted dot" = "wdot", "Dot" = "dot",
                                  "Reverse dot" = "rdot"),
                      selected = "wdot")
        ),

        actionButton("q_run", "Search library", class = "btn-primary",
                     icon = bsicons::bs_icon("search")),
        actionButton("q_clear", "Clear", class = "btn-outline-secondary btn-sm",
                     icon = icon("rotate-left"))
      ),

      layout_column_wrap(
        width = 1/3, fill = FALSE, heights_equal = "row",
        value_box(
          title = "Query spectra", value = textOutput("q_n_queries", inline = TRUE),
          showcase = bsicons::bs_icon("list-ol"), theme = "primary",
          textOutput("q_peak_note", inline = TRUE)
        ),
        value_box(
          title = "Candidates searched",
          value = textOutput("q_n_cand", inline = TRUE),
          showcase = bsicons::bs_icon("funnel"), theme = "secondary",
          textOutput("q_cand_note", inline = TRUE)
        ),
        value_box(
          title = "Best match",
          value = tags$span(style = "font-size: 1.2rem; line-height: 1.3;",
                            textOutput("q_best", inline = TRUE)),
          showcase = bsicons::bs_icon("trophy"), theme = "success",
          textOutput("q_best_sub", inline = TRUE)
        )
      ),

      # Only earns its space once there is more than one spectrum to pick from;
      # a one-row table above a pasted spectrum's own results is just noise.
      conditionalPanel(
        condition = "output.q_is_batch === true",
        card(
          # The card sits in the panel's flex column, where the results row
          # below grows and would otherwise squeeze this one down to two
          # visible rows -- no use as the batch's navigation. Pinning it at
          # its height keeps five or six.
          full_screen = TRUE, height = 320,
          class = "flex-grow-0 flex-shrink-0",
          card_header(
            "Uploaded spectra",
            tooltip(
              bsicons::bs_icon("info-circle", title = "About the query list"),
              paste(
                "One row per record in the MSP file; click a row to see its",
                "hits and mirror plot. Best match and score are for the",
                "spectrum's own top hit under the current ranking.",
                "An adduct the library does not contain is shown but not used",
                "as a filter, so a typo cannot silently return nothing.")
            ),
            class = "d-flex justify-content-between align-items-center"
          ),
          DTOutput("q_queries")
        )
      ),

      layout_columns(
        col_widths = breakpoints(sm = 12, lg = c(8, 4)),

        card(
          full_screen = TRUE, min_height = 340,
          card_header(
            textOutput("q_hits_title", inline = TRUE),
            tags$div(
              class = "d-flex align-items-center gap-3",
              downloadLink("q_download", bsicons::bs_icon(
                "download", title = "Download all hits as CSV")),
              popover(
                bsicons::bs_icon("gear", title = "Scoring options"),
                title = "Scoring",
                numericInput("q_top_n", "Hits returned per spectrum",
                             value = 200, min = 10, max = 2000, step = 10),
                tags$hr(),
                tags$p(class = "small text-muted mb-2",
                       "Weighted dot product uses w = m/z", tags$sup("a"),
                       " × I", tags$sup("b"), ". Stein & Scott's a = 3, ",
                       "b = 0.6 is the NIST/MS-DIAL default; a = 2, b = 0.5 ",
                       "is the MassBank convention."),
                layout_columns(
                  col_widths = c(6, 6), gap = "0.5rem",
                  numericInput("q_mz_pow", "m/z power (a)", value = MZ_POWER,
                               min = 0, max = 5, step = 0.1),
                  numericInput("q_int_pow", "Intensity power (b)",
                               value = INT_POWER, min = 0, max = 2, step = 0.1)
                )
              ),
              tooltip(
                bsicons::bs_icon("info-circle", title = "About the scores"),
                paste(
                  "All three scores are cosine similarities on 0-1.",
                  "Dot: on relative intensities; unmatched peaks on either",
                  "side count against the score.",
                  "W.dot (weighted): the same on m/z- and intensity-weighted",
                  "peaks, which favours specific high-mass fragments.",
                  "R.dot (reverse): query peaks absent from the library record",
                  "are ignored, so co-isolated contamination does not penalise",
                  "a hit - but records with very few peaks score high easily.",
                  "Match is matched peaks out of the query's peaks.")
              )
            ),
            class = "d-flex justify-content-between align-items-center"
          ),
          DTOutput("q_hits")
        ),

        card(
          full_screen = TRUE, min_height = 340,
          card_header(
            textOutput("q_spec_title", inline = TRUE),
            popover(
              bsicons::bs_icon("gear", title = "Plot options"),
              title = "Plot options",
              checkboxInput("q_label_peaks", "Label peaks", TRUE),
              sliderInput("q_label_n", "Max labels", min = 1, max = 20,
                          value = 8),
              checkboxInput("q_show_precursor", "Mark precursor m/z", TRUE)
            ),
            class = "d-flex justify-content-between align-items-center"
          ),
          plotlyOutput("q_spectrum"),
          card_footer(uiOutput("q_meta"))
        )
      )
    )
  ),

  nav_spacer(),
  nav_item(
    actionLink("about", "About", class = "nav-link",
               icon = bsicons::bs_icon("info-circle"))
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

      tags$h6("Spectrum search"),
      tags$p(
        class = "small mb-3",
        "The second tab matches MS/MS spectra against every record the ",
        "polarity, adduct and precursor filters allow \u2014 one pasted peak ",
        "list, or a whole ", tags$code(".msp"), " file of them, in which case ",
        "each record's own ", tags$code("PRECURSORMZ"), ", ",
        tags$code("IONMODE"), " and ", tags$code("PRECURSORTYPE"), " drive its ",
        "search. Hits are scored by ",
        "dot product, weighted dot product (w = m/z", tags$sup("a"), " \u00d7 I",
        tags$sup("b"), ", after Stein & Scott) and reverse dot product, all as ",
        "cosine similarities on 0\u20131. The first search of a session spends ",
        "about ten seconds building an in-memory index of all ",
        format(N_PEAKS, big.mark = ","), " library peaks; every search after ",
        "that is near-instant."
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
      return(spectrum_placeholder(
        "Select a lipid from the results to view its spectrum"))
    }
    cmp <- comparison()
    spectrum_plot(
      top = s, ref = reference(),
      match_a = if (is.null(cmp)) NULL else cmp$match$a,
      match_b = if (is.null(cmp)) integer(0) else cmp$match$b,
      ann_top = explained(), ann_ref = explained_ref(),
      label_peaks = isTRUE(input$label_peaks), label_n = input$label_n,
      show_precursor = isTRUE(input$show_precursor), source = SPEC_SOURCE
    )
  })

  # Labels are laid out server-side at render, so a zoom would otherwise keep
  # the spacing computed for the whole spectrum and never resolve a crowded
  # region. Re-lay them out for the visible window on every zoom/pan, pushing
  # the result through a proxy so the user's zoom survives the update.
  observeEvent(event_data("plotly_relayout", source = SPEC_SOURCE), {
    s <- selected()
    if (is.null(s)) return()
    upd <- spectrum_relayout(event_data("plotly_relayout", source = SPEC_SOURCE),
                             s, reference(), input$label_n,
                             isTRUE(input$label_peaks))
    if (is.null(upd)) return()
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

  # ------------------------------------------------- spectrum search --------

  updateSelectizeInput(
    session, "q_adduct", server = TRUE,
    choices = setNames(ADDUCTS$adduct,
                       sprintf("%s  (%s)", ADDUCTS$adduct,
                               format(ADDUCTS$n, big.mark = ",", trim = TRUE)))
  )

  # --- query spectra: an uploaded MSP file, or one pasted peak list ---------

  q_result   <- reactiveVal(NULL)   # the last search, or NULL
  q_upload   <- reactiveVal(NULL)   # the parsed MSP file, or NULL
  file_epoch <- reactiveVal(0L)

  # Shiny cannot clear a fileInput's selection, so the control is rebuilt.
  output$q_file_ui <- renderUI({
    file_epoch()
    fileInput("q_file", "MSP file", accept = c(".msp", ".txt", "text/plain"),
              buttonLabel = "Browse…", placeholder = "no file selected")
  })

  observeEvent(input$q_file, {
    f <- input$q_file
    if (is.null(f)) return()
    parsed <- tryCatch(
      parse_msp(readLines(f$datapath, warn = FALSE, encoding = "UTF-8")),
      error = function(e) NULL)

    if (is.null(parsed) || nrow(parsed$meta) == 0) {
      showNotification(
        paste0("No MS/MS spectra found in '", f$name,
               "'. Expected MSP: header lines such as PRECURSORMZ and ",
               "IONMODE, then tab-separated m/z and intensity."),
        type = "error", duration = 10)
      # Clear the control too, or it keeps reading "Upload complete" under a
      # file the app just refused.
      file_epoch(file_epoch() + 1L)
      q_upload(NULL); q_result(NULL)
      return()
    }
    if (nrow(parsed$meta) > MAX_QUERY_SPECTRA) {
      showNotification(
        sprintf("'%s' holds %s spectra; keeping the first %s.", f$name,
                format(nrow(parsed$meta), big.mark = ","),
                format(MAX_QUERY_SPECTRA, big.mark = ",")),
        type = "warning", duration = 10)
      keep <- seq_len(MAX_QUERY_SPECTRA)
      parsed$meta <- parsed$meta[keep, , drop = FALSE]
      parsed$peaks <- parsed$peaks[keep]
    }
    if (parsed$n_dropped > 0) {
      showNotification(
        sprintf("%d record%s in '%s' carried no peaks and %s skipped.",
                parsed$n_dropped, if (parsed$n_dropped == 1) "" else "s",
                f$name, if (parsed$n_dropped == 1) "was" else "were"),
        type = "message")
    }
    parsed$file <- f$name
    q_upload(parsed)
    q_result(NULL)
  })

  observeEvent(input$q_clear, {
    updateTextAreaInput(session, "q_text", value = "")
    updateNumericInput(session, "q_prec", value = NA)
    updateSelectizeInput(session, "q_adduct", selected = character(0))
    updateRadioButtons(session, "q_mode", selected = "")
    file_epoch(file_epoch() + 1L)
    q_upload(NULL)
    q_result(NULL)
  })

  # Parsed as the user types so the peak count is live feedback on the paste.
  q_parsed <- reactive(parse_peak_text(input$q_text)) |> debounce(300)

  # A numeric input the user has cleared reads back as NA; every one of these
  # would otherwise turn a search into a silent no-match.
  num_or <- function(x, default) {
    if (isTruthy(x) && !is.na(x) && is.finite(x)) x else default
  }
  q_opts <- reactive(list(
    tol_da    = num_or(input$q_tol_da, 0.01),
    tol_ppm   = num_or(input$q_tol_ppm, 0),
    mz_power  = num_or(input$q_mz_pow, MZ_POWER),
    int_power = num_or(input$q_int_pow, INT_POWER),
    min_match = as.integer(num_or(input$q_min_match, 1)),
    top_n     = as.integer(num_or(input$q_top_n, 200)),
    prec      = if (isTruthy(input$q_prec) && !is.na(input$q_prec))
                  input$q_prec else NULL
  ))

  # The two input routes collapse into one shape here, so everything
  # downstream -- the query table, the scoring loop, the plot -- is written
  # once against a list of spectra rather than twice.
  q_input <- reactive({
    up <- q_upload()
    if (!is.null(up)) return(up)
    pk <- q_parsed()$peaks
    if (nrow(pk) == 0) return(NULL)
    list(
      meta = data.frame(
        name = "Pasted spectrum",
        precursor_mz = if (is.null(q_opts()$prec)) NA_real_ else q_opts()$prec,
        adduct = "", ion_mode = "", retention_time = NA_real_,
        n_peaks = nrow(pk), stringsAsFactors = FALSE),
      peaks = list(pk), n_dropped = 0L, file = NULL)
  })

  output$q_is_batch <- reactive({
    q <- q_input()
    !is.null(q) && nrow(q$meta) > 1
  })
  outputOptions(output, "q_is_batch", suspendWhenHidden = FALSE)

  # --- the search -----------------------------------------------------------

  observeEvent(input$q_run, {
    qi <- q_input()
    if (is.null(qi)) {
      showNotification("Upload an MSP file or paste a peak list first.",
                       type = "warning")
      q_result(NULL)
      return()
    }
    o <- q_opts()
    use_meta <- isTRUE(input$q_file_meta) && !is.null(q_upload())
    n <- nrow(qi$meta)

    # Peaks are capped per spectrum, not across the file: a profile-mode record
    # would otherwise expand into tens of millions of candidate pairs.
    peaks <- lapply(qi$peaks, function(pk) {
      if (nrow(pk) <= MAX_QUERY_PEAKS) return(pk)
      pk <- pk[order(-pk$rel)[seq_len(MAX_QUERY_PEAKS)], , drop = FALSE]
      pk[order(pk$mz), , drop = FALSE]
    })
    n_capped <- sum(vapply(qi$peaks, nrow, 1L) > MAX_QUERY_PEAKS)
    if (n_capped > 0) {
      showNotification(
        sprintf("%d spectr%s searched on %d most intense peaks.", n_capped,
                if (n_capped == 1) "um" else "a", MAX_QUERY_PEAKS),
        type = "message")
    }

    ix <- if (peak_index_ready()) peak_index(con) else {
      withProgress(
        message = "Indexing the library", value = 0,
        detail = "one-off, about ten seconds",
        peak_index(con, function(f, msg) setProgress(value = f, detail = msg))
      )
    }

    t0 <- Sys.time()
    hits <- withProgress(
      message = if (n > 1) sprintf("Searching %s spectra", format(n, big.mark = ","))
                else "Searching the library",
      value = 0,
      lapply(seq_len(n), function(i) {
        m <- qi$meta[i, ]
        # A record's own header wins where it has one; the sidebar covers the
        # rest. An adduct the library has never heard of is dropped rather
        # than applied, or one stray PRECURSORTYPE would return nothing at all
        # with no hint as to why.
        mode <- if (use_meta && nzchar(m$ion_mode)) m$ion_mode else input$q_mode
        add  <- if (use_meta && nzchar(m$adduct) && m$adduct %in% ADDUCTS$adduct) {
          m$adduct
        } else input$q_adduct
        prec <- if (use_meta) {
          if (is.finite(m$precursor_mz)) m$precursor_mz else NULL
        } else o$prec

        h <- search_spectrum(
          ix, peaks[[i]], tol_da = o$tol_da, tol_ppm = o$tol_ppm,
          mode = mode, adduct = add, precursor = prec,
          min_match = o$min_match, mz_power = o$mz_power,
          int_power = o$int_power, rank_by = input$q_rank, top_n = o$top_n)

        if (n > 1) setProgress(value = i / n, detail = m$name)
        attr(h, "prec") <- prec
        h
      })
    )
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")

    # Names and classes are not in the index -- it holds only what scoring
    # needs -- so they are fetched once for every id actually on show.
    ids <- unique(unlist(lapply(hits, `[[`, "id")))
    meta <- if (length(ids)) DBI::dbGetQuery(
      con, paste0("SELECT id, name, lipid_class, formula, retention_time
                     FROM lipid WHERE id IN (",
                  paste(ids, collapse = ","), ")")) else NULL
    hits <- lapply(hits, function(h) {
      if (nrow(h) == 0) {
        h$name <- character(0); h$lipid_class <- character(0)
        h$formula <- character(0); h$retention_time <- numeric(0)
        return(h)
      }
      k <- match(h$id, meta$id)
      h$name <- meta$name[k]; h$lipid_class <- meta$lipid_class[k]
      h$formula <- meta$formula[k]; h$retention_time <- meta$retention_time[k]
      h
    })

    q_result(list(meta = qi$meta, peaks = peaks, hits = hits,
                  elapsed = elapsed, opts = o, batch = n > 1,
                  file = qi$file))
  })

  # --- which query spectrum is on show -------------------------------------

  q_i <- reactive({
    qi <- q_input()
    if (is.null(qi)) return(NA_integer_)
    n <- nrow(qi$meta)
    if (n == 1L) return(1L)
    i <- input$q_queries_rows_selected
    if (is.null(i) || length(i) == 0 || i > n) 1L else as.integer(i)
  })

  # What is actually on screen: the spectra the last search ran on, or -- before
  # any search -- whatever is loaded. Keeping the peak note, the plot and the
  # scores on one source stops them describing different spectra after an edit.
  q_shown <- reactive({
    r <- q_result()
    if (!is.null(r)) list(meta = r$meta, peaks = r$peaks) else q_input()
  })

  q_hit_table <- reactive({
    r <- q_result(); i <- q_i()
    if (is.null(r) || is.na(i) || i > length(r$hits)) return(NULL)
    r$hits[[i]]
  })

  # --- value boxes ----------------------------------------------------------

  output$q_n_queries <- renderText({
    qi <- q_input()
    if (is.null(qi)) "0" else format(nrow(qi$meta), big.mark = ",")
  })
  output$q_peak_note <- renderText({
    qs <- q_shown(); i <- q_i()
    if (is.null(qs) || is.na(i) || i > length(qs$peaks)) {
      p <- q_parsed()
      if (p$n_skipped > 0) "no numeric peaks found" else "upload or paste a spectrum"
    } else {
      pk <- qs$peaks[[i]]
      sprintf("%d peaks, m/z %.4f – %.4f in the selected one",
              nrow(pk), min(pk$mz), max(pk$mz))
    }
  })
  output$q_n_cand <- renderText({
    h <- q_hit_table()
    if (is.null(h) || is.na(attr(h, "n_filtered"))) "—"
    else format(attr(h, "n_filtered"), big.mark = ",", scientific = FALSE)
  })
  output$q_cand_note <- renderText({
    r <- q_result(); h <- q_hit_table()
    if (is.null(h)) "not searched yet"
    else sprintf("%s shared a peak — %.2f s for %s spectr%s",
                 format(attr(h, "n_hit"), big.mark = ",", scientific = FALSE),
                 r$elapsed, format(length(r$hits), big.mark = ","),
                 if (length(r$hits) == 1) "um" else "a")
  })
  output$q_best <- renderText({
    h <- q_hit_table()
    if (is.null(h) || nrow(h) == 0) "—" else h$name[1]
  })
  output$q_best_sub <- renderText({
    h <- q_hit_table()
    if (is.null(h)) "run a search"
    else if (nrow(h) == 0) "no record matched"
    else sprintf("dot %.3f · weighted %.3f · reverse %.3f",
                 h$dot[1], h$wdot[1], h$rdot[1])
  })

  # --- the uploaded-spectra table ------------------------------------------

  output$q_queries <- renderDT({
    qi <- q_input()
    if (is.null(qi)) return(NULL)
    r <- q_result()
    # Only pair the hits with the query list they were produced from: a new
    # upload must not inherit the previous file's best matches.
    have <- !is.null(r) && length(r$hits) == nrow(qi$meta) &&
            identical(r$meta, qi$meta)

    # `i` is the query row; `empty` both types the result and fills the rows
    # that have no hits yet.
    top_hit <- function(f, empty) vapply(seq_len(nrow(qi$meta)), function(i) {
      if (!have) return(empty)
      h <- r$hits[[i]]
      if (nrow(h) == 0) empty else f(h, i)
    }, empty)

    d <- data.frame(
      `#`       = seq_len(nrow(qi$meta)),
      Name      = qi$meta$name,
      `m/z`     = qi$meta$precursor_mz,
      Polarity  = ifelse(nzchar(qi$meta$ion_mode), qi$meta$ion_mode, "—"),
      Adduct    = ifelse(nzchar(qi$meta$adduct), qi$meta$adduct, "—"),
      Peaks     = qi$meta$n_peaks,
      `Best match` = top_hit(function(h, i) h$name[1], NA_character_),
      Score     = top_hit(function(h, i) h[[input$q_rank]][1], NA_real_),
      # Counted against each spectrum's own peak total, which is the capped
      # one the search actually used.
      Match     = top_hit(function(h, i) sprintf("%d / %d", h$n_match[1],
                                                 nrow(r$peaks[[i]])),
                          NA_character_),
      check.names = FALSE
    )

    datatable(
      d, selection = list(mode = "single", selected = 1),
      rownames = FALSE, fillContainer = TRUE,
      options = list(
        paging = FALSE, scrollY = "100%", scrollCollapse = TRUE,
        scrollX = FALSE, dom = "ti", autoWidth = FALSE,
        columnDefs = list(
          list(className = "dt-right",
               targets = which(names(d) %in%
                 c("#", "m/z", "Peaks", "Score", "Match")) - 1L),
          list(width = "4%",  targets = 0),
          list(width = "22%", targets = 1),
          list(width = "10%", targets = 2),
          list(width = "9%",  targets = 3),
          list(width = "11%", targets = 4),
          list(width = "7%",  targets = 5),
          list(width = "22%", targets = 6),
          list(width = "8%",  targets = 7),
          list(width = "7%",  targets = 8))
      ),
      callback = JS(
        "table.on('draw.dt', function() {",
        "  table.cells().every(function() {",
        "    var n = this.node(); n.title = n.textContent;",
        "  });",
        "});"
      )
    ) |>
      formatRound("m/z", 4, mark = "") |>
      formatRound("Score", 3)
  }, server = TRUE)

  # --- the hit table --------------------------------------------------------

  output$q_hits_title <- renderText({
    r <- q_result()
    if (is.null(r) || !isTRUE(r$batch)) "Library matches"
    else paste("Library matches —", r$meta$name[q_i()])
  })

  output$q_hits <- renderDT({
    d <- q_hit_table()
    prec <- if (is.null(d)) NULL else attr(d, "prec")
    n_query <- {
      r <- q_result()
      if (is.null(r)) 0L else nrow(r$peaks[[q_i()]])
    }
    if (is.null(d) || nrow(d) == 0) {
      d <- data.frame(Name = character(0), `m/z` = numeric(0),
                      Adduct = character(0), Match = character(0),
                      Dot = numeric(0), `W.dot` = numeric(0),
                      `R.dot` = numeric(0), check.names = FALSE)
    } else {
      tbl <- data.frame(
        Name    = d$name,
        `m/z`   = d$precursor_mz,
        Adduct  = d$adduct,
        Match   = sprintf("%d / %d", d$n_match, n_query),
        Dot     = d$dot,
        `W.dot` = d$wdot,
        `R.dot` = d$rdot,
        check.names = FALSE
      )
      # The precursor error only means anything when a precursor was given.
      if (!is.null(prec)) {
        tbl <- cbind(tbl[1:2],
                     `Δppm` = 1e6 * (d$precursor_mz - prec) / prec,
                     tbl[3:7])
      }
      d <- tbl
    }
    # table-layout is fixed (see the CSS above), so every column needs an
    # explicit share or the three score columns get squeezed off the card.
    widths <- c(Name = 23, `m/z` = 13, `Δppm` = 11, Adduct = 16,
                Match = 12, Dot = 9, `W.dot` = 10, `R.dot` = 10)[names(d)]
    widths <- 100 * widths / sum(widths)
    right <- c("m/z", "Δppm", "Match", "Dot", "W.dot", "R.dot")

    tab <- datatable(
      d,
      # Preselecting the top hit means the mirror plot is populated the moment
      # a search finishes, rather than showing an empty card until a click.
      selection = list(mode = "single", selected = if (nrow(d) > 0) 1 else NULL),
      rownames = FALSE, fillContainer = TRUE,
      options = list(
        paging = FALSE, scrollY = "100%", scrollCollapse = TRUE,
        scrollX = FALSE, dom = "ti", autoWidth = FALSE,
        columnDefs = c(
          list(list(className = "dt-right",
                    targets = which(names(d) %in% right) - 1L)),
          lapply(seq_along(widths), function(i)
            list(width = paste0(round(widths[i], 1), "%"), targets = i - 1L))
        )
      ),
      callback = JS(
        "table.on('draw.dt', function() {",
        "  table.cells().every(function() {",
        "    var n = this.node(); n.title = n.textContent;",
        "  });",
        "});"
      )
    )
    # formatRound() errors on an empty column set, which is exactly what the
    # "no precursor given" and "no hits" cases hand it.
    round_cols <- function(t, cols, digits, mark = ",") {
      cols <- intersect(cols, names(d))
      if (length(cols)) formatRound(t, cols, digits, mark = mark) else t
    }
    tab <- round_cols(tab, "m/z", 4, mark = "")
    tab <- round_cols(tab, "Δppm", 1)
    round_cols(tab, c("Dot", "W.dot", "R.dot"), 3)
  }, server = TRUE)

  # --- the mirror plot ------------------------------------------------------

  q_selected <- reactive({
    d <- q_hit_table()
    i <- input$q_hits_rows_selected
    if (is.null(d) || nrow(d) == 0 || is.null(i) || length(i) == 0) return(NULL)
    get_spectrum(con, d$id[i])
  })

  # Falls back to whatever is loaded but not yet searched, so an upload or a
  # paste is plotted straight away -- the quickest way to see that it parsed
  # the way the user meant it to.
  q_top <- reactive({
    src <- q_shown(); i <- q_i()
    if (is.null(src) || is.na(i) || i > length(src$peaks)) return(NULL)
    pk <- src$peaks[[i]]
    if (nrow(pk) == 0) return(NULL)
    m <- src$meta[i, ]
    as_query_spectrum(pk, m$precursor_mz, m$name)
  })

  q_pair <- reactive({
    top <- q_top(); lib <- q_selected()
    if (is.null(top) || is.null(lib)) return(NULL)
    o <- q_result()$opts
    score_pair(top$peaks, lib$peaks, o$tol_da, o$tol_ppm,
               o$mz_power, o$int_power)
  })

  output$q_spec_title <- renderText({
    s <- q_selected()
    if (is.null(s)) "Query spectrum" else paste("Query vs", s$meta$name)
  })

  output$q_spectrum <- renderPlotly({
    top <- q_top()
    if (is.null(top)) {
      return(spectrum_placeholder(
        "Upload an MSP file or paste an MS/MS peak list"))
    }
    m <- q_pair()
    spectrum_plot(
      top = top, ref = q_selected(),
      match_a = if (is.null(m)) NULL else m$match$a,
      match_b = if (is.null(m)) integer(0) else m$match$b,
      label_peaks = isTRUE(input$q_label_peaks), label_n = input$q_label_n,
      show_precursor = isTRUE(input$q_show_precursor), source = QSPEC_SOURCE
    )
  })

  observeEvent(event_data("plotly_relayout", source = QSPEC_SOURCE), {
    top <- q_top()
    if (is.null(top)) return()
    upd <- spectrum_relayout(
      event_data("plotly_relayout", source = QSPEC_SOURCE),
      top, q_selected(), input$q_label_n, isTRUE(input$q_label_peaks))
    if (is.null(upd)) return()
    plotlyProxyInvoke(plotlyProxy("q_spectrum", session), "relayout", upd)
  })

  output$q_meta <- renderUI({
    s <- q_selected(); m <- q_pair(); top <- q_top()
    if (is.null(s) || is.null(m) || is.null(top)) return(NULL)
    o <- q_result()$opts
    field <- function(label, value) {
      if (is.null(value) || is.na(value) || !nzchar(as.character(value))) return(NULL)
      tags$div(class = "me-4 d-inline-block",
               tags$small(class = "text-muted", label), " ",
               tags$span(as.character(value)))
    }
    tags$div(
      tags$div(
        class = "small",
        field("Class", s$meta$lipid_class),
        field("Formula", s$meta$formula),
        field("Adduct", s$meta$adduct),
        field("Mode", s$meta$ion_mode),
        field("RT", if (!is.na(s$meta$retention_time))
                      sprintf("%.2f min", s$meta$retention_time)),
        field("Peaks", s$meta$n_peaks)
      ),
      tags$div(
        class = "small mt-1 pt-1 border-top",
        tags$span(class = "me-4",
                  tags$small(class = "text-muted", "matched"), " ",
                  sprintf("%d of %d / %d peaks", length(m$match$a),
                          nrow(top$peaks), nrow(s$peaks))),
        tags$span(class = "me-4",
                  tags$small(class = "text-muted", "dot"), " ",
                  sprintf("%.3f", m$dot)),
        tags$span(class = "me-4",
                  tags$small(class = "text-muted", "weighted"), " ",
                  sprintf("%.3f", m$wdot)),
        tags$span(class = "me-4",
                  tags$small(class = "text-muted", "reverse"), " ",
                  sprintf("%.3f", m$rdot)),
        tags$span(class = "text-muted",
                  sprintf("± %g Da / %g ppm", o$tol_da, o$tol_ppm))
      )
    )
  })

  output$q_example <- downloadHandler(
    filename = function() basename(EXAMPLE_MSP),
    contentType = "text/plain",
    content = function(file) file.copy(EXAMPLE_MSP, file, overwrite = TRUE)
  )

  # --- CSV export -----------------------------------------------------------

  output$q_download <- downloadHandler(
    filename = function() sprintf("spectrum-search-%s.csv",
                                  format(Sys.time(), "%Y%m%d-%H%M%S")),
    content = function(file) {
      r <- q_result()
      if (is.null(r)) {
        write.csv(data.frame(), file, row.names = FALSE)
        return()
      }
      # Every spectrum's hits in one frame, each row carrying the query it came
      # from and its rank, so a batch run drops straight into a pivot table.
      out <- do.call(rbind, lapply(seq_along(r$hits), function(i) {
        d <- r$hits[[i]]
        if (nrow(d) == 0) return(NULL)
        prec <- attr(d, "prec")
        o <- data.frame(
          query          = r$meta$name[i],
          query_index    = i,
          query_mz       = r$meta$precursor_mz[i],
          query_adduct   = r$meta$adduct[i],
          query_polarity = r$meta$ion_mode[i],
          n_query_peaks  = nrow(r$peaks[[i]]),
          rank           = seq_len(nrow(d)),
          name           = d$name,
          lipid_class    = d$lipid_class,
          formula        = d$formula,
          precursor_mz   = round(d$precursor_mz, 4),
          adduct         = d$adduct,
          ion_mode       = d$ion_mode,
          retention_time = round(d$retention_time, 3),
          n_matched      = d$n_match,
          dot            = round(d$dot, 5),
          weighted_dot   = round(d$wdot, 5),
          reverse_dot    = round(d$rdot, 5),
          stringsAsFactors = FALSE
        )
        if (!is.null(prec)) {
          o$delta_ppm <- round(1e6 * (d$precursor_mz - prec) / prec, 2)
        }
        o
      }))
      if (is.null(out)) out <- data.frame()
      write.csv(out, file, row.names = FALSE, na = "")
    }
  )

}

shinyApp(ui, server)
