box::use(
  .. / analysis / fit_controls[
    fitControlsUI
  ],
  .. / analysis / fit_data[
    fitDataUI
  ],
  .. / analysis / fit_results[
    fitResultsUI
  ],
  .. / analysis / plot_box[
    plottingBoxFitUI
  ],
  shiny[
    column,
    fluidRow,
    NS,
    tagList
  ]
)

#' @export
menuAnalyseUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        5,
        fitDataUI("fitData")
      ),
      column(
        7,
        fitControlsUI("fitControls")
      ),
      column(
        12,
        fitResultsUI("fitResults"),
        plottingBoxFitUI("plottingBoxFit")
      )
    )
  )
}
