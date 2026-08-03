box::use(
  .. / import / load_data[
    loadDataUI
  ],
  .. / import / plot_box[
    plottingBoxUI
  ],
  .. / import / plot_input[
    plotInputUI
  ],
  .. / import / processing[
    processingUI
  ],
  shiny[
    column,
    fluidRow,
    NS,
    tagList
  ]
)

#' @export
menuImportUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        4,
        loadDataUI("load"),
        processingUI("processing")
      ),
      column(
        8,
        plotInputUI("plotInput")
      )
    ),
    fluidRow(
      column(
        12,
        plottingBoxUI("plottingBox")
      )
    )
  )
}
