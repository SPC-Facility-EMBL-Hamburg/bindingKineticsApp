box::use(
  .. / export / data[
    exportDataUI
  ],
  .. / export / logbook[
    exportLogbookUI
  ],
  shiny[
    fluidRow,
    NS,
    tagList
  ]
)

#' @export
menuExportUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      exportDataUI("exportData"),
      exportLogbookUI("exportLogbook")
    )
  )
}
