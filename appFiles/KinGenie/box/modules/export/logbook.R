box::use(
  shiny[
    column,
    downloadHandler,
    downloadLink,
    fluidRow,
    HTML,
    moduleServer,
    NS,
    p,
    tagList
  ],
  shinydashboard[
    box
  ]
)

#' @export
exportLogbookUI <- function(id) {
  ns <- NS(id)
  tagList(
    box(
      title = "Logbook", width = 2, solidHeader = T, status = "primary",
      fluidRow(
        column(
          12, p(
            style = "font-size: 120%", HTML(""),
            downloadLink(ns("download_log_book"), "Logbook")
          )
        )
      )
    )
  )
}

#' @export
exportLogbookServer <- function(id, logbook) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$download_log_book <- downloadHandler(
      filename = function() {
        paste0("logbook_KinGenie_", Sys.Date(), ".txt")
      },
      content = function(file) {
        lines <- trimws(unlist(logbook$logbook()))

        lines <- paste(lines, collapse = "\n")

        cat(lines, file = file, sep = "\n")
      }
    )
  })
}
