box::use(
  .. / .. / helpers[
    pandas_to_r
  ],
  shiny[
    column,
    downloadHandler,
    downloadLink,
    fluidRow,
    HTML,
    moduleServer,
    NS,
    p,
    req,
    tagList
  ],
  shinydashboard[
    box
  ],
  utils[
    write.csv
  ]
)

#' @export
exportDataUI <- function(id) {
  ns <- NS(id)
  tagList(
    box(
      title = "Data", width = 2, solidHeader = T, status = "primary",
      fluidRow(
        column(
          12,
          p(
            style = "font-size: 120%", HTML(""),
            downloadLink(ns("download_curves"), "Raw data (kinetics)")
          )
        ),
        column(
          12,
          p(
            style = "font-size: 120%", HTML(""),
            downloadLink(ns("download_fitted_curves"), "Fitted data (kinetics)")
          )
        ),
      )
    )
  )
}

#' @export
exportDataServer <- function(id, state, pyKinetics) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$download_curves <- downloadHandler(
      filename = function() {
        paste0("raw_curves_KinGenie_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(state$fit_dataset_loaded)

        py_df <- pyKinetics$create_export_df("raw")
        df <- pandas_to_r(py_df)

        write.csv(df, file = file, row.names = FALSE)
      }
    )

    output$download_fitted_curves <- downloadHandler(
      filename = function() {
        paste0("fitted_curves_KinGenie_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(state$fit_dataset_loaded)
        req(state$kinetics_fit_done)

        py_df <- pyKinetics$create_export_df("fit")
        df <- pandas_to_r(py_df)

        write.csv(df, file = file, row.names = FALSE)
      }
    )
  })
}
