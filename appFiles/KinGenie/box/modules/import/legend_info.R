box::use(
  .. / .. / tables[
    get_plotting_df,
    get_rtable_legend
  ],
  rhandsontable[
    hot_to_r,
    renderRHandsontable,
    rhandsontable,
    rHandsontableOutput
  ],
  shiny[
    column,
    fluidRow,
    moduleServer,
    NS,
    observeEvent,
    reactiveVal,
    req,
    tagList
  ]
)

#' @export
legendUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(column(width = 12, rHandsontableOutput(ns("legendInfo"))))
  )
}

#' @export
legendServer <- function(id, state, pyKinetics) {
  moduleServer(id, function(input, output, session) {
    legend_df <- reactiveVal()

    observeEvent(state$legend_version,
      {
        req(state$traces_loaded)

        if (state$surface_based_binding) {
          labels <- unlist(pyKinetics$get_experiment_properties("sensor_names"))
          ids <- unlist(pyKinetics$get_experiment_properties("sensor_names_unique"))
        } else {
          labels <- unlist(pyKinetics$get_experiment_properties("traces_names"))
          ids <- unlist(pyKinetics$get_experiment_properties("traces_names_unique"))
        }

        legend_df(get_plotting_df(ids, labels))
      },
      ignoreInit = TRUE
    )

    output$legendInfo <- renderRHandsontable({
      req(legend_df())
      get_rtable_legend(legend_df())
    })

    observeEvent(input$legendInfo, {
      legend_df(hot_to_r(input$legendInfo))
    })

    list(
      df = legend_df
    )
  })
}
