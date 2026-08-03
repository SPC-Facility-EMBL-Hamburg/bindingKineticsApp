box::use(
  shiny[
    moduleServer,
    reactiveValues
  ]
)

#' @export
simResultsServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    reactiveValues(
      available = FALSE,
      data = NULL,
      relaxation_plot_done = FALSE
    )
  })
}
