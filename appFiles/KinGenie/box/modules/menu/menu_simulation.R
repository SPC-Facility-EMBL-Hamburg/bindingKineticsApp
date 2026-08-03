box::use(
  .. / simulate / binding_params[
    bindingParamsUI
  ],
  .. / simulate / experimental_params[
    experimentalParamsUI
  ],
  .. / simulate / model_selection_and_run[
    modelSelectionAndRunUI
  ],
  .. / simulate / plot_result[
    simPlotsTabBoxUI
  ],
  shiny[
    fluidRow,
    NS,
    tagList
  ]
)

#' @export
menuSimulationUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      modelSelectionAndRunUI("modelSelectionAndRun"),
      experimentalParamsUI("experimentalParams"),
      bindingParamsUI("bindingParams"),
      simPlotsTabBoxUI("simPlotTabBox")
    )
  )
}
