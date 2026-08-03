box::use(
  .. / plotting / appearance[
    appearanceConfigUI
  ],
  .. / plotting / axis[
    axisConfigUI
  ],
  .. / plotting / export[
    exportConfigUI
  ],
  .. / plotting / visualization[
    visualizationConfigUI
  ],
  shiny[
    column,
    NS,
    tagList
  ],
  shinydashboard[
    box
  ]
)

#' @export
plottingBoxFitUI <- function(id, state) {
  ns <- NS(id)

  tagList(
    box(
      title = "3. Plotting", width = 12, solidHeader = T, status = "primary",
      column(
        12,
        appearanceConfigUI("appearanceConfigFit", width = 3),
        axisConfigUI("axisConfigFit"),
        exportConfigUI("exportConfigFit"),
        visualizationConfigUI("visualizationConfigFit")
      )
    )
  )
}
