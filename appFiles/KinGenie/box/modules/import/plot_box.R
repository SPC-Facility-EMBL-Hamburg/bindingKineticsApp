box::use(
  . / legend_info[
    legendUI
  ],
  .. / plotting / appearance[
    appearanceConfigUI
  ],
  .. / plotting / axis[
    axisConfigUI
  ],
  .. / plotting / export[
    exportConfigUI
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
plottingBoxUI <- function(id, state) {
  ns <- NS(id)

  tagList(
    box(
      title = "3. Plotting", width = 12, solidHeader = T, status = "primary",
      column(
        6,
        appearanceConfigUI("appearanceConfig"),
        axisConfigUI("axisConfig"),
        exportConfigUI("exportConfig")
      ),
      column(
        6,
        legendUI("legend")
      )
    )
  )
}
