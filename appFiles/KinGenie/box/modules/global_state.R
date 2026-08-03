box::use(
  shiny[
    moduleServer,
    reactiveValues
  ]
)

#' @export
appStateServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    reactiveValues(
      traces_loaded = FALSE,
      surface_based_binding = TRUE,
      selected_experiment = NULL,
      legend_version = 0,
      ligand_info_version = 0,
      selected_trace = NULL,
      selected_color = "#E41A1C",
      kinetics_table_shown = FALSE,
      ss_plot_shown = FALSE,
      ss_table_shown = FALSE,
      kinetics_ci_95_table_shown = FALSE,
      diagnostics_plot_shown = FALSE,
      fit_dataset_loaded = FALSE,
      ss_fit_done = FALSE,
      kinetics_fit_done = FALSE,
      show_ligand_info = FALSE,
      is_single_cycle = FALSE,
      k_obs_plot_shown = FALSE
    )
  })
}
