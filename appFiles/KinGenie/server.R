options(shiny.maxRequestSize = 100 * 1024^2)
options(stringsAsFactors = FALSE)

box::use(
  . / box / dialogs[
    welcome_message
  ],
  . / box / modules / analysis / fit_controls[
    fitControlsServer
  ],
  . / box / modules / analysis / fit_data[
    fitDataServer
  ],
  . / box / modules / analysis / fit_results[
    fitResultsServer
  ],
  . / box / modules / export / data[
    exportDataServer
  ],
  . / box / modules / export / logbook[
    exportLogbookServer
  ],
  . / box / modules / global_logbook[
    logbookServer
  ],
  . / box / modules / global_plot_config[
    plotConfigServer
  ],
  . / box / modules / global_state[
    appStateServer
  ],
  . / box / modules / import / legend_info[
    legendServer
  ],
  . / box / modules / import / load_data[
    loadDataServer
  ],
  . / box / modules / import / plot_input[
    plotInputServer
  ],
  . / box / modules / import / processing[
    processingServer
  ],
  . / box / modules / plotting / appearance[
    appearanceConfigServer
  ],
  . / box / modules / plotting / axis[
    axisConfigServer
  ],
  . / box / modules / plotting / export[
    exportConfigServer
  ],
  . / box / modules / plotting / visualization[
    visualizationConfigServer
  ],
  . / box / modules / simulate / binding_params[
    bindingParamsServer
  ],
  . / box / modules / simulate / experimental_params[
    experimentalParamsServer
  ],
  . / box / modules / simulate / export[
    simResultsExportServer
  ],
  . / box / modules / simulate / model_selection_and_run[
    modelSelectionAndRunServer
  ],
  . / box / modules / simulate / plot_result[
    simPlotsTabBoxServer
  ],
  . / box / modules / simulate / results_data[
    simResultsServer
  ],
  . / box / modules / simulate / state[
    simStateServer
  ],
  reticulate[
    import
  ]
)

# Import the Python package
pykingenie <- import("pykingenie")

server <- function(input, output, session) {
  user <- Sys.info()[["user"]]

  if (user != "os") {
    welcome_message(appName) # helpers.R
  }

  # To handle the general processing, the unfolding models,
  # the secondary structure calculation, and the custom models
  pyKinetics <- pykingenie$KineticsAnalyzer()

  state <- appStateServer("state")
  plot_config <- plotConfigServer("plot_config")
  logbook <- logbookServer("logbook")

  loadDataServer("load", state, pyKinetics, pykingenie, logbook)

  legend <- legendServer("legend", state, pyKinetics)
  dataset <- fitDataServer("fitData", state, pyKinetics)

  processingServer("processing", state, pyKinetics, legend$df, logbook)
  plotInputServer("plotInput", state, pyKinetics, legend$df, plot_config)

  ### --- Plot configuration for import data panel --- ###
  exportConfigServer("exportConfig", plot_config)
  axisConfigServer("axisConfig", plot_config)
  appearanceConfigServer(
    "appearanceConfig", state, plot_config,
    show_colour_controls = TRUE, legend_df_reac = legend$df
  )
  ### --- End of plot configuration for import data panel --- ###


  ### --- Plot configuration for analysis panel --- ###
  exportConfigServer("exportConfigFit", plot_config)
  axisConfigServer("axisConfigFit", plot_config)
  appearanceConfigServer(
    "appearanceConfigFit", state, plot_config,
    show_colour_controls = FALSE
  )
  visualizationConfigServer("visualizationConfigFit", state, plot_config)
  ### --- End of plot configuration for analysis panel --- ###

  fitControlsServer("fitControls", state, dataset$df, pyKinetics, logbook)
  fitResultsServer("fitResults", state, pyKinetics, plot_config, logbook)

  ### --- Export data modules --- ###
  exportDataServer("exportData", state, pyKinetics)
  exportLogbookServer("exportLogbook", logbook)
  ### --- End of export data modules --- ###

  #### --- Simulation modules --- ###
  sim_state <- simStateServer("sim_state")
  sim_results <- simResultsServer("sim_results")

  modelSelectionAndRunServer("modelSelectionAndRun", sim_state, sim_results, pykingenie)
  experimentalParamsServer("experimentalParams", sim_state)
  bindingParamsServer("bindingParams", sim_state)
  simPlotsTabBoxServer("simPlotTabBox", sim_state, sim_results, plot_config)
  simResultsExportServer("simResultsExport", sim_results)
  ###  --- End of simulation modules --- ###
}
