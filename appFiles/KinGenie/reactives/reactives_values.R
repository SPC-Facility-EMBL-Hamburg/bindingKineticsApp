# Set of initial reactiveValues that are useful to manipulate the UI elements

reactives <- reactiveValues(

    logbook              = list(), # record data manipulation steps
    sample_plate_loaded  = FALSE,   # show plots/tables
    traces_loaded        = FALSE,  # show plots/tables
    fit_dataset_loaded   = FALSE,
    ss_fit_done          = FALSE,
    kinetics_fit_done    = FALSE,
    simulation_run       = FALSE,
    simulation_results   = NULL,
    relaxation_plot_done = FALSE,
    smax_update_allowed  = FALSE,
    ss_table_shown       = FALSE,
    ss_plot_shown        = FALSE,
    kinetics_table_shown = FALSE,
    kinetics_ci95_table_shown = FALSE,
    surface_based_binding     = TRUE,
    diagnostic_plots_done     = FALSE,
    update_flag               = TRUE,
    is_single_cycle           = FALSE,
    solution_model            = NULL,
    k_obs_plot_shown          = FALSE,
    kinetics_table_shown_sol  = FALSE # Show the kinetics table (solution-based)
)

# TO - DO Make the reactiveValues available to the ui

# Expose surface_based_binding to the ui code
output$surface_based_binding   <- reactive( { return( reactives$surface_based_binding  ) } )
outputOptions(output, "surface_based_binding" , suspendWhenHidden = FALSE)

