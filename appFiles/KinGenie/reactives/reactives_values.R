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
    kinetics_table_shown_sol  = FALSE, # Show the kinetics table (solution-based)

    selected_color = "#E41A1C",
    selected_trace = NULL,

    plot_config = list(

        x_label = "Time (s)", # x-axis label
        y_label = "Signal (a.u.)", # y-axis label
        width = 1000, # plot width in pixels
        height = 600, # plot height in pixels
        type = 'png', # plot type
        axis_size = 16, # axis size
        x_legend_pos = 1, # x legend position
        y_legend_pos = 1, # y legend position
        color_bar_length = 0.5, # color bar length
        show_colorbar = TRUE, # show colorbar
        show_grid_x = FALSE, # show x grid
        show_grid_y = FALSE, # show y grid
        marker_size = 5, # marker size
        line_width = 2, # line width
        max_points = 4000, # max points to display per subplot
        n_xticks = 4, # number of ticks on x-axis
        n_yticks = 3, # number of ticks on y-axis
        tick_length = 8, # tick length
        tick_width = 2, # tick width
        split_by_smax = TRUE,
        smooth_curves = FALSE,
        rolling_window_size = 0.5

    )

)

# TO - DO Make the reactiveValues available to the ui

# Expose surface_based_binding to the ui code
output$surface_based_binding   <- reactive( { return( reactives$surface_based_binding  ) } )
outputOptions(output, "surface_based_binding" , suspendWhenHidden = FALSE)

