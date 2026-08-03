box::use(
  .. / .. / plot_functions[
    plot_relaxation_rates,
    plot_simulation,
    plotlyOutput600px
  ],
  plotly[
    renderPlotly
  ],
  shiny[
    moduleServer,
    NS,
    renderUI,
    req,
    tabPanel,
    tagList,
    uiOutput
  ],
  shinydashboard[
    tabBox
  ]
)

#' @export
simPlotsTabBoxUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("sim_plots_tab_box"))
  )
}

#' @export
simPlotsTabBoxServer <- function(id, sim_state, sim_results, plot_config_reac) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$sim_plots_tab_box <- renderUI({
      tabBox(
        title = "", width = 12, id = "tabset_sim",
        tabPanel("Signal", plotlyOutput600px(ns("signal_sim_plot"))),
        if (sim_results$relaxation_plot_done) {
          tabPanel("Relaxation rates", plotlyOutput600px(ns("plot_relaxation_rates")))
        }
      )
    })

    output$plot_relaxation_rates <- renderPlotly({
      req(sim_results$available)
      req(sim_results$data)
      req(sim_results$relaxation_plot_done)

      k_obs_dominant_per_pc <- sim_results$data$k_obs_dominant_per_pc
      k_obs_non_dominant_per_pc <- sim_results$data$k_obs_non_dominant_per_pc
      prot_concs <- sim_results$data$prot_concs
      lig_concs <- sim_results$data$lig_concs

      fig <- plot_relaxation_rates(
        k_obs_dominant_per_pc,
        k_obs_non_dominant_per_pc,
        prot_concs,
        lig_concs,
        plot_config = plot_config_reac()
      )
    })

    output$signal_sim_plot <- renderPlotly({
      req(sim_results$available)
      req(sim_results$data)

      if (sim_state$model_type_sim == "surface") {
        association_time_per_cycle <- sim_results$data$association_time_per_cycle
        dissociation_time_per_cycle <- sim_results$data$dissociation_time_per_cycle
        signal_per_cycle_assoc <- sim_results$data$signal_per_cycle_assoc
        signal_per_cycle_disso <- sim_results$data$signal_per_cycle_disso
        smax_all <- sim_results$data$smax_all
        lig_concs <- sim_results$data$lig_concs

        fig <- plot_simulation(
          association_time_per_cycle,
          dissociation_time_per_cycle,
          signal_per_cycle_assoc,
          signal_per_cycle_disso,
          smax_all,
          lig_concs,
          plot_width = plot_config_reac()$plot_width,
          plot_height = plot_config_reac()$plot_height,
          plot_type = plot_config_reac()$plot_type,
          font_size = plot_config_reac()$font_size,
          show_grid_x = plot_config_reac()$show_grid_x,
          show_grid_y = plot_config_reac()$show_grid_y,
          marker_size = plot_config_reac()$marker_size,
          is_single_cycle = sim_state$is_single_cycle_sim
        )
      } else {
        interaction_time <- sim_results$data$interaction_time
        signal_per_pc <- sim_results$data$signal_per_pc
        prot_concs <- sim_results$data$prot_concs
        lig_concs <- sim_results$data$lig_concs

        fig <- plot_simulation(list(interaction_time), # Convert to list (we only have one cycle here)
          NULL, # No dissociation time
          list(signal_per_pc), # Convert to list (we only have one cycle here)
          NULL, # No dissociation signal
          prot_concs,
          lig_concs,
          plot_width  = plot_config_reac()$plot_width,
          plot_height = plot_config_reac()$plot_height,
          plot_type   = plot_config_reac()$plot_type,
          font_size   = plot_config_reac()$font_size,
          show_grid_x = plot_config_reac()$show_grid_x,
          show_grid_y = plot_config_reac()$show_grid_y,
          marker_size = plot_config_reac()$marker_size
        )
      }

      fig
    })
  })
}
