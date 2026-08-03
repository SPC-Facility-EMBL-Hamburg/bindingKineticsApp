box::use(
  .. / .. / helpers[
    pandas_to_r
  ],
  .. / .. / plot_functions[
    diagnostic_plot,
    plot_association_dissociation,
    plot_association_dissociation_residuals,
    plot_interactions,
    plot_interactions_residuals,
    plot_many_relaxation_rates,
    plot_steady_state,
    plot_steady_state_residuals,
    plotlyOutput600px
  ],
  DT[
    datatable,
    DTOutput,
    formatStyle,
    renderDT,
    styleEqual
  ],
  plotly[
    renderPlotly
  ],
  shiny[
    column,
    fluidRow,
    moduleServer,
    NS,
    renderTable,
    renderUI,
    req,
    tableOutput,
    tabPanel,
    tagList,
    uiOutput
  ],
  shinydashboard[
    tabBox
  ]
)

#' @export
fitResultsUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(
        12,
        uiOutput(ns("fitting_tabpanel"))
      )
    )
  )
}

#' @export
fitResultsServer <- function(id, state, pyKinetics, plot_config_reac, logbook) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$fitting_tabpanel <- renderUI({
      req(state$traces_loaded)

      if (state$surface_based_binding) {
        tabBox(
          title = "", width = 12, id = "tabBoxFit",
          tabPanel("Assoc. and diss. traces", plotlyOutput600px(ns("tracesAssDissFit"))),
          tabPanel("Residuals", plotlyOutput600px(ns("residuals"))),
          if (state$kinetics_table_shown) {
            tabPanel("Fitted params (Kinetics)", value = "kineticsParamsTable", tableOutput(ns("fittingInfoKinetics")))
          },
          if (state$ss_table_shown) {
            tabPanel("Fitted params (Steady-state)", tableOutput(ns("fittingInfoSS")))
          },
          tabPanel("Fitted params (boundaries)", DTOutput(ns("fittedParamsBoundaries"))),
          if (state$diagnostics_plot_shown) {
            tabPanel("Observed constants", plotlyOutput600px(ns("diagnostic_plot")))
          },
          if (state$ss_plot_shown) {
            tabPanel("Steady-state", plotlyOutput600px(ns("steady_state")))
          },
          if (state$kinetics_ci_95_table_shown) {
            tabPanel("Asymmetric error", tableOutput(ns("fittingInfoKineticsCI95")))
          }
        )
      } else {
        tabBox(
          title = "", width = 12, id = "tabBoxFitInSolution",
          tabPanel("Traces", plotlyOutput600px(ns("tracesFitSolution"))),
          tabPanel("Residuals", plotlyOutput600px(ns("residualsSolution"))),
          if (state$kinetics_fit_done) {
            tabPanel("Fitted params (Kinetics)", tableOutput(ns("fittingInfoKinetics")))
          },
          if (state$k_obs_plot_shown) {
            tabPanel("Observed constants", plotlyOutput600px(ns("kobs_plot")))
          }
        )
      }
    })

    output$tracesAssDissFit <- renderPlotly({
      req(state$fit_dataset_loaded)
      req(state$surface_based_binding)
      req(plot_config_reac())

      fig <- plot_association_dissociation(
        pyKinetics$fittings,
        plot_config = plot_config_reac(),
        plot_assoc  = TRUE,
        plot_disso  = TRUE,
        plot_fit    = state$kinetics_fit_done
      )

      fig
    })

    output$fittingInfoSS <- renderTable(
      {
        req(state$ss_fit_done)
        req(state$fit_dataset_loaded)

        dfs <- list()

        for (name in pyKinetics$fittings_names) {
          py_df <- pyKinetics$fittings[[name]]$fit_params_ss
          df <- pandas_to_r(py_df)
          dfs[[length(dfs) + 1]] <- df
        }

        do.call(rbind, dfs)
      },
      digits = 5
    )

    output$fittingInfoKinetics <- renderTable(
      {
        req(state$kinetics_fit_done)
        req(state$fit_dataset_loaded)

        dfs <- list()

        for (name in pyKinetics$fittings_names) {
          py_df <- pyKinetics$fittings[[name]]$fit_params_kinetics

          df <- pandas_to_r(py_df)
          df$Name <- name

          dfs[[length(dfs) + 1]] <- df
        }

        do.call(rbind, dfs)
      },
      digits = 5
    )


    output$fittedParamsBoundaries <- renderDT({
      req(state$fit_dataset_loaded)
      # Ask either for state$kinetics_fit_done or state$ss_fit_done
      req(state$kinetics_fit_done || state$ss_fit_done)

      dfs <- list()

      for (name in pyKinetics$fittings_names) {
        py_df <- pyKinetics$fittings[[name]]$fitted_params_boundaries
        df <- pandas_to_r(py_df)
        df$Name <- name

        dfs[[length(dfs) + 1]] <- df
      }

      # Combine all data frames into one
      # The first column is the fitted parameter, the second is the lower boundary, and the third is the upper boundary

      df <- do.call(rbind, dfs)

      # Round to the first three significative digits (of the first three columns)
      df[, 1:3] <- signif(df[, 1:3], 3)

      # Find if the difference between the fitted and the boundaries is too small
      close_values1 <- abs(df[, 2] - df[, 1]) / df[, 1] < 0.05
      close_values2 <- abs(df[, 3] - df[, 1]) / df[, 3] < 0.05

      # Check if the fitted value is too close to the lower or upper boundary
      df$highlight <- !(close_values1 | close_values2)

      df <- datatable(
        df,
        options = list(
          columnDefs = list(
            list(visible = FALSE, targets = 5) # hide the helper column if needed
          )
        )
      ) |>
        formatStyle(
          "Fitted_parameter_value",
          backgroundColor = styleEqual(c(TRUE, FALSE), c("lightgreen", "red")),
          valueColumns = "highlight"
        )

      df
    })

    output$fittingInfoKineticsCI95 <- renderTable(
      {
        req(state$kinetics_fit_done)
        req(state$fit_dataset_loaded)

        dfs <- list()

        for (name in pyKinetics$fittings_names) {
          py_df <- pyKinetics$fittings[[name]]$fit_params_kinetics_ci95
          df <- pandas_to_r(py_df)
          df$Name <- name

          dfs[[length(dfs) + 1]] <- df
        }

        do.call(rbind, dfs)
      },
      digits = 5
    )

    output$diagnostic_plot <- renderPlotly({
      req(state$traces_loaded)
      req(state$fit_dataset_loaded)
      req(state$surface_based_binding)
      req(state$diagnostics_plot_shown)

      experiments <- pyKinetics$fittings

      k_obs <- lapply(experiments, function(x) unlist(x$k_obs))
      ligand_conc <- lapply(experiments, function(x) unlist(x$lig_conc_lst))

      dfs <- list()

      for (i in seq(1, length(experiments))) {
        df <- data.frame(
          "ligand_conc"     = ligand_conc[[i]],
          "k_obs"           = k_obs[[i]],
          "experiment_name" = names(experiments)[i]
        )

        dfs[[i]] <- df
      }

      df <- do.call(rbind, dfs)

      df$experiment_name <- factor(df$experiment_name, levels = names(experiments))

      fig <- diagnostic_plot(
        df,
        plot_config = plot_config_reac()
      )

      fig
    })


    output$steady_state <- renderPlotly({
      req(state$fit_dataset_loaded)
      req(state$surface_based_binding)

      fig <- plot_steady_state(
        pyKinetics$fittings,
        plot_config = plot_config_reac(),
        plot_fit    = state$ss_fit_done
      )

      fig
    })

    output$residuals <- renderPlotly({
      req(state$fit_dataset_loaded)
      req(state$surface_based_binding)

      fig <- NULL

      if (state$kinetics_fit_done) {
        fig <- plot_association_dissociation_residuals(
          pyKinetics$fittings,
          plot_config = plot_config_reac(),
          plot_assoc  = TRUE,
          plot_disso  = TRUE,
          plot_fit    = state$kinetics_fit_done
        )
      }

      if (state$ss_fit_done) {
        fig <- plot_steady_state_residuals(
          pyKinetics$fittings,
          plot_config = plot_config_reac(),
          plot_fit    = state$ss_fit_done
        )
      }

      fig
    })

    output$tracesFitSolution <- renderPlotly({
      req(!state$surface_based_binding)
      req(state$fit_dataset_loaded)

      fig <- plot_interactions(
        pyKinetics$fittings,
        plot_config = plot_config_reac(),
        plot_fit    = state$kinetics_fit_done
      )

      fig
    })


    output$residualsSolution <- renderPlotly({
      req(!state$surface_based_binding)
      req(state$fit_dataset_loaded)
      req(state$kinetics_fit_done)

      fig <- plot_interactions_residuals(
        pyKinetics$fittings,
        plot_config = plot_config_reac()
      )

      fig
    })

    output$kobs_plot <- renderPlotly({
      req(state$traces_loaded)
      req(plot_config_reac())
      req(state$kinetics_fit_done)
      req(!state$surface_based_binding)

      # Require that we fitted single or double exponentials
      req(state$solution_model == "double")

      model <- state$solution_model

      fittings <- pyKinetics$fittings

      k_obs_dominant_per_pc_lst <- lapply(fittings, function(f) (f$k_obs_1_per_prot))
      k_obs_non_dominant_per_pc_lst <- lapply(fittings, function(f) (f$k_obs_2_per_prot))
      protein_conc_lst <- lapply(fittings, function(f) (f$unq_prot_conc))
      ligand_conc_per_pc_lst <- lapply(fittings, function(f) (f$lig_conc_per_protein))

      print(k_obs_dominant_per_pc_lst)
      print(k_obs_non_dominant_per_pc_lst)
      print(protein_conc_lst)
      print(ligand_conc_per_pc_lst)

      fig <- plot_many_relaxation_rates(
        k_obs_dominant_per_pc_lst,
        k_obs_non_dominant_per_pc_lst,
        protein_conc_lst,
        ligand_conc_per_pc_lst,
        plot_config = plot_config_reac()
      )

      fig
    })
  })
}
