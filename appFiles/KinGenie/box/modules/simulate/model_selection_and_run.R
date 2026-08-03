box::use(
  . / export[
    simResultsExportUI
  ],
  .. / .. / dialogs[
    pop_up_warning,
    print_model_message
  ],
  reticulate[
    np_array
  ],
  shiny[
    actionButton,
    column,
    fluidRow,
    HTML,
    icon,
    moduleServer,
    NS,
    numericInput,
    observeEvent,
    p,
    req,
    selectInput,
    span,
    tagList,
    updateSelectInput
  ],
  shinydashboard[
    box
  ],
  tippy[
    tippy_this
  ],
  utils[
    tail
  ]
)

#' @export
modelSelectionAndRunUI <- function(id) {
  ns <- NS(id)
  tagList(
    box(
      title = "1. Model", width = 3, solidHeader = T, status = "primary",
      fluidRow(
        column(6, p(
          HTML("<b>Type</b>"),
          span(icon("info-circle"), id = "info_uu_sim_2-0"),
          selectInput(ns("model_type_sim"), NULL,
            c(
              "Surface binding" = "surface",
              "In-solution" = "solution"
            ),
            selectize = FALSE
          ), # Set selectize to FALSE to avoid text overflowing outside the box
          tippy_this(
            elementId = "info_uu_sim_2-0",
            tooltip = "Choose if the target protein is bound to a surface (e.g., BLI) or in solution (e.g., NMR).
                    In the case of surface binding, the ligand concentration remains constant.
                    In the case of in-solution binding, ligand depletion takes place.",
            placement = "right"
          )
        )),
        column(6, p(
          HTML("<b>Model</b>"),
          span(icon("info-circle"), id = "info_uu_sim_2-1"),
          selectInput(ns("model_selected_sim"), NULL,
            c(
              "1:1" = "one_site",
              "1:1 (induced fit)" = "one_site_induced_fit",
              "1:1 (conformational selection)" = "one_site_conformational_selection" # ,
              # "1:1 (two ligands)"  = "heteregeneous_ligand"
            ),
            selectize = FALSE
          ), # Set selectize to FALSE to avoid text overflowing outside the box

          tippy_this(
            elementId = "info_uu_sim_2-1",
            tooltip = "Select the model for the simulation.
                    More information in the User Guide.", placement = "right"
          )
        ))
      ),
      fluidRow(
        column(7, p(
          HTML('<p style="margin-bottom:0px;"><br></p>'),
          actionButton(
            inputId = ns("btn_cal_simulation"), label = "Run simulation",
            icon("meteor"),
            style = "color: #fff; background-color: #337ab7;
                border-color: #2e6da4"
          )
        )),
        column(5, p(
          HTML("<b>Time step (s)</b>"),
          span(icon("info-circle"), id = "info_uu_sim_2-2"),
          numericInput(ns("time_step"), NULL, 0.5, min = 0, max = 10),
          tippy_this(
            elementId = "info_uu_sim_2-2",
            tooltip = "Time step to simulate the association and dissociation curves.",
            placement = "right"
          )
        )),
        simResultsExportUI("simResultsExport")
      )
    )
  )
}

#' @export
modelSelectionAndRunServer <- function(id, sim_state, sim_results, pykingenie) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$model_type_sim, {
      model_type_sim <- input$model_type_sim
      sim_state$model_type_sim <- model_type_sim

      if (model_type_sim == "surface") {
        updateSelectInput(session, "model_selected_sim", choices = c(
          "1:1" = "one_site",
          "1:1 (mass transport limitation)" = "one_site_mtl",
          "1:1 (induced fit)" = "one_site_induced_fit",
          "1:1 (conformational selection)" = "one_site_conformational_selection",
          "Heterogeneous analyte" = "heterogeneous_analyte",
          "2:1 (ligand has two binding sites)" = "ligand_has_two_sites"
        ))
      } else {
        updateSelectInput(session, "model_selected_sim", choices = c(
          "1:1" = "one_site",
          "1:1 (induced fit)" = "one_site_induced_fit",
          "1:1 (conformational selection)" = "one_site_conformational_selection"
        ))
      }
    })

    observeEvent(input$model_selected_sim, {
      sim_state$model_selected_sim <- input$model_selected_sim
    })

    observeEvent(input$time_step, {
      sim_state$time_step <- input$time_step
    })


    observeEvent(list(input$model_type_sim, input$model_selected_sim), {
      model_type <- input$model_type_sim
      model_sel <- input$model_selected_sim

      sim_results$available <- FALSE
      sim_results$data <- NULL

      if (model_type == "solution") {
        if (model_sel == "one_site_conformational_selection") {
          # Example values base on Fabian Paul ,Thomas R. Weikl, 2016:
          sim_state$total_time <- 1
          sim_state$kon_sim_1to1_adv <- 100
          sim_state$koff_sim_1to1_adv <- 1
          sim_state$kc_sim_1to1_adv <- 10
          sim_state$krev_sim_1to1_adv <- 100
          sim_state$protein_conc_sim <- 1
          sim_state$init_lig_sim_solution <- 1.8
          sim_state$lig_dil_factor_sim <- 1.6
          sim_state$time_step <- 0.02
        }

        if (model_sel == "one_site_induced_fit") {
          # Example values based on Fabian Paul ,Thomas R. Weikl, 2016:
          sim_state$total_time <- 1
          sim_state$kon_sim_1to1_adv <- 100
          sim_state$koff_sim_1to1_adv <- 100
          sim_state$kc_sim_1to1_adv <- 1
          sim_state$krev_sim_1to1_adv <- 10
          sim_state$protein_conc_sim <- 0.5
          sim_state$init_lig_sim_solution <- 1.8
          sim_state$lig_dil_factor_sim <- 1.6
          sim_state$time_step <- 0.02
        }

        if (model_sel == "one_site") {
          # Example values base on Fabian Paul ,Thomas R. Weikl, 2016:
          sim_state$total_time <- 1
          sim_state$kd_sim_1to1 <- 0.1
          sim_state$koff_sim_1to1 <- 0.01
          sim_state$protein_conc_sim <- 1
          sim_state$init_lig_sim_solution <- 5
          sim_state$lig_dil_factor_sim <- 1.6
          sim_state$time_step <- 0.02
        }
      } else {
        sim_state$time_step <- 0.5
        sim_state$protein_smax_sim <- 5
        sim_state$association_time <- 300
        sim_state$dissociation_time <- 600
        sim_state$numb_cycles_sim <- 1
        sim_state$kd_sim_1to1 <- 0.5
        sim_state$koff_sim_1to1 <- 0.01
        sim_state$kon_sim_1to1_adv <- 0.5
        sim_state$koff_sim_1to1_adv <- 0.01
        sim_state$kc_sim_1to1_adv <- 1
        sim_state$krev_sim_1to1_adv <- 10
      }
    })

    observeEvent(input$btn_cal_simulation, {
      sim_results$available <- FALSE
      sim_results$data <- NULL

      model_type_sim <- input$model_type_sim
      model_selected_sim <- input$model_selected_sim

      print_model_message(model_type_sim, model_selected_sim)

      # 0. Default number of cycles to one, it will change if we have single-cycle-kinetics
      numb_cycles <- 1

      # 1. Obtain the ligand concentration
      if (model_type_sim == "surface") {
        init_lig_sim <- sim_state$init_lig_sim_surface
      } else {
        init_lig_sim <- sim_state$init_lig_sim_solution
      }

      lig_concs <- init_lig_sim / (sim_state$lig_dil_factor_sim^(0:floor(sim_state$numb_dil_sim)))

      # Sort the lig_concs from lowest to highest
      lig_concs <- sort(lig_concs)

      # 2. Obtain the equilibrium dissociation constant and dissociation rate
      if (model_selected_sim %in% c("one_site", "one_site_mtl")) {
        Kd <- sim_state$kd_sim_1to1
        k_off <- sim_state$koff_sim_1to1
      }

      if (model_selected_sim == "one_site_mtl") {
        k_tr <- sim_state$ktr_sim_1to1
      }

      if (model_selected_sim == "ligand_has_two_sites") {
        k_on <- sim_state$kon_sim_ligand_two_sites
        k_off <- sim_state$koff_sim_ligand_two_sites
        coop_factor <- sim_state$coop_sim_ligand_two_sites
      }

      # 2. Obtain the (dis)association rates, and the conformational change rates
      if (model_selected_sim %in% c("one_site_induced_fit", "one_site_conformational_selection")) {
        k_on <- sim_state$kon_sim_1to1_adv
        k_off <- sim_state$koff_sim_1to1_adv
        k_c <- sim_state$kc_sim_1to1_adv
        k_rev <- sim_state$krev_sim_1to1_adv
      }

      if (model_selected_sim == "heterogeneous_analyte") {
        Kd1 <- sim_state$kd1_sim_hetero
        k_off1 <- sim_state$koff1_sim_hetero
        Kd2 <- sim_state$kd2_sim_hetero
        k_off2 <- sim_state$koff2_sim_hetero
        pop1 <- sim_state$pop1_sim_hetero / 100
        pop2 <- 1 - pop1

        smax1 <- sim_state$pop1_sim_smax
        smax2 <- sim_state$total_smax_sim - smax1

        if (smax1 < 0 || smax2 < 0) {
          pop_up_warning("The Smax of each population must be greater than 0.
                    Change the values of the Total Smax and the Smax of Population 1.")
          req(FALSE)
        }

        np_pop_fractions <- np_array(c(pop1, pop2))
        np_max_signal <- np_array(c(smax1, smax2))
        np_Kds <- np_array(c(Kd1, Kd2))
        np_k_offs <- np_array(c(k_off1, k_off2))

        if (pop1 + pop2 != 1) {
          pop_up_warning("The sum of the populations must be equal to 100%")
          req(FALSE)
        }

        if ((smax1 + smax2) != sim_state$total_smax_sim) {
          pop_up_warning("The sum of the populations Smax must be equal to the total Smax")
          req(FALSE)
        }
      }

      # Surface binding
      if (model_type_sim == "surface") {
        # Create a copy of lig_concs for surface-based binding
        lig_concs_init <- lig_concs

        # Models that require saving the signal of individual species
        multi_state_models <- c(
          "one_site_induced_fit",
          "one_site_conformational_selection",
          "heterogeneous_analyte",
          "ligand_has_two_sites"
        )

        # 3. Obtain the association and dissociation times
        association_time <- seq(0, sim_state$association_time, sim_state$time_step)
        dissociation_time <- seq(0, sim_state$dissociation_time, sim_state$time_step) + max(association_time)

        association_time_per_cycle <- list()
        dissociation_time_per_cycle <- list()

        association_time_per_cycle[[1]] <- association_time
        dissociation_time_per_cycle[[1]] <- dissociation_time

        dissociation_time_shift <- dissociation_time - min(dissociation_time)

        # Create the signal list. One sublist per cylce,
        # Each sublist has as many subsublists as protein concentrations (or smaxs)
        # Then, we will have, in each subsublist element, subsubsublists for each ligand concentration
        # For example, signal_per_cycle_assoc[[1]][[1]][[2]] will be the association signal
        # for the first cycle, the first protein concentration and the second ligand concentration

        signal_per_cycle_assoc <- list()
        signal_per_cycle_disso <- list()

        signal_per_pc_assoc <- list()
        signal_per_pc_disso <- list()

        if (model_selected_sim %in% multi_state_models) {
          # To store the signal of each species
          # For example, the signal produced by the species ES of the induced fit model
          signal_per_cycle_disso_adv <- list()
          signal_per_pc_disso_adv <- list()
        }

        # Obtain the smax
        if (model_selected_sim == "heterogeneous_analyte") {
          smax_all <- sim_state$total_smax_sim / (sim_state$prot_dil_factor_sim^(0:floor(sim_state$numb_dil_sim_prot)))
        } else if (model_selected_sim == "ligand_has_two_sites") {
          pl_rmax_all <- sim_state$pl_rmax_sim_two_sites / (sim_state$prot_dil_factor_sim^(0:floor(sim_state$numb_dil_sim_prot)))
          pl2_rmax_all <- sim_state$pl2_rmax_sim_two_sites / (sim_state$prot_dil_factor_sim^(0:floor(sim_state$numb_dil_sim_prot)))

          # Just to use the for-loop later
          smax_all <- 10 / (sim_state$prot_dil_factor_sim^(0:floor(sim_state$numb_dil_sim_prot)))
        } else {
          smax_all <- sim_state$protein_smax_sim / (sim_state$prot_dil_factor_sim^(0:floor(sim_state$numb_dil_sim_prot)))
        }

        np_association_time <- np_array(association_time)
        np_dissociation_time <- np_array(dissociation_time)
        np_dissociation_time_shift <- np_array(dissociation_time_shift)

        # Use only the first concentration if we have single-cycle kinetics
        if (sim_state$is_single_cycle_sim) {
          # lig_concs2 is lig_concs without the first concentration
          lig_concs2 <- lig_concs_init[-1]
          # Set lig_concs to the first concentration only
          lig_concs <- lig_concs_init[1]

          # Set the number of cycles as the number of ligand concentrations
          numb_cycles <- length(lig_concs_init)
        }

        smax_counter <- 0
        for (smax in smax_all) {
          smax_counter <- smax_counter + 1

          # Create the signal sublist. One element per ligand concentration
          signal_a <- list()
          signal_d <- list()

          if (model_selected_sim %in% multi_state_models) {
            # To store the signal of each species, for this particular smax
            signal_d_adv <- list()
          }

          for (lc in lig_concs) {
            if (model_selected_sim == "one_site") {
              signal <- pykingenie$one_site_association_analytical(np_association_time, 0, smax, k_off, Kd, lc)
            }

            if (model_selected_sim == "one_site_mtl") {
              signal <- pykingenie$solve_ode_one_site_mass_transport_association(
                np_association_time,
                0, lc / 2, lc, Kd, k_off, k_tr, smax
              )
            }

            if (model_selected_sim == "heterogeneous_analyte") {
              np_signal_start <- np_array(c(0, 0))

              signal_matrix <- pykingenie$solve_ode_mixture_analyte_association(
                np_association_time,
                np_signal_start, lc, np_pop_fractions, np_max_signal,
                np_k_offs, np_Kds
              )

              signal_matrix <- as.matrix(signal_matrix)

              # Sum the signals of the two populations
              signal <- apply(signal_matrix, 2, sum)
            }

            if (model_selected_sim == "one_site_induced_fit") {
              signal_matrix <- pykingenie$solve_induced_fit_association(
                np_association_time, lc, k_on, k_off, k_c, k_rev, 0, 0, smax
              )

              signal_df <- as.data.frame(signal_matrix)

              colnames(signal_df) <- c("signal", "sP1L", "sP2L")

              signal <- signal_df$signal
            }

            if (model_selected_sim == "one_site_conformational_selection") {
              signal_matrix <- pykingenie$solve_conformational_selection_association(
                np_association_time, lc, k_on, k_off, k_c, k_rev,
                smax = smax
              )

              signal_df <- as.data.frame(signal_matrix)

              colnames(signal_df) <- c("signal", "sP1", "sP2")

              signal <- signal_df$signal
            }

            if (model_selected_sim == "ligand_has_two_sites") {
              signal <- pykingenie$solve_two_site_cooperative_association(
                np_association_time,
                lc,
                k_on,
                k_off,
                coop_factor,
                Rmax_PL = pl_rmax_all[smax_counter],
                Rmax_LPL = pl2_rmax_all[smax_counter],
                fPL_0 = 0,
                fLPL_0 = 0
              )

              signal_df <- as.data.frame(signal)

              colnames(signal_df) <- c("signal", "PL", "PL2")

              signal <- signal_df$signal
            }

            signal_a[[length(signal_a) + 1]] <- signal

            s0 <- signal[length(signal)] # Get the last signal value, after the association phase

            if (model_selected_sim == "one_site") {
              signal <- pykingenie$one_site_dissociation_analytical(np_dissociation_time_shift, s0, k_off)
            }

            if (model_selected_sim == "one_site_mtl") {
              signal <- pykingenie$solve_ode_one_site_mass_transport_dissociation(np_dissociation_time, s0, Kd, k_off, k_tr, smax)
            }

            if (model_selected_sim == "heterogeneous_analyte") {
              # Signal of the first and second populations
              signal1 <- tail(signal_matrix[1, ], 1)
              signal2 <- tail(signal_matrix[2, ], 1)

              np_signal_start <- np_array(c(signal1, signal2))

              signal_matrix <- pykingenie$solve_ode_mixture_analyte_dissociation(
                np_dissociation_time, np_signal_start, np_k_offs
              )

              signal_matrix <- as.matrix(signal_matrix)

              signal <- apply(signal_matrix, 2, sum)

              # Store the last column of the signal matrix
              signal_d_adv[[length(signal_d_adv) + 1]] <- signal_matrix[, ncol(signal_matrix)]
            }

            if (model_selected_sim == "ligand_has_two_sites") {
              # Signal PL and LPL at the end of the association phase
              signal_pl <- tail(signal_df$PL, 1)
              signal_lpl <- tail(signal_df$PL2, 1)

              rmax_pl <- pl_rmax_all[smax_counter]
              rmax_lpl <- pl2_rmax_all[smax_counter]

              fraction_pl <- signal_pl / rmax_pl
              fraction_lpl <- signal_lpl / rmax_lpl

              signal <- pykingenie$solve_two_site_cooperative_dissociation(
                np_dissociation_time,
                k_off,
                coop_factor,
                Rmax_PL = rmax_pl,
                Rmax_LPL = rmax_lpl,
                fPL_0 = fraction_pl,
                fLPL_0 = fraction_lpl
              )

              signal_df <- as.data.frame(signal)

              colnames(signal_df) <- c("signal", "PL", "PL2")

              signal <- signal_df$signal

              signal_d_adv[[length(signal_d_adv) + 1]] <- tail(signal_df[, c("PL", "PL2")], 1)
            }

            if (model_selected_sim == "one_site_induced_fit") {
              sP2L <- tail(signal_df$sP2L, 1)

              signal_mat <- pykingenie$solve_induced_fit_dissociation(np_dissociation_time, k_off, k_c, k_rev, s0, sP2L, smax)
              signal_df <- as.data.frame(signal_mat)
              colnames(signal_df) <- c("signal", "sP1L", "sP2L")

              signal <- signal_df$signal

              signal_d_adv[[length(signal_d_adv) + 1]] <- tail(signal_df, 1)
            }

            if (model_selected_sim == "one_site_conformational_selection") {
              sP1 <- tail(signal_df$sP1, 1)

              signal_mat <- pykingenie$solve_conformational_selection_dissociation(
                np_dissociation_time, k_off, k_c, k_rev,
                smax = smax, sP1 = sP1, sP2L = s0
              )
              signal_df <- as.data.frame(signal_mat)
              colnames(signal_df) <- c("signal", "sP1", "sP2")

              signal <- signal_df$signal

              signal_d_adv[[length(signal_d_adv) + 1]] <- tail(signal_df, 1)
            }

            signal_d[[length(signal_d) + 1]] <- signal
          }

          signal_per_pc_assoc[[length(signal_per_pc_assoc) + 1]] <- signal_a
          signal_per_pc_disso[[length(signal_per_pc_disso) + 1]] <- signal_d

          if (model_selected_sim %in% multi_state_models) {
            signal_per_pc_disso_adv[[length(signal_per_pc_disso_adv) + 1]] <- signal_d_adv
          }
        }

        signal_per_cycle_assoc[[1]] <- signal_per_pc_assoc
        signal_per_cycle_disso[[1]] <- signal_per_pc_disso

        if (model_selected_sim %in% multi_state_models) {
          signal_per_cycle_disso_adv[[1]] <- signal_per_pc_disso_adv
        }

        # Repeat the association and dissociation steps if we have more than one cycles
        if (numb_cycles > 1) {
          for (cycle in 2:numb_cycles) {
            signal_per_pc_assoc <- list()
            signal_per_pc_disso <- list()

            if (model_selected_sim %in% multi_state_models) {
              signal_per_pc_disso_adv <- list()
            }

            last_disso_time <- dissociation_time_per_cycle[[cycle - 1]]

            association_time_cycle <- association_time + max(last_disso_time)
            dissociation_time_cycle <- dissociation_time + max(last_disso_time)

            association_time_per_cycle[[cycle]] <- association_time_cycle
            dissociation_time_per_cycle[[cycle]] <- dissociation_time_cycle

            association_time_cycle_shifted <- association_time_cycle - min(association_time_cycle)
            dissociation_time_cycle_shifted <- dissociation_time_cycle - min(dissociation_time_cycle)

            np_association_time <- np_array(association_time_cycle)
            np_dissociation_time <- np_array(dissociation_time_cycle)

            np_dissociation_time_shifted <- np_array(dissociation_time_cycle_shifted)
            np_association_time_shifted <- np_array(association_time_cycle_shifted)

            smax_counter <- 0

            signal_prev_cycle <- signal_per_cycle_disso[[cycle - 1]]

            if (model_selected_sim %in% multi_state_models) {
              signal_prev_cycle_adv <- signal_per_cycle_disso_adv[[cycle - 1]]
            }

            for (smax in smax_all) {
              smax_counter <- smax_counter + 1
              disso_signal_smax <- signal_prev_cycle[[smax_counter]]

              if (model_selected_sim %in% multi_state_models) {
                disso_signal_adv <- signal_prev_cycle_adv[[smax_counter]]
              }

              signal_a <- list()
              signal_d <- list()

              if (model_selected_sim %in% multi_state_models) {
                signal_d_adv <- list()
              }

              lc_counter <- 0

              lig_concs <- lig_concs2[cycle - 1]

              for (lc in lig_concs) {
                lc_counter <- lc_counter + 1

                disso_signal_lc <- disso_signal_smax[[lc_counter]]
                last_signal <- disso_signal_lc[length(disso_signal_lc)]

                if (model_selected_sim == "one_site") {
                  signal <- pykingenie$one_site_association_analytical(np_association_time_shifted, last_signal, smax, k_off, Kd, lc)
                }

                if (model_selected_sim == "one_site_mtl") {
                  signal <- pykingenie$solve_ode_one_site_mass_transport_association(
                    np_association_time, last_signal, lc / 2, lc, Kd, k_off, k_tr, smax
                  )
                }

                if (model_selected_sim == "heterogeneous_analyte") {
                  last_col <- c(disso_signal_adv[[lc_counter]])

                  np_signal_start <- np_array(last_col)

                  signal_matrix <- pykingenie$solve_ode_mixture_analyte_association(
                    np_association_time,
                    np_signal_start, lc, np_pop_fractions, np_max_signal,
                    np_k_offs, np_Kds
                  )

                  signal_matrix <- as.matrix(signal_matrix)

                  # Sum the signals of the two populations
                  signal <- apply(signal_matrix, 2, sum)
                }

                if (model_selected_sim == "one_site_induced_fit") {
                  last_row <- disso_signal_adv[[lc_counter]]

                  # Initial conditions for St - P1L - sP2L, and sP2L
                  # where St is the max signal, sP1L is the signal produced by P1L,
                  # and sP2L is the signal produced by P2L (complex after the conf. change)

                  sP1L <- last_row$sP1L[1]
                  sP2L <- last_row$sP2L[1]

                  signal_matrix <- pykingenie$solve_induced_fit_association(
                    np_association_time, lc, k_on, k_off, k_c, k_rev, sP1L, sP2L, smax
                  )

                  signal_df <- as.data.frame(signal_matrix)
                  colnames(signal_df) <- c("signal", "sP1L", "sP2L")

                  signal <- signal_df$signal
                }

                if (model_selected_sim == "one_site_conformational_selection") {
                  last_row <- disso_signal_adv[[lc_counter]]

                  # Initial conditions for sP1, and sP2L (signal)
                  # where sP1 is proportional to the amount of protein in state '1'
                  # Remember: P1 <-> P2 ; P2 + L <-> P2L

                  sP1 <- last_row$sP1[1]
                  sP2L <- last_row$signal[1]

                  signal_matrix <- pykingenie$solve_conformational_selection_association(
                    np_association_time, lc, k_on, k_off, k_c, k_rev,
                    smax = smax, sP1 = sP1, sP2L = sP2L
                  )

                  signal_df <- as.data.frame(signal_matrix)
                  colnames(signal_df) <- c("signal", "sP1", "sP2")

                  signal <- signal_df$signal
                }

                if (model_selected_sim == "ligand_has_two_sites") {
                  last_row <- disso_signal_adv[[lc_counter]]

                  signal_pl <- last_row$PL[1]
                  signal_lpl <- last_row$PL2[1]

                  rmax_pl <- pl_rmax_all[smax_counter]
                  rmax_lpl <- pl2_rmax_all[smax_counter]

                  fraction_pl <- signal_pl / rmax_pl
                  fraction_lpl <- signal_lpl / rmax_lpl

                  signal <- pykingenie$solve_two_site_cooperative_association(
                    np_association_time,
                    lc,
                    k_on,
                    k_off,
                    coop_factor,
                    Rmax_PL = rmax_pl,
                    Rmax_LPL = rmax_lpl,
                    fPL_0 = fraction_pl,
                    fLPL_0 = fraction_lpl
                  )

                  signal_df <- as.data.frame(signal)

                  colnames(signal_df) <- c("signal", "PL", "PL2")

                  signal <- signal_df$signal
                }

                signal_a[[length(signal_a) + 1]] <- signal

                s0 <- signal[length(signal)]

                if (model_selected_sim == "one_site") {
                  signal <- pykingenie$one_site_dissociation_analytical(np_dissociation_time_shifted, s0, k_off)
                }

                if (model_selected_sim == "one_site_mtl") {
                  signal <- pykingenie$solve_ode_one_site_mass_transport_dissociation(
                    np_dissociation_time, s0, Kd, k_off, k_tr, smax
                  )
                }

                if (model_selected_sim == "heterogeneous_analyte") {
                  # Signal of the first and second populations
                  signal1 <- tail(signal_matrix[1, ], 1)
                  signal2 <- tail(signal_matrix[2, ], 1)

                  np_signal_start <- np_array(c(signal1, signal2))

                  signal_matrix <- pykingenie$solve_ode_mixture_analyte_dissociation(
                    np_dissociation_time, np_signal_start, np_k_offs
                  )

                  signal_matrix <- as.matrix(signal_matrix)

                  signal <- apply(signal_matrix, 2, sum)
                }

                if (model_selected_sim == "one_site_induced_fit") {
                  sP2L <- tail(signal_df$sP2L, 1)

                  signal_mat <- pykingenie$solve_induced_fit_dissociation(np_dissociation_time, k_off, k_c, k_rev, s0, sP2L, smax)
                  signal_df <- as.data.frame(signal_mat)
                  colnames(signal_df) <- c("signal", "sP1L", "sP2L")

                  signal <- signal_df$signal
                }

                if (model_selected_sim == "one_site_conformational_selection") {
                  sP1 <- tail(signal_df$sP1, 1)

                  signal_mat <- pykingenie$solve_conformational_selection_dissociation(
                    np_dissociation_time, k_off, k_c, k_rev,
                    smax = smax, sP1 = sP1, sP2L = s0
                  )
                  signal_df <- as.data.frame(signal_mat)
                  colnames(signal_df) <- c("signal", "sP1", "sP2")

                  signal <- signal_df$signal
                }

                if (model_selected_sim == "ligand_has_two_sites") {
                  fraction_pl <- tail(signal_df$PL, 1) / rmax_pl
                  fraction_lpl <- tail(signal_df$PL2, 1) / rmax_lpl

                  signal <- pykingenie$solve_two_site_cooperative_dissociation(
                    np_dissociation_time,
                    k_off,
                    coop_factor,
                    Rmax_PL = rmax_pl,
                    Rmax_LPL = rmax_lpl,
                    fPL_0 = fraction_pl,
                    fLPL_0 = fraction_lpl
                  )

                  signal_df <- as.data.frame(signal)

                  colnames(signal_df) <- c("signal", "PL", "PL2")

                  signal <- signal_df$signal
                }

                signal_d[[length(signal_d) + 1]] <- signal

                if (model_selected_sim %in% c(
                  "one_site_induced_fit",
                  "one_site_conformational_selection",
                  "ligand_has_two_sites"
                )) {
                  signal_d_adv[[length(signal_d_adv) + 1]] <- tail(signal_df, 1)
                }

                if (model_selected_sim == "heterogeneous_analyte") {
                  signal_d_adv[[length(signal_d_adv) + 1]] <- signal_matrix[, ncol(signal_matrix)]
                }
              }

              signal_per_pc_assoc[[length(signal_per_pc_assoc) + 1]] <- signal_a
              signal_per_pc_disso[[length(signal_per_pc_disso) + 1]] <- signal_d

              if (model_selected_sim %in% multi_state_models) {
                signal_per_pc_disso_adv[[length(signal_per_pc_disso_adv) + 1]] <- signal_d_adv
              }
            }

            signal_per_cycle_assoc[[cycle]] <- signal_per_pc_assoc
            signal_per_cycle_disso[[cycle]] <- signal_per_pc_disso

            if (model_selected_sim %in% multi_state_models) {
              signal_per_cycle_disso_adv[[cycle]] <- signal_per_pc_disso_adv
            }
          }
        }

        sim_results$data <- list(
          "association_time_per_cycle"  = association_time_per_cycle,
          "dissociation_time_per_cycle" = dissociation_time_per_cycle,
          "signal_per_cycle_assoc"      = signal_per_cycle_assoc,
          "signal_per_cycle_disso"      = signal_per_cycle_disso,
          "smax_all"                    = smax_all,
          "lig_concs"                   = lig_concs_init
        )
      }

      # Surface binding
      if (model_type_sim == "solution") {
        model_is_conf_sel <- model_selected_sim == "one_site_conformational_selection"
        model_is_if <- model_selected_sim == "one_site_induced_fit"

        bound_signal_is_equal <- model_is_conf_sel || (model_is_if && sim_state$signal_ES_int == sim_state$signal_ES)

        time_step <- sim_state$time_step
        total_time <- sim_state$total_time

        k_obs_dominant_per_pc <- list()
        k_obs_non_dominant_per_pc <- list()

        if (total_time / time_step > 1e5) {
          pop_up_warning(
            "The total number of steps is larger than 10000.
                        Please reduce the interaction time or increase the time step."
          )

          req(FALSE)
        }

        # 3. Obtain the association and dissociation times
        interaction_time <- seq(0, total_time, time_step)

        signal_per_pc <- list()

        # Obtain the smax
        prot_concs <- sim_state$protein_conc_sim / (sim_state$prot_dil_factor_sim^(0:sim_state$numb_dil_sim_prot))

        for (pc in prot_concs) {
          k_obs_dominant <- c()
          k_obs_non_dominant <- c()

          signal <- list()

          for (lc in lig_concs) {
            if (model_selected_sim == "one_site") {
              signal_single <- pykingenie$signal_ode_one_site_insolution(
                np_array(interaction_time), k_off, Kd, pc, lc,
                signal_a = sim_state$signal_E,
                signal_b = sim_state$signal_S,
                signal_complex = sim_state$signal_ES_simple
              )
            }

            if (model_is_if) {
              # Initial concentration of E·S, and ES
              y <- list(0, 0)

              signal_single <- pykingenie$signal_ode_induced_fit_insolution(
                np_array(interaction_time), y, k_on, k_off, k_c, k_rev,
                E_tot = pc, S_tot = lc,
                t0 = 0, signal_E = sim_state$signal_E, signal_S = sim_state$signal_S,
                signal_ES_int = sim_state$signal_ES_int, signal_ES = sim_state$signal_ES
              )

              k_obs <- pykingenie$get_kobs_induced_fit(lc, pc, k_rev, k_c, k_on, k_off)
              k_obs_dominant <- c(k_obs_dominant, k_obs)

              k_obs2 <- pykingenie$get_kobs_induced_fit(lc, pc, k_rev, k_c, k_on, k_off, dominant = FALSE)
              k_obs_non_dominant <- c(k_obs_non_dominant, k_obs2)
            }

            if (model_is_conf_sel) {
              # Concentration of E1, E2, S, and E2S
              e_concs <- pykingenie$get_initial_concentration_conformational_selection(pc, k_c, k_rev)

              y <- list(e_concs[2], 0) # Initial concentrations of E2 and E2S

              signal_single <- pykingenie$signal_ode_conformational_selection_insolution(
                np_array(interaction_time), y, k_c, k_rev, k_on, k_off,
                E_tot = pc,
                S_tot = lc,
                t0 = 0,
                signal_E1 = sim_state$signal_E,
                signal_E2 = sim_state$signal_E,
                signal_S = sim_state$signal_S,
                signal_E2S = sim_state$signal_ES_simple
              )

              k_obs <- pykingenie$get_kobs_conformational_selection(lc, pc, k_rev, k_c, k_on, k_off)
              k_obs_dominant <- c(k_obs_dominant, k_obs)

              k_obs2 <- pykingenie$get_kobs_conformational_selection(lc, pc, k_rev, k_c, k_on, k_off, dominant = FALSE)
              k_obs_non_dominant <- c(k_obs_non_dominant, k_obs2)
            }

            signal[[length(signal) + 1]] <- signal_single
          }

          signal_per_pc[[length(signal_per_pc) + 1]] <- signal

          if (model_is_if || model_is_conf_sel) {
            k_obs_dominant_per_pc[[length(k_obs_dominant_per_pc) + 1]] <- k_obs_dominant
            k_obs_non_dominant_per_pc[[length(k_obs_non_dominant_per_pc) + 1]] <- k_obs_non_dominant
          }
        }

        sim_results$relaxation_plot_done <- (model_is_if || model_is_conf_sel) && bound_signal_is_equal

        sim_results$data <- list(
          "interaction_time"  = interaction_time,
          "signal_per_pc"     = signal_per_pc,
          "prot_concs"        = prot_concs,
          "lig_concs"         = lig_concs
        )
      }

      sim_results$is_single_cycle <- sim_state$is_single_cycle_sim
      sim_results$model_type <- model_type_sim
      sim_results$available <- TRUE

    })
  })
}
