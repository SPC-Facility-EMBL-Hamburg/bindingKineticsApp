observeEvent(input$is_single_cycle_sim, {

    output$signal_sim_plot <- NULL

})

observeEvent(input$model_type_sim, {

    model_type_sim <- input$model_type_sim

    if (model_type_sim == "surface") {

        updateSelectInput(session, "model_selected_sim", choices = c(

            "1:1"                              = "one_site",
            "1:1 (mass transport limitation)"  = "one_site_mtl",
            "1:1 (induced fit)"                = "one_site_induced_fit",
            "1:1 (conformational selection)"   = "one_site_conformational_selection",
            "Heterogeneous analyte"            = "heterogeneous_analyte"

        ))

    } else {

        updateSelectInput(session, "model_selected_sim", choices = c(
            "1:1"  = "one_site",
            "1:1 (induced fit)"  = "one_site_induced_fit",
            "1:1 (conformational selection)"  = "one_site_conformational_selection"
        ))

    }
})


output$pop2_sim_hetero <- renderText({
    paste0(100 - input$pop1_sim_hetero)
})

output$pop2_sim_smax <- renderText({
    paste0(input$total_smax_sim - input$pop1_sim_smax)
})

# Update smax of each population to be between 0 and the global smax
observeEvent(input$total_smax_sim,{

    reactives$smax_update_allowed <- FALSE

    smax_global <- input$total_smax_sim
    smax_1      <- input$pop1_sim_smax
    smax_2      <- smax_global - smax_1
    diff        <- smax_global - smax_1 - smax_2

    updateNumericInput(session,'pop1_sim_smax',value=smax_1+diff,max=smax_global)

})


observeEvent(input$btn_cal_simulation, {

    reactives$simulation_run            <- FALSE
    reactives$simulation_results        <- NULL
    output$signal_sim_plot              <- NULL

    # Check if the relaxation data tab exists and remove it if it does
    if (reactives$relaxation_plot_done) {
        removeTab(inputId = "tabset_sim", target = "Relaxation rates")
        reactives$relaxation_plot_done <- FALSE
    }

    model_type_sim     <- input$model_type_sim
    model_selected_sim <- input$model_selected_sim

    print_model_message(model_type_sim,model_selected_sim)

    # 0. Default number of cycles to one, it will change if we have single-cycle-kinetics
    numb_cycles <- 1

    # 1. Obtain the ligand concentration
    lig_concs       <- input$init_lig_sim / (input$lig_dil_factor_sim ^ (0:floor(input$numb_dil_sim)))

    # Sort the lig_concs from lowest to highest
    lig_concs       <- sort(lig_concs)

    # 2. Obtain the equilibrium dissociation constant and dissociation rate
    if (model_selected_sim %in% c('one_site','one_site_mtl')) {

        Kd    <- input$kd_sim_1to1
        k_off <- input$koff_sim_1to1

    }

    if (model_selected_sim == 'one_site_mtl') {

        k_tr <- input$ktr_sim_1to1

    }

    # 2. Obtain the (dis)association rates, and the conformational change rates
    if (model_selected_sim %in% c('one_site_induced_fit','one_site_conformational_selection')) {

        k_on   <- input$kon_sim_1to1_adv
        k_off  <- input$koff_sim_1to1_adv
        k_c    <- input$kc_sim_1to1_adv
        k_rev  <- input$krev_sim_1to1_adv

    }

    if (model_selected_sim == 'heterogeneous_analyte') {

        Kd1    <- input$kd1_sim_hetero
        k_off1 <- input$koff1_sim_hetero
        Kd2    <- input$kd2_sim_hetero
        k_off2 <- input$koff2_sim_hetero
        pop1   <- input$pop1_sim_hetero / 100
        pop2   <- 1 - pop1

        smax1  <- input$pop1_sim_smax
        smax2  <- input$total_smax_sim - smax1

        if (smax1 < 0 || smax2 < 0) {
            popUpWarning('The Smax of each population must be greater than 0.
            Change the values of the Total Smax and the Smax of Population 1.')
            return(NULL)
        }

        np_pop_fractions <- np_array(c(pop1,pop2))
        np_max_signal    <- np_array(c(smax1,smax2))
        np_Kds           <- np_array(c(Kd1,Kd2))
        np_k_offs        <- np_array(c(k_off1,k_off2))

        if (pop1 + pop2 != 1) {
            popUpWarning('The sum of the populations must be equal to 100%')
            return(NULL)
        }

        if ((smax1 + smax2) != input$total_smax_sim) {
            popUpWarning('The sum of the populations Smax must be equal to the total Smax')
            return(NULL)
        }

    }

    # Surface binding
    if (model_type_sim == "surface") {

        # Create a copy of lig_concs for surface-based binding
        lig_concs_init  <- lig_concs

        # Models that require saving the signal of individual species
        multi_state_models <- c(
            'one_site_induced_fit',
            'one_site_conformational_selection',
            'heterogeneous_analyte')

        # 3. Obtain the association and dissociation times
        association_time  <- seq(0,input$association_time,input$time_step)
        dissociation_time <- seq(0,input$dissociation_time,input$time_step) + max(association_time)

        association_time_per_cycle  <- list()
        dissociation_time_per_cycle <- list()

        association_time_per_cycle[[1]]  <- association_time
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
            signal_per_pc_disso_adv    <- list()

        }

        # Obtain the smax
        if (model_selected_sim == 'heterogeneous_analyte') {

            smax1_all <- smax1 / (input$prot_dil_factor_sim ^ (0:floor(input$numb_dil_sim_prot)))
            smax2_all <- smax1 / (input$prot_dil_factor_sim ^ (0:floor(input$numb_dil_sim_prot)))

            smax_all <- input$total_smax_sim / (input$prot_dil_factor_sim ^ (0:floor(input$numb_dil_sim_prot)))

        } else {

            smax_all <- input$protein_smax_sim / (input$prot_dil_factor_sim ^ (0:floor(input$numb_dil_sim_prot)))

        }

        np_association_time  <- np_array(association_time)
        np_dissociation_time <- np_array(dissociation_time)
        np_dissociation_time_shift <- np_array(dissociation_time_shift)
        # Use only the first concentration if we have single-cycle kinetics
        if (input$is_single_cycle_sim) {

            # lig_concs2 is lig_concs without the first concentration
            lig_concs2 <- lig_concs_init[-1]
            # Set lig_concs to the first concentration only
            lig_concs  <- lig_concs_init[1]

            # Set the number of cycles as the number of ligand concentrations
            numb_cycles <- length(lig_concs_init)

        }

        for (smax in smax_all) {

            # Create the signal sublist. One element per ligand concentration
            signal_a <- list()
            signal_d <- list()

            if (model_selected_sim %in% multi_state_models) {

                # To store the signal of each species, for this particular smax
                signal_d_adv <- list()

            }

            for (lc in lig_concs) {

                if (model_selected_sim == 'one_site') {

                    signal <- pykingenie$one_site_association_analytical(np_association_time,0,smax,k_off,Kd,lc)

                }

                if (model_selected_sim == 'one_site_mtl') {

                    signal <- pykingenie$solve_ode_one_site_mass_transport_association(np_association_time,
                    0,lc/2,lc,Kd,k_off,k_tr,smax)

                }

                if (model_selected_sim == 'heterogeneous_analyte') {

                    np_signal_start  <- np_array(c(0,0))

                    signal_matrix <- pykingenie$solve_ode_mixture_analyte_association(np_association_time,
                    np_signal_start,lc,np_pop_fractions,np_max_signal,
                    np_k_offs,np_Kds)

                    signal_matrix <- as.matrix(signal_matrix)

                    # Sum the signals of the two populations
                    signal <- apply(signal_matrix,2,sum)

                }

                if (model_selected_sim == 'one_site_induced_fit') {

                    signal_matrix <- pykingenie$solve_induced_fit_association(
                        np_association_time, lc, k_on, k_off, k_c, k_rev,0,0,smax
                    )

                    signal_df <- as.data.frame(signal_matrix)

                    colnames(signal_df) <- c('signal','sP1L','sP2L')

                    signal    <- signal_df$signal

                }

                if (model_selected_sim == 'one_site_conformational_selection') {

                    signal_matrix <- pykingenie$solve_conformational_selection_association(
                        np_association_time, lc, k_on, k_off, k_c, k_rev,smax=smax
                    )

                    signal_df <- as.data.frame(signal_matrix)

                    colnames(signal_df) <- c('signal','sP1','sP2')

                    signal    <- signal_df$signal

                }


                signal_a[[length(signal_a)+1]] <- signal

                s0 <- signal[length(signal)] # Get the last signal value, after the association phase

                if (model_selected_sim == 'one_site') {

                    signal <- pykingenie$one_site_dissociation_analytical(np_dissociation_time_shift,s0,k_off)

                }

                if (model_selected_sim == 'one_site_mtl') {

                    signal <- pykingenie$solve_ode_one_site_mass_transport_dissociation(np_dissociation_time,s0,Kd,k_off,k_tr,smax)

                }

                if (model_selected_sim == 'heterogeneous_analyte') {

                    # Signal of the first and second populations
                    signal1 <- tail(signal_matrix[1,],1)
                    signal2 <- tail(signal_matrix[2,],1)

                    np_signal_start  <- np_array(c(signal1,signal2))

                    signal_matrix <- pykingenie$solve_ode_mixture_analyte_dissociation(
                        np_dissociation_time,np_signal_start,np_k_offs)

                    signal_matrix <- as.matrix(signal_matrix)

                    signal <- apply(signal_matrix,2,sum)

                    # Store the last column of the signal matrix
                    signal_d_adv[[length(signal_d_adv)+1]] <- signal_matrix[,ncol(signal_matrix)]

                }

                if (model_selected_sim == 'one_site_induced_fit') {

                    sP2L <- tail(signal_df$sP2L,1)

                    signal_mat <- pykingenie$solve_induced_fit_dissociation(np_dissociation_time, k_off, k_c, k_rev,s0,sP2L,smax)
                    signal_df <- as.data.frame(signal_mat)
                    colnames(signal_df) <- c('signal','sP1L','sP2L')

                    signal    <- signal_df$signal

                    signal_d_adv[[length(signal_d_adv)+1]] <- tail(signal_df,1)
                }

                if (model_selected_sim == 'one_site_conformational_selection') {

                    sP1 <- tail(signal_df$sP1,1)

                    signal_mat <- pykingenie$solve_conformational_selection_dissociation(
                        np_dissociation_time, k_off, k_c, k_rev,smax=smax,sP1=sP1,sP2L=s0
                    )
                    signal_df <- as.data.frame(signal_mat)
                    colnames(signal_df) <- c('signal','sP1','sP2')

                    signal    <- signal_df$signal

                    signal_d_adv[[length(signal_d_adv)+1]] <- tail(signal_df,1)
                }


                signal_d[[length(signal_d)+1]] <- signal

            }

            signal_per_pc_assoc[[length(signal_per_pc_assoc)+1]] <- signal_a
            signal_per_pc_disso[[length(signal_per_pc_disso)+1]] <- signal_d

            if (model_selected_sim %in% multi_state_models) {

                signal_per_pc_disso_adv[[length(signal_per_pc_disso_adv)+1]] <- signal_d_adv

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

                last_disso_time <- dissociation_time_per_cycle[[cycle-1]]

                association_time_cycle  <- association_time  + max(last_disso_time)
                dissociation_time_cycle <- dissociation_time + max(last_disso_time)

                association_time_per_cycle[[cycle]]  <- association_time_cycle
                dissociation_time_per_cycle[[cycle]] <- dissociation_time_cycle

                association_time_cycle_shifted  <- association_time_cycle  - min(association_time_cycle)
                dissociation_time_cycle_shifted <- dissociation_time_cycle - min(dissociation_time_cycle)

                np_association_time  <- np_array(association_time_cycle)
                np_dissociation_time <- np_array(dissociation_time_cycle)

                np_dissociation_time_shifted <- np_array(dissociation_time_cycle_shifted)
                np_association_time_shifted  <- np_array(association_time_cycle_shifted)

                smax_counter <- 0

                signal_prev_cycle <- signal_per_cycle_disso[[cycle-1]]

                if (model_selected_sim %in% multi_state_models) {

                    signal_prev_cycle_adv <- signal_per_cycle_disso_adv[[cycle-1]]

                }

                for (smax in smax_all) {

                    smax_counter        <- smax_counter + 1
                    disso_signal_smax   <- signal_prev_cycle[[smax_counter]]

                    if (model_selected_sim %in% multi_state_models) {

                        disso_signal_adv <- signal_prev_cycle_adv[[smax_counter]]

                    }

                    signal_a <- list()
                    signal_d <- list()

                    if (model_selected_sim %in% multi_state_models) {

                        signal_d_adv <- list()

                    }

                    lc_counter <- 0

                    lig_concs  <- lig_concs2[cycle-1]

                    for (lc in lig_concs) {

                        lc_counter <- lc_counter + 1

                        disso_signal_lc <- disso_signal_smax[[lc_counter]]
                        last_signal     <- disso_signal_lc[length(disso_signal_lc)]

                        if (model_selected_sim == 'one_site') {

                            signal <- pykingenie$one_site_association_analytical(np_association_time_shifted,last_signal,smax,k_off,Kd,lc)

                        }

                        if (model_selected_sim == 'one_site_mtl') {

                            signal <- pykingenie$solve_ode_one_site_mass_transport_association(
                            np_association_time,last_signal,lc/2,lc,Kd,k_off,k_tr,smax)
                        }

                        if (model_selected_sim == 'heterogeneous_analyte') {

                            last_col            <- c(disso_signal_adv[[lc_counter]])

                            np_signal_start <- np_array(last_col)

                            signal_matrix <- pykingenie$solve_ode_mixture_analyte_association(np_association_time,
                            np_signal_start,lc,np_pop_fractions,np_max_signal,
                            np_k_offs,np_Kds)

                            signal_matrix <- as.matrix(signal_matrix)

                            # Sum the signals of the two populations
                            signal <- apply(signal_matrix,2,sum)

                        }

                        if (model_selected_sim == 'one_site_induced_fit') {

                            last_row            <- disso_signal_adv[[lc_counter]]

                            # Initial conditions for St - P1L - sP2L, and sP2L
                            # where St is the max signal, sP1L is the signal produced by P1L,
                            # and sP2L is the signal produced by P2L (complex after the conf. change)

                            sP1L <- last_row$sP1L[1]
                            sP2L <- last_row$sP2L[1]

                            signal_matrix <- pykingenie$solve_induced_fit_association(
                                np_association_time, lc, k_on, k_off, k_c, k_rev,sP1L,sP2L,smax
                            )

                            signal_df <- as.data.frame(signal_matrix)
                            colnames(signal_df) <- c('signal','sP1L','sP2L')

                            signal    <- signal_df$signal

                        }

                        if (model_selected_sim == 'one_site_conformational_selection') {

                            last_row            <- disso_signal_adv[[lc_counter]]

                            # Initial conditions for sP1, and sP2L (signal)
                            # where sP1 is proportional to the amount of protein in state '1'
                            # Remember: P1 <-> P2 ; P2 + L <-> P2L

                            sP1  <- last_row$sP1[1]
                            sP2L <- last_row$signal[1]

                            signal_matrix <- pykingenie$solve_conformational_selection_association(
                                np_association_time, lc, k_on, k_off, k_c, k_rev,smax=smax,sP1=sP1,sP2L=sP2L
                            )

                            signal_df <- as.data.frame(signal_matrix)
                            colnames(signal_df) <- c('signal','sP1','sP2')

                            signal    <- signal_df$signal

                        }

                        signal_a[[length(signal_a)+1]] <- signal

                        s0 <- signal[length(signal)]

                        if (model_selected_sim == 'one_site') {

                            signal <- pykingenie$one_site_dissociation_analytical(np_dissociation_time_shifted,s0,k_off)

                        }

                        if (model_selected_sim == 'one_site_mtl') {

                            signal <- pykingenie$solve_ode_one_site_mass_transport_dissociation(
                            np_dissociation_time,s0,Kd,k_off,k_tr,smax)

                        }

                        if (model_selected_sim == 'heterogeneous_analyte') {

                            # Signal of the first and second populations
                            signal1 <- tail(signal_matrix[1,],1)
                            signal2 <- tail(signal_matrix[2,],1)

                            np_signal_start  <- np_array(c(signal1,signal2))

                            signal_matrix <- pykingenie$solve_ode_mixture_analyte_dissociation(
                                np_dissociation_time,np_signal_start,np_k_offs)

                            signal_matrix <- as.matrix(signal_matrix)

                            signal <- apply(signal_matrix,2,sum)

                        }

                        if (model_selected_sim == 'one_site_induced_fit') {

                            sP2L <- tail(signal_df$sP2L,1)

                            signal_mat <- pykingenie$solve_induced_fit_dissociation(np_dissociation_time, k_off, k_c, k_rev,s0,sP2L,smax)
                            signal_df <- as.data.frame(signal_mat)
                            colnames(signal_df) <- c('signal','sP1L','sP2L')

                            signal    <- signal_df$signal

                        }

                        if (model_selected_sim == 'one_site_conformational_selection') {

                            sP1 <- tail(signal_df$sP1,1)

                            signal_mat <- pykingenie$solve_conformational_selection_dissociation(
                                np_dissociation_time, k_off, k_c, k_rev,smax=smax,sP1=sP1,sP2L=s0
                            )
                            signal_df <- as.data.frame(signal_mat)
                            colnames(signal_df) <- c('signal','sP1','sP2')

                            signal    <- signal_df$signal

                        }

                        signal_d[[length(signal_d)+1]] <- signal

                        if (model_selected_sim %in% c(
                            'one_site_induced_fit',
                            'one_site_conformational_selection')) {

                            signal_d_adv[[length(signal_d_adv)+1]] <- tail(signal_df,1)

                        }

                        if (model_selected_sim == 'heterogeneous_analyte') {

                            signal_d_adv[[length(signal_d_adv)+1]] <- signal_matrix[,ncol(signal_matrix)]

                        }

                    }

                    signal_per_pc_assoc[[length(signal_per_pc_assoc)+1]] <- signal_a
                    signal_per_pc_disso[[length(signal_per_pc_disso)+1]] <- signal_d

                    if (model_selected_sim %in% multi_state_models) {

                        signal_per_pc_disso_adv[[length(signal_per_pc_disso_adv)+1]] <- signal_d_adv

                    }

                }

                signal_per_cycle_assoc[[cycle]] <- signal_per_pc_assoc
                signal_per_cycle_disso[[cycle]] <- signal_per_pc_disso

                if (model_selected_sim %in% multi_state_models) {

                    signal_per_cycle_disso_adv[[cycle]] <- signal_per_pc_disso_adv

                }
            }
        }

        reactives$simulation_results <- list(
            'association_time_per_cycle'  = association_time_per_cycle,
            'dissociation_time_per_cycle' = dissociation_time_per_cycle,
            'signal_per_cycle_assoc'      = signal_per_cycle_assoc,
            'signal_per_cycle_disso'      = signal_per_cycle_disso,
            'smax_all'                    = smax_all,
            'lig_concs'                   = lig_concs_init)

        output$signal_sim_plot <- renderPlotly({

            fig <- plot_simulation(association_time_per_cycle,
                                   dissociation_time_per_cycle,
                                   signal_per_cycle_assoc,
                                   signal_per_cycle_disso,
                                   smax_all,
                                   lig_concs_init,
                                   plot_width  = input$plot_width_sim,
                                   plot_height = input$plot_height_sim,
                                   plot_type   = input$plot_type_sim,
                                   font_size   = input$plot_axis_size_sim,
                                   show_grid_x = input$plot_show_grid_x_sim,
                                   show_grid_y = input$plot_show_grid_y_sim,
                                   marker_size = input$plot_marker_size_sim,
                                   is_single_cycle = input$is_single_cycle_sim
                                   )

            return(fig)

        })

    }

    # Surface binding
    if (model_type_sim == "solution") {

        model_is_conf_sel <- model_selected_sim == 'one_site_conformational_selection'
        model_is_if       <- model_selected_sim == 'one_site_induced_fit'

        if (model_is_conf_sel || model_is_if) {

            bound_signal_is_equal <- input$signal_ES_int == input$signal_ES

        } else {

            bound_signal_is_equal <- TRUE

        }

        time_step  <- input$time_step
        total_time <- input$total_time

        k_obs_dominant_per_pc     <- list()
        k_obs_non_dominant_per_pc <- list()

        if (total_time / time_step  > 1e5) {

            popUpWarning(
                'The total number of steps is larger than 10000.
                Please reduce the interaction time or increase the time step.'
            )

            return(NULL)
        }

        # 3. Obtain the association and dissociation times
        interaction_time <- seq(0,total_time,time_step)

        signal_per_pc <- list()

        # Obtain the smax
        prot_concs <- input$protein_conc_sim / (input$prot_dil_factor_sim ^ (0:input$numb_dil_sim_prot))

        for (pc in prot_concs) {

            k_obs_dominant     <- c()
            k_obs_non_dominant <- c()

            signal <- list()

            for (lc in lig_concs) {

                if (model_selected_sim == 'one_site') {

                    signal_single <- pykingenie$signal_ode_one_site_insolution(
                        np_array(interaction_time),k_off,Kd,pc,lc,
                        signal_a=input$signal_E,
                        signal_b=input$signal_S,
                        signal_complex=input$signal_ES_simple
                        )

                }

                if (model_is_if) {

                    # Initial concentration of E·S, and ES
                    y <- list(0,0)

                    signal_single <- pykingenie$signal_ode_induced_fit_insolution(
                        np_array(interaction_time), y, k_on, k_off, k_c, k_rev,
                        E_tot=pc, S_tot=lc,
                        t0=0,signal_E=input$signal_E,signal_S=input$signal_S,
                        signal_ES_int=input$signal_ES_int,signal_ES=input$signal_ES
                    )

                    k_obs <- pykingenie$get_kobs_induced_fit(lc, pc, k_rev, k_c, k_on, k_off)
                    k_obs_dominant <- c(k_obs_dominant,k_obs)

                    k_obs2 <- pykingenie$get_kobs_induced_fit(lc, pc, k_rev, k_c, k_on, k_off,dominant=FALSE)
                    k_obs_non_dominant <- c(k_obs_non_dominant,k_obs2)

                }

                if (model_is_conf_sel) {

                    # Concentration of E1, E2, S, and E2S
                    e_concs <- pykingenie$get_initial_concentration_conformational_selection(pc,k_c,k_rev)

                    y <- list(e_concs[2],0) # Initial concentrations of E2 and E2S

                    signal_single <- pykingenie$signal_ode_conformational_selection_insolution(
                        np_array(interaction_time), y, k_c, k_rev,k_on, k_off,
                        E_tot=pc,
                        S_tot=lc,
                        t0=0,
                        signal_E1=input$signal_E,
                        signal_E2=input$signal_E,
                        signal_S=input$signal_S,
                        signal_E2S=input$signal_ES_simple
                    )

                    k_obs <- pykingenie$get_kobs_conformational_selection(lc, pc, k_rev, k_c, k_on, k_off)
                    k_obs_dominant <- c(k_obs_dominant,k_obs)

                    k_obs2 <- pykingenie$get_kobs_conformational_selection(lc, pc, k_rev, k_c, k_on, k_off,dominant=FALSE)
                    k_obs_non_dominant <- c(k_obs_non_dominant,k_obs2)

                }

                signal[[length(signal)+1]] <- signal_single

            }

            signal_per_pc[[length(signal_per_pc)+1]] <- signal

            if (model_is_if || model_is_conf_sel) {

                k_obs_dominant_per_pc[[length(k_obs_dominant_per_pc)+1]]         <- k_obs_dominant
                k_obs_non_dominant_per_pc[[length(k_obs_non_dominant_per_pc)+1]] <- k_obs_non_dominant

            }

        }

        if (model_is_if || (model_is_conf_sel && bound_signal_is_equal)) {

            # Append plot tab with the relaxation rates
            new_tab <- tabPanel("Relaxation rates", plotlyOutput("plot_relaxation_rates"))
            insertTab(inputId = "tabset_sim", tab = new_tab, target = "Signal", position = "after")

            output$plot_relaxation_rates <- renderPlotly({

                fig <- plot_relaxation_rates(k_obs_dominant_per_pc,
                                             k_obs_non_dominant_per_pc,
                                             prot_concs,
                                             lig_concs,
                                             plot_width  = input$plot_width_sim,
                                             plot_height = input$plot_height_sim,
                                             plot_type   = input$plot_type_sim,
                                             font_size   = input$plot_axis_size_sim,
                                             show_grid_x = input$plot_show_grid_x_sim,
                                             show_grid_y = input$plot_show_grid_y_sim,
                                             marker_size = input$plot_marker_size_sim
                                             )
            })

            reactives$relaxation_plot_done <- TRUE

        }

        output$signal_sim_plot <- renderPlotly({

            fig <- plot_simulation(list(interaction_time), # Convert to list (we only have one cycle here)
                                   NULL, # No dissociation time
                                   list(signal_per_pc), # Convert to list (we only have one cycle here)
                                   NULL, # No dissociation signal
                                   prot_concs,
                                   lig_concs,
                                   plot_width  = input$plot_width_sim,
                                   plot_height = input$plot_height_sim,
                                   plot_type   = input$plot_type_sim,
                                   font_size   = input$plot_axis_size_sim,
                                   show_grid_x = input$plot_show_grid_x_sim,
                                   show_grid_y = input$plot_show_grid_y_sim,
                                   marker_size = input$plot_marker_size_sim
                                   )

            return(fig)

        })

        reactives$simulation_results <- list(
            'interaction_time'  = interaction_time,
            'signal_per_pc'     = signal_per_pc,
            'prot_concs'        = prot_concs,
            'lig_concs'         = lig_concs
            )

    }

    reactives$simulation_results[['model_type']] <- model_type_sim

    reactives$simulation_results[['is_single_cycle']] <- input$is_single_cycle_sim

    reactives$simulation_run <- TRUE

})

observeEvent(list(input$model_type_sim,input$model_selected_sim),{

    model_type <- input$model_type_sim
    model_sel  <- input$model_selected_sim

    if (model_type == "solution") {

        if (model_sel == 'one_site_conformational_selection') {

            #Example values base on Fabian Paul ,Thomas R. Weikl, 2016:
            updateNumericInput(session,'total_time',NULL,value=1)
            updateNumericInput(session,'kon_sim_1to1_adv',NULL,value=100)
            updateNumericInput(session,'koff_sim_1to1_adv',NULL,value=1)
            updateNumericInput(session,'kc_sim_1to1_adv',NULL,value=10)
            updateNumericInput(session,'krev_sim_1to1_adv',NULL,value=100)
            updateNumericInput(session,'protein_conc_sim',NULL,value=1)
            updateNumericInput(session,'init_lig_sim',NULL,value=1.8)
            updateNumericInput(session,'lig_dil_factor_sim',NULL,value=1.6)
            updateNumericInput(session,'time_step',NULL,value=0.02)

        }

        if (model_sel == 'one_site_induced_fit') {

            #Example values based on Fabian Paul ,Thomas R. Weikl, 2016:
            updateNumericInput(session,'total_time',NULL,value=1)
            updateNumericInput(session,'kon_sim_1to1_adv',NULL,value=100)
            updateNumericInput(session,'koff_sim_1to1_adv',NULL,value=100)
            updateNumericInput(session,'kc_sim_1to1_adv',NULL,value=1)
            updateNumericInput(session,'krev_sim_1to1_adv',NULL,value=10)
            updateNumericInput(session,'protein_conc_sim',NULL,value=0.5)
            updateNumericInput(session,'init_lig_sim',NULL,value=1.8)
            updateNumericInput(session,'lig_dil_factor_sim',NULL,value=1.6)
            updateNumericInput(session,'time_step',NULL,value=0.02)

        }

        if (model_sel == 'one_site') {

            #Example values base on Fabian Paul ,Thomas R. Weikl, 2016:
            updateNumericInput(session,'total_time',NULL,value=1)
            updateNumericInput(session,'kd_sim_1to1',NULL,value=0.1)
            updateNumericInput(session,'koff_sim_1to1',NULL,value=0.01)
            updateNumericInput(session,'protein_conc_sim',NULL,value=1)
            updateNumericInput(session,'init_lig_sim',NULL,value=5)
            updateNumericInput(session,'lig_dil_factor_sim',NULL,value=1.6)
            updateNumericInput(session,'time_step',NULL,value=0.02)

        }

    } else {

        updateNumericInput(session,'time_step',NULL,value=0.5)
        updateNumericInput(session,'protein_smax_sim',NULL,value=5)
        updateNumericInput(session,'association_time',NULL,value=300)
        updateNumericInput(session,'dissociation_time',NULL,value=600)
        updateNumericInput(session,'numb_cycles_sim',NULL,value=1)
        updateNumericInput(session,'kd_sim_1to1',NULL,value=0.5)
        updateNumericInput(session,'koff_sim_1to1',NULL,value=0.01)
        updateNumericInput(session,'kon_sim_1to1_adv',NULL,value=0.5)
        updateNumericInput(session,'koff_sim_1to1_adv',NULL,value=0.01)
        updateNumericInput(session,'kc_sim_1to1_adv',NULL,value=1)
        updateNumericInput(session,'krev_sim_1to1_adv',NULL,value=10)

    }

})





