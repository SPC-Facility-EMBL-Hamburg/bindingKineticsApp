output$download_log_book <-   downloadHandler(

    filename = function() {
        paste0("logbook_KinGenie_",Sys.Date(),".txt")
    },

    content  = function(file) {

        lines <- trimws(unlist(reactives$logbook))

        lines <- paste(lines, collapse = "\n")

        cat(lines,file = file, sep = "\n")
    }

)

output$download_curves <- downloadHandler(

    filename = function() {
        paste0("raw_curves_KinGenie_",Sys.Date(),".csv")
    },

    content = function(file) {

        req(reactives$fit_dataset_loaded)

        dfs <- list()

        for (fit in pyKinetics$fittings) {

            lig_conc <- fit$lig_conc

            for (i in 1:length(lig_conc)) {

                l   <- as.numeric(lig_conc[[i]])

                x1   <- fit$time_assoc[[i]]
                y1   <- fit$assoc[[i]]

                x2   <- fit$time_disso[[i]]
                y2   <- fit$disso[[i]]

                for (j in 1:ncol(y1)) {

                    df1 <- data.frame(
                        'Time' = x1[,j],
                        'Signal' = y1[,j],
                        'Analyte_concentration_micromolar_constant' = l[j],
                        'Type' = 'Association',
                        'ID'   = fit$names[i]
                    )

                    df2 <- data.frame(
                        'Time' = x2[,j],
                        'Signal' = y2[,j],
                        'Analyte_concentration_micromolar_constant' = 0,
                        'Type' = 'Dissociation',
                        'ID'   = fit$names[i]
                    )

                    dfs[[length(dfs)+1]] <- df1
                    dfs[[length(dfs)+1]] <- df2

                }
            }

        }

        df <- do.call(rbind,dfs)

        write.csv(df,file = file,row.names = FALSE)

    }
)

output$download_fitted_curves <- downloadHandler(

    filename = function() {
        paste0("fitted_curves_KinGenie_",Sys.Date(),".csv")
    },

    content = function(file) {

        req(reactives$fit_dataset_loaded)
        req(reactives$kinetics_fit_done)

        dfs <- list()

        for (fit in pyKinetics$fittings) {

            # Counter that iterates over the fits, which are stored as list of lists
            fc_counter <- 1

            lig_conc <- fit$lig_conc

            fitted_curves_assoc         <- fit$signal_assoc_fit
            fitted_curves_disso         <- fit$signal_disso_fit

            we_have_fitted_assoc <- !is.null(fitted_curves_assoc)
            we_have_fitted_disso <- !is.null(fitted_curves_disso)

            # if we don't have fitted data, skip to next iteration
            if (!we_have_fitted_assoc & !we_have_fitted_disso) next

            for (i in 1:length(lig_conc)) {

                l   <- as.numeric(lig_conc[[i]])

                x1   <- fit$time_assoc[[i]]
                x2   <- fit$time_disso[[i]]

                n_traces <- ncol(x1)

                for (j in 1:n_traces) {

                    if (we_have_fitted_assoc) {

                        df1 <- data.frame(
                            'Time' = x1[,j],
                            'Fitted_signal' = fitted_curves_assoc[[fc_counter]],
                            'Analyte_concentration_micromolar_constant' = l[j],
                            'Type' = 'Association',
                            'ID'   = fit$names[i]
                        )

                        dfs[[length(dfs)+1]] <- df1

                    }

                    if (we_have_fitted_disso) {

                        df2 <- data.frame(
                            'Time' = x2[,j],
                            'Fitted_signal' = fitted_curves_disso[[fc_counter]],
                            'Analyte_concentration_micromolar_constant' = 0,
                            'Type' = 'Dissociation',
                            'ID'   = fit$names[i]
                        )

                        dfs[[length(dfs)+1]] <- df2

                    }

                    # increase the fc_counter, because we have fitted data
                    fc_counter <- fc_counter + 1

                }
            }
        }

        df <- do.call(rbind,dfs)

        write.csv(df,file = file,row.names = FALSE)

    }
)

output$btn_export_simulation <-   downloadHandler(

    filename = function() {
        paste0("simulation_KinGenie_",Sys.Date(),".csv")

    },

    content  = function(file) {

        req(reactives$simulation_run)

        sim_results <- reactives$simulation_results

        sim_model <- sim_results[['model_type']]

        dfs <- list()

        if (sim_model == "solution") {

            interaction_time <- sim_results[['interaction_time']]
            signal_per_pc    <- sim_results[['signal_per_pc']]
            prot_concs       <- sim_results[['prot_concs']]
            lig_concs        <- sim_results[['lig_concs']]

            for (i in 1:length(prot_concs)) {
                for (j in 1:length(lig_concs)) {
                    df <- data.frame(
                        'Time' = interaction_time,
                        'Signal' = signal_per_pc[[i]][[j]],
                        'Protein_concentration_micromolar' = prot_concs[i],
                        'Ligand_concentration_micromolar' = lig_concs[j]
                    )

                    dfs[[length(dfs)+1]] <- df
                }
            }

        } else {

            is_single_cycle            <- sim_results[['is_single_cycle']]

            association_time_per_cycle  <- sim_results[['association_time_per_cycle']]
            dissociation_time_per_cycle <- sim_results[['dissociation_time_per_cycle']]
            signal_per_cycle_assoc      <- sim_results[['signal_per_cycle_assoc']]
            signal_per_cycle_disso      <- sim_results[['signal_per_cycle_disso']]
            smax_all                    <- sim_results[['smax_all']]
            lig_concs                   <- sim_results[['lig_concs']]

            # Option 1 - multi-cycle kinetics, all ligand concentration measured in parallel
            if (!is_single_cycle) {

                for (i in 1:length(smax_all)) {

                    for (j in 1:length(lig_concs)) {

                        df1 <- data.frame(
                            'Time'   = association_time_per_cycle[[1]],
                            'Signal' = signal_per_cycle_assoc[[1]][[i]][[j]],
                            'Smax'   = smax_all[[i]],
                            'Analyte_concentration_micromolar_constant' = lig_concs[j]
                        )

                        df2 <- data.frame(
                            'Time'   = dissociation_time_per_cycle[[1]],
                            'Signal' = signal_per_cycle_disso[[1]][[i]][[j]],
                            'Smax'   = smax_all[[i]],
                            'Analyte_concentration_micromolar_constant' = 0
                        )

                        dfs[[length(dfs)+1]] <- df1
                        dfs[[length(dfs)+1]] <- df2

                    }
                }

            # Option 2 - single-cycle kinetics, all ligand concentration measured one after another
            } else {

                for (cycle in 1:length(association_time_per_cycle)) {

                    for (i in 1:length(smax_all)) {

                        df1 <- data.frame(
                            'Time'   = association_time_per_cycle[[cycle]],
                            'Signal' = signal_per_cycle_assoc[[cycle]][[i]][[1]],
                            'Smax'   = smax_all[[i]],
                            'Analyte_concentration_micromolar_constant' = lig_concs[cycle],
                            'Cycle'  = cycle
                        )

                        df2 <- data.frame(
                            'Time'   = dissociation_time_per_cycle[[cycle]],
                            'Signal' = signal_per_cycle_disso[[cycle]][[i]][[1]],
                            'Smax'   = smax_all[[i]],
                            'Analyte_concentration_micromolar_constant' = 0,
                            'Cycle'  = cycle
                        )

                        dfs[[length(dfs)+1]] <- df1
                        dfs[[length(dfs)+1]] <- df2

                    }
                }

            }
        }

        df <- do.call(rbind,dfs)
        write.csv(df,file = file,row.names = FALSE)

    }
)



