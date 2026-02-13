append_sample_plate_plot <- function(experiment_name) {

    plot_output <- paste0('Plate_',experiment_name)

    # Add plot tab to the tabBoxImport tabBox
    appendTab(inputId = "tabBoxImport",
              tabPanel(plot_output,plotlyOutput(plot_output))
              )

    output[[plot_output]] <- renderPlotly({

        pySingleExp <- pyKinetics$experiments[[experiment_name]]

        sample_row    <- pySingleExp$sample_row
        sample_column <- pySingleExp$sample_column
        sample_type   <- pySingleExp$sample_type
        sample_id     <- pySingleExp$sample_id

        sample_conc_labeled <- pySingleExp$sample_conc_labeled

        fig <- plot_plate_info(
            sample_row,sample_column,sample_type,
            sample_id,sample_conc_labeled,
            experiment_name,
            plot_config = reactives$plot_config
        )

        return(fig)

    })

}

output$traces <- renderPlotly({

    req(reactives$traces_loaded)
    req(input$legendInfo)
    req(reactives$surface_based_binding)

    legend_df <- hot_to_r(input$legendInfo)

    xs_all <- lapply(pyKinetics$experiments,function(x) x$xs)
    ys_all <- lapply(pyKinetics$experiments,function(x) x$ys)

    fig <- plot_traces_all(
        xs_all,
        ys_all,
        legend_df$Legend,
        legend_df$Color,
        legend_df$Show,
        plot_config = reactives$plot_config
    )

    return(fig)

})

output$tracesInSolution <- renderPlotly({

    req(reactives$traces_loaded)
    req(input$legendInfo)
    req(!reactives$surface_based_binding)

    legend_df <- hot_to_r(input$legendInfo)

    xs_all <- lapply(pyKinetics$experiments,function(x) x$xs)
    ys_all <- lapply(pyKinetics$experiments,function(x) x$ys)

    fig <- plot_traces_all(
        xs_all,
        ys_all,
        legend_df$Legend,
        legend_df$Color,
        legend_df$Show,
        plot_config = reactives$plot_config
    )

    return(fig)

})

output$tracesAssDiss <- renderPlotly({

    req(reactives$traces_loaded)
    req(input$legendInfo)
    req(reactives$surface_based_binding)

    legend_df <- hot_to_r(input$legendInfo)

    xs_all <- list()
    ys_all <- list()

    cnt <- 1
    for (exp in pyKinetics$experiments) {

        py_df <- exp$df_steps
        df <- pandas_to_r(py_df)

        selected <- grepl('ASS',df$Type)

        xs <- exp$xs
        ys <- exp$ys

        xs_subset <- lapply(xs,function(x) x[selected])
        ys_subset <- lapply(ys,function(y) y[selected])

        xs_all[[cnt]] <- xs_subset
        ys_all[[cnt]] <- ys_subset

        cnt <- cnt + 1

    }

    fig <- plot_traces_all(
        xs_all,
        ys_all,
        legend_df$Legend,
        legend_df$Color,
        legend_df$Show,
        plot_config = reactives$plot_config
    )

    return(fig)

})

output$steady_state <- renderPlotly({

    req(reactives$fit_dataset_loaded)
    req(reactives$surface_based_binding)

    fig <- plot_steady_state(
        pyKinetics$fittings,
        plot_config = reactives$plot_config,
        plot_fit    = reactives$ss_fit_done
    )

    return(fig)

})

output$diagnostic_plot <- renderPlotly({

    req(reactives$traces_loaded)
    req(input$legendInfo)
    req(reactives$surface_based_binding)
    req(reactives$diagnostic_plots_done)

    experiments <- pyKinetics$fittings

    k_obs       <- lapply(experiments,function(x) unlist(x$k_obs))
    ligand_conc <- lapply(experiments,function(x) unlist(x$lig_conc_lst))

    dfs <- list()

    for (i in seq(1,length(experiments))) {

        df <- data.frame(
            'ligand_conc'     = ligand_conc[[i]],
            'k_obs'           = k_obs[[i]],
            'experiment_name' = names(experiments)[i]
        )

        dfs[[i]] <- df

    }

    df <- do.call(rbind,dfs)

    df$experiment_name <- factor(df$experiment_name,levels = names(experiments))

    fig <- diagnostic_plot(
        df,
        plot_config = reactives$plot_config
    )

    return(fig)

})

output$tracesAssDissFit <- renderPlotly({

    req(reactives$fit_dataset_loaded)
    req(reactives$surface_based_binding)

    fig <- plot_association_dissociation(
        pyKinetics$fittings,
        plot_config = reactives$plot_config,
        plot_assoc  = TRUE,
        plot_disso  = TRUE,
        plot_fit    = reactives$kinetics_fit_done
    )

    return(fig)

})

output$residuals <- renderPlotly({

    req(reactives$fit_dataset_loaded)
    req(reactives$surface_based_binding)

    fig <- NULL

    if (reactives$kinetics_fit_done) {

        fig <- plot_association_dissociation_residuals(

            pyKinetics$fittings,
            plot_config = reactives$plot_config,
            plot_assoc  = TRUE,
            plot_disso  = TRUE,
            plot_fit    = reactives$kinetics_fit_done
        )
    }

    if (reactives$ss_fit_done) {

        fig <- plot_steady_state_residuals(
            pyKinetics$fittings,
            plot_config = reactives$plot_config,
            plot_fit    = reactives$ss_fit_done
        )
    }

    return(fig)

})

output$tracesFitSolution <- renderPlotly({

    req(!reactives$surface_based_binding)
    req(reactives$fit_dataset_loaded)

    fig <- plot_interactions(
        pyKinetics$fittings,
        plot_config = reactives$plot_config,
        plot_fit    = reactives$kinetics_fit_done
    )

    return(fig)

})


output$residualsSolution <- renderPlotly({

    req(!reactives$surface_based_binding)
    req(reactives$fit_dataset_loaded)
    req(reactives$kinetics_fit_done)

    fig <- plot_interactions_residuals(
        pyKinetics$fittings,
        plot_config = reactives$plot_config
    )

    return(fig)

})

output$kobs_plot <- renderPlotly({

    req(reactives$traces_loaded)
    req(input$legendInfo)
    req(reactives$kinetics_fit_done)
    req(!reactives$surface_based_binding)

    # Require that we fitted single or double exponentials
    req(reactives$solution_model == "double")

    model <- reactives$solution_model

    fittings <- pyKinetics$fittings

    k_obs_dominant_per_pc_lst     <- lapply(fittings,function(f) (f$k_obs_1_per_prot))
    k_obs_non_dominant_per_pc_lst <- lapply(fittings,function(f) (f$k_obs_2_per_prot))
    protein_conc_lst              <- lapply(fittings,function(f) (f$unq_prot_conc))
    ligand_conc_per_pc_lst        <- lapply(fittings,function(f) (f$lig_conc_per_protein))

    fig <- plot_many_relaxation_rates(
        k_obs_dominant_per_pc_lst,
        k_obs_non_dominant_per_pc_lst,
        protein_conc_lst,
        ligand_conc_per_pc_lst,
        plot_config = reactives$plot_config
    )

    return(fig)

})