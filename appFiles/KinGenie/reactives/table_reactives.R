render_ligand_info_df <- function(){

    output$ligandInfo <- renderRHandsontable({

        if (reactives$surface_based_binding) {

            pyKinetics$merge_ligand_conc_df()

            py_df <- pyKinetics$combined_ligand_conc_df

            df <- pandas_to_r(py_df)

            return(render_combined_ligand_conc_df(df))

        } else {

            pyKinetics$merge_conc_df_solution()

            py_df <- pyKinetics$combined_conc_df
            df <- pandas_to_r(py_df)

            return(render_combined_ligand_conc_df_solution(df))

        }

    })

    return(NULL)

}

output$stepsInfo <- renderDT({

    req(reactives$traces_loaded)
    req(reactives$surface_based_binding)

    n <- length(pyKinetics$experiment_names)

    dfs <- lapply(pyKinetics$experiments, function(exp) {

        py_df <- exp$df_steps
        df <- pandas_to_r(py_df)

        if (n > 1) df['Name'] <- exp$name

        return(df)

    })

    df <- do.call(rbind,dfs)

    DT::datatable(df, options = list(scrollY = '400px', paging = FALSE), rownames = FALSE)

})

output$tableSubtraction <- renderRHandsontable({

    req(reactives$traces_loaded)
    req(input$inputBaseline)

    names <- pyKinetics$experiments[[input$selectedExperiment]]$sensor_names

    rdf   <- get_rtable_processing( names[names != input$inputBaseline] )

    return(rdf)

})

output$fittingInfoSS <- renderTable({

    req(reactives$ss_fit_done)
    req(reactives$fit_dataset_loaded)

    dfs <- list()

    for (name in pyKinetics$fittings_names) {

        py_df <- pyKinetics$fittings[[name]]$fit_params_ss
        df <- pandas_to_r(py_df)
        dfs[[length(dfs)+1]] <- df

    }

    return(do.call(rbind,dfs))

},digits = 5)

output$fittingInfoKinetics <- renderTable({

    req(reactives$kinetics_fit_done)
    req(reactives$fit_dataset_loaded)

    dfs <- list()

    for (name in pyKinetics$fittings_names) {

        py_df <- pyKinetics$fittings[[name]]$fit_params_kinetics
        df <- pandas_to_r(py_df)
        df$Name <- name

        dfs[[length(dfs)+1]] <- df

    }

    return(do.call(rbind,dfs))

},digits = 5)


output$fittedParamsBoundaries <- renderDT({

    req(reactives$fit_dataset_loaded)
    # Ask either for reactives$kinetics_fit_done or reactives$ss_fit_done
    req(reactives$kinetics_fit_done || reactives$ss_fit_done)

    dfs <- list()

    for (name in pyKinetics$fittings_names) {

        py_df <- pyKinetics$fittings[[name]]$fitted_params_boundaries
        df <- pandas_to_r(py_df)
        df$Name <- name

        dfs[[length(dfs)+1]] <- df

    }

    # Combine all data frames into one
    # The first column is the fitted parameter, the second is the lower boundary, and the third is the upper boundary

    df <- do.call(rbind,dfs)

    # Round to the first three significative digits (of the first three columns)
    df[,1:3] <- signif(df[,1:3], 3)

    # Find if the difference between the fitted and the boundaries is too small
    close_values1 <- abs(df[,2] - df[,1]) / df[,1] < 0.05
    close_values2 <- abs(df[,3] - df[,1]) / df[,3] < 0.05

    # Check if the fitted value is too close to the lower or upper boundary
    df$highlight <- !(close_values1 | close_values2)

    df <-  datatable(
        df,
        options = list(
            columnDefs = list(
                list(visible = FALSE, targets = 5)  # hide the helper column if needed
        ))) %>%
        formatStyle(
            'Fitted_parameter_value',
            backgroundColor = styleEqual(c(TRUE, FALSE), c("lightgreen", "red")),
            valueColumns = 'highlight'
        )

    return(df)

})

output$fittingInfoKineticsCI95 <- renderTable({

    req(reactives$kinetics_fit_done)
    req(reactives$fit_dataset_loaded)

    dfs <- list()

    for (name in pyKinetics$fittings_names) {

        py_df <- pyKinetics$fittings[[name]]$fit_params_kinetics_ci95
        df <- pandas_to_r(py_df)
        df$Name <- name

        dfs[[length(dfs)+1]] <- df

    }

    return(do.call(rbind,dfs))

},digits = 5)
