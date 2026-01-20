append_record_to_logbook <- function(record_str,include_time=FALSE,add_empty_line=FALSE) {

    if (include_time) {
        record_str <- paste0(as.character(format(Sys.time(),usetz = TRUE)),' ',record_str)
    }

    if (add_empty_line) record_str <- c('',record_str)

    reactives$logbook <- append(reactives$logbook, record_str)

    return(NULL)
}

updateSelectedExperiment <- function(newExperimentName) {

    possible_exps <- add_all_option(pyKinetics$experiment_names)
    updateSelectInput(session,'selectedExperiment',NULL,possible_exps,newExperimentName)

}

render_legend_df <- function() {

    if (reactives$surface_based_binding) {

        labels <- unlist(pyKinetics$get_experiment_properties('sensor_names'))
        ids    <- unlist(pyKinetics$get_experiment_properties('sensor_names_unique'))

    } else {

        labels <- unlist(pyKinetics$get_experiment_properties('traces_names'))
        ids    <- unlist(pyKinetics$get_experiment_properties('traces_names_unique'))

    }

    output$legendInfo <- render_RHandsontable(get_plotting_df(ids,labels))

    return(NULL)
}

# Load the files in ./www/test_bli_folder/
observeEvent(input$loadExampleData,{

    reactives$traces_loaded         <- FALSE
    reactives$surface_based_binding <- TRUE
    # Clear all experiments
    pyKinetics$delete_experiment(pyKinetics$experiment_names)

    # list files in test_bli_folder
    files <- list.files("./www/test_bli_folder",full.names = TRUE)
    names <- list.files("./www/test_bli_folder")

    newExperimentName <- 'Example'

    bli_experiment <- pykingenie$OctetExperiment(newExperimentName)

    bli_experiment$read_sample_plate_info(files,names)
    bli_experiment$read_sensor_data(files,names)

    pyKinetics$add_experiment(bli_experiment,newExperimentName)

    render_legend_df()

    updateSelectedExperiment(newExperimentName)
    updateSelectInput(session,"experiment2delete",choices     = c("ALL",pyKinetics$experiment_names))

    reactives$traces_loaded         <- TRUE
    render_ligand_info_df()

})

observeEvent(input$kineticFiles, {

    req(input$kineticFiles)

    files     <- input$kineticFiles$datapath

    frd_files <- files[grepl('.frd',files)]
    csv_files <- files[grepl('.csv',files)]

    n_frd <- length(frd_files)
    n_csv <- length(csv_files)

    if (n_frd == 0 && n_csv == 0) {

        popUpWarning('No .frd or .csv files found.')
        return(NULL)

    }

    if (n_frd > 0) {

        newExperimentName <- pykingenie$guess_experiment_name(frd_files[1])

    } else {

        newExperimentName <- 'Experiment'

    }

    showModal(modalDialog(

        tags$h3('Set the experiment name:'),
        textInput("newExperimentName", NULL,newExperimentName),

        footer=tagList(
          actionButton('submitLoadExperiment', 'Submit'),
          modalButton('Cancel')
        )

    ))
})

observeEvent(input$submitLoadExperiment,{

    removeModal()
    traces_loaded <- FALSE

    updateCheckboxInput(session, "showLigandInfo", value = TRUE)

    files <- input$kineticFiles$datapath
    names <- input$kineticFiles$name

    sorted_indices <- order(names)

    files <- files[sorted_indices]
    names <- names[sorted_indices]

    # Check if the type of experiment
    exp_type <- pykingenie$guess_experiment_type(files)

    # Check if it is surface-based
    exp_is_surface_based <- exp_type == 'surface'

    # Check if we have data loaded and if it is of the same type
    if (reactives$traces_loaded) {

        if (exp_is_surface_based && !reactives$surface_based_binding) {

            popUpWarning('Oops! We can not load the file,
            because the type of experiment (surface-based) is different from the one loaded.')
            return(NULL)

        }

        if (!exp_is_surface_based && reactives$surface_based_binding) {

            popUpWarning('Oops! We can not load the file,
            because the type of experiment (solution-based) is different from the one loaded.')
            return(NULL)

        }

    }

    reactives$surface_based_binding <- exp_is_surface_based

    reactives$traces_loaded <- FALSE

    frd_files <- files[grepl('.frd',files)]
    n_frds    <- length(frd_files)

    if (n_frds>0) {

        bli_experiment <- pykingenie$OctetExperiment(input$newExperimentName)
        bli_experiment$read_sample_plate_info(files,names)

        result <- tryCatch(
            {

                bli_experiment$read_sensor_data(files,names)

            }, error = function(e) {

                if (inherits(e, "python.builtin.ValueError")) {
                    err <- py_last_error()
                    popUpWarning(
                        paste0("⚠ Processing error: ", err$value)
                    )
                    return('Error')
                } else {
                    stop(e) # rethrow non-Python errors
                }
            }
        )

        if (!is.null(result)) return(NULL)

        pyKinetics$add_experiment(bli_experiment,input$newExperimentName)

        if (bli_experiment$traces_loaded) {

            render_legend_df()

            updateSelectedExperiment(input$newExperimentName)

        }

        traces_loaded <- any(unlist(pyKinetics$get_experiment_properties('traces_loaded')))

        if (bli_experiment$sample_plate_loaded) {

            append_sample_plate_plot(input$newExperimentName)

        }

    } else if (any(grepl('simulation_KinGenie_',names))) {

        nFiles <- length(files)

        for (i in 1:nFiles) {

            if (nFiles > 1) {
                newExperimentName <- paste0(input$newExperimentName,i)
            } else {
                newExperimentName <- input$newExperimentName
            }

            if (exp_is_surface_based) {

                kingenie_simulation <- pykingenie$KinGenieCsv(newExperimentName)

            } else {

                kingenie_simulation <- pykingenie$KinGenieCsvSolution(newExperimentName)

            }

                kingenie_simulation$read_csv(files[i])
                pyKinetics$add_experiment(kingenie_simulation,newExperimentName)

        }

        render_legend_df()
        updateSelectedExperiment(input$newExperimentName)

        traces_loaded <- TRUE

    } else if (any(grepl('.ini',files))) {

        gt <- pykingenie$GatorExperiment(input$newExperimentName)
        gt$read_all_gator_data(files,names)
        pyKinetics$add_experiment(gt,input$newExperimentName)

        Sys.sleep(0.5)
        render_legend_df()
        updateSelectedExperiment(input$newExperimentName)

        traces_loaded <- any(unlist(pyKinetics$get_experiment_properties('traces_loaded')))

    } else if (all(grepl('.csv',files))) {

        showModal(modalDialog(

            tags$h3('Please select if the data corresponds to surface or in-solution binding:'),
            selectInput("bindingType", NULL,choices = c(
                'Surface-based binding'  = 'surface',
                'In-solution binding'    = 'solution'
            )),

            tags$script(HTML("
              $(document).on('keypress', function(e) {
                if(e.which == 13) {
                  e.preventDefault();
                  $('#importCSVs').click();
                }
              });
            ")),

            footer=tagList(
              actionButton('importCSVs', 'Submit'),
              modalButton('Cancel')
            )

        ))

        return(NULL)

    }

    Sys.sleep(0.6)

    reactives$traces_loaded       <- traces_loaded
    render_ligand_info_df()

    if (!traces_loaded) {

        popUpWarning('No traces were loaded. Please check the input files.')
        return(NULL)
    }

    append_record_to_logbook(paste0('Loading experiment: ',input$newExperimentName),include_time = TRUE,add_empty_line = TRUE)

    # Append files to logbook
    append_record_to_logbook(paste0('Files: ',paste(names,collapse = ', ')),include_time = FALSE)

    updateSelectInput(session,"experiment2delete",choices     = c("ALL",pyKinetics$experiment_names))

})


observeEvent(input$importCSVs,{

    removeModal()

    exp_is_surface_based <- input$bindingType == 'surface'

    reactives$surface_based_binding <- exp_is_surface_based

    files <- input$kineticFiles$datapath
    names <- input$kineticFiles$name

    sorted_indices <- order(names)

    files <- files[sorted_indices]
    names <- names[sorted_indices]

    nFiles <- length(files)

    for (i in 1:nFiles) {

        if (nFiles > 1) {
            newExperimentName <- paste0(input$newExperimentName,i)
        } else {
            newExperimentName <- input$newExperimentName
        }

        if (exp_is_surface_based) {

            kingenie_csv <- pykingenie$KinGenieCsv(newExperimentName)

        } else {

            kingenie_csv <- pykingenie$KinGenieCsvSolution(newExperimentName)

        }

            kingenie_csv$read_csv(files[i])
            pyKinetics$add_experiment(kingenie_csv,newExperimentName)

    }

    pyKinetics$collapse_solution_experiments(input$newExperimentName)

    render_legend_df()
    updateSelectedExperiment(input$newExperimentName)

    traces_loaded <- TRUE

    Sys.sleep(0.6)

    reactives$traces_loaded       <- traces_loaded
    render_ligand_info_df()

    if (!traces_loaded) {

        popUpWarning('No traces were loaded. Please check the input files.')
        return(NULL)
    }

    append_record_to_logbook(paste0('Loading experiment: ',input$newExperimentName),include_time = TRUE,add_empty_line = TRUE)

    # Append files to logbook
    append_record_to_logbook(paste0('Files: ',paste(names,collapse = ', ')),include_time = FALSE)

    updateSelectInput(session,"experiment2delete",choices     = c("ALL",pyKinetics$experiment_names))

})

observeEvent(input$triggerDeletion,{

    req(reactives$traces_loaded)

    experiment2delete <- input$experiment2delete

    if (experiment2delete == 'None') return(NULL)

    reactives$traces_loaded   <- FALSE

    if (experiment2delete == 'ALL') {

        pyKinetics$delete_experiment(pyKinetics$experiment_names)

        output$legendInfo <- NULL

        updateSelectInput(session,"experiment2delete",NULL, 'None')
        updateSelectInput(session,'selectedExperiment',NULL,'None')

        append_record_to_logbook('All experiments were deleted.',include_time = TRUE,add_empty_line = TRUE)

    } else {

        pyKinetics$delete_experiment(experiment2delete)

        append_record_to_logbook(paste0('Experiment ',experiment2delete,' was deleted.'),include_time = TRUE,add_empty_line = TRUE)

        current_experiments <- pyKinetics$experiment_names

        if (length(current_experiments) == 0) {

            updateSelectInput(session,"experiment2delete",NULL, 'None')
            updateSelectInput(session,'selectedExperiment',NULL,'None')
            output$legendInfo <- NULL

        } else {

            updateSelectInput(session,'experiment2delete',NULL,c("ALL",pyKinetics$experiment_names))
            updateSelectInput(session,'selectedExperiment',NULL,current_experiments)

            render_legend_df()
            render_ligand_info_df()

            reactives$traces_loaded   <- TRUE

        }

    }

})

observeEvent(input$legendInfo,{

    req(reactives$traces_loaded)

    ids <- as.character(hot_to_r(input$legendInfo)$Internal_ID)

    updateSelectInput(session,"mol2changeColor",NULL,ids,input$mol2changeColor)

})

observeEvent(input$colorForLegend,{

    req(reactives$traces_loaded)

    legend_df <- hot_to_r(input$legendInfo)

    idx                  <- which(legend_df$Internal_ID == input$mol2changeColor)
    legend_df$Color[idx] <- input$colorForLegend

    output$legendInfo    <- render_RHandsontable(legend_df)

})

# Modify what processing options can be done depending on the selected experiment
observeEvent(input$selectedExperiment,{

    sel_exp <- input$selectedExperiment

    # IF we have a BLI experiments, leave all the options
    req(reactives$traces_loaded)

    if (sel_exp == 'All') {

      types <- unlist(pyKinetics$get_experiment_properties('type'))

      unq_types <- unique(types)

      if (length(unq_types) > 1) {

          popUpWarning('Processing multiple experiments of different types is not supported. Please select a single experiment.')
          return(NULL)

      } else {

          py_type <- types[1]

      }

    } else {

        experiment   <- pyKinetics$experiments[[sel_exp]]
        py_type      <- experiment$type

    }

    if (py_type %in% c('BLI_experiment','Gator_experiment')) {

        choices <-  c(
            'Align association phase' = 'align_association',
            'Inter-step correction (dissociation)'= 'correct_dissociation'
        )

        if (sel_exp != 'All') {

            choices <- c(choices,
                'Subtract baseline'       = 'subtract',
                'Average'                 = 'average'
            )

        }

        # If we have two experiments, add an option to subtract a whole experiment
        n_exps <- length(pyKinetics$experiment_names)

        if (n_exps > 1 && sel_exp != 'All') {
            choices <- c(choices,'Subtract experiment' = 'subtract_experiment')
        }

        if (py_type == 'BLI_experiment' && sel_exp != 'All') {
            choices <-  c(choices,"Merge steps" = 'merge_steps')
        }

    }  else {
        choices <- c('Set time cutoff' = 'set_time_cutoff')
    }

    updateSelectInput(session,'operation',choices = choices)

})

# Modal dialog to ask for the input and desired units
observeEvent(input$triggerProcessing,{

    req(reactives$traces_loaded)

    sel_exp <- input$selectedExperiment

    if (sel_exp != 'All') {
        experiment  <- pyKinetics$experiments[[sel_exp]]
    }

    operation <- input$operation

    if (operation == 'set_time_cutoff') {

        showModal(modalDialog(

            tags$h3('Please select the time threshold:'),
            numericInput("inputTimeCutOff", NULL,value=1,min=0,max=1e6),

            tags$script(HTML("
              $(document).on('keypress', function(e) {
                if(e.which == 13) {
                  e.preventDefault();
                  $('#submitTimeCutOff').click();
                }
              });
            ")),

            footer=tagList(
              actionButton('submitTimeCutOff', 'Submit'),
              modalButton('Cancel')
            )

        ))

    } else {

        # Operations for surface-based experiments require the sensor names
        if (sel_exp != 'All') {
            sensor_names <- experiment$sensor_names
        }

    }

      if (operation == 'subtract_experiment') {

        # Guess which sensor to subtract using the one with the lowest max signal
       all_experiments   <- pyKinetics$experiment_names
       other_experiments <- all_experiments[all_experiments != sel_exp]

        showModal(modalDialog(

            tags$h3('Please select the experiment to use as reference:'),
            selectInput("baselineExperiment", NULL,other_experiments,selectize=FALSE),
            checkboxInput('expSubtractionIsInPlace',"Apply in-place subtraction",TRUE),

            tags$script(HTML("
              $(document).on('keypress', function(e) {
                if(e.which == 13) {
                  e.preventDefault();
                  $('#submitExpSubtraction').click();
                }
              });
            ")),

            footer=tagList(
              actionButton('submitExpSubtraction', 'Submit'),
              modalButton('Cancel')
            )

        ))

    }

    if (operation == 'subtract') {

        # Guess which sensor to subtract using the one with the lowest max signal
        ys  <- experiment$ys
        idx <- find_probable_baseline(ys)

        sel_sensor_name  <- sensor_names[idx]
        sel_sensor_names <- c(sel_sensor_name,sensor_names[sensor_names != sel_sensor_name])

        showModal(modalDialog(

            tags$h3('Please select the baseline:'),
            selectInput("inputBaseline", NULL,sel_sensor_names,selectize=FALSE),

            tags$script(HTML("
              $(document).on('keypress', function(e) {
                if(e.which == 13) {
                  e.preventDefault();
                  $('#submitSubtraction1').click();
                }
              });
            ")),

            footer=tagList(
              actionButton('submitSubtraction1', 'Submit'),
              modalButton('Cancel')
            )

        ))

    }

    if (operation %in% c('average','align_association','correct_dissociation')) {

        if (sel_exp != 'All') {
            output$tableSelection  <- renderRHandsontable({get_rtable_processing(sensor_names)})
        } else {

            dfs <- list()

            for (i in 1:length(pyKinetics$experiment_names)) {

              exp_name     <- pyKinetics$experiment_names[i]
              experiment   <- pyKinetics$experiments[[exp_name]]
              sensor_names <- experiment$sensor_names

              df_temp <- get_sensor_df(sensor_names,exp_name)
              dfs[[i]] <- df_temp

            }

            df <- do.call(rbind,dfs)

            rdf <- rhandsontable(df) %>%
              hot_col('Select') %>%
              hot_table(stretchH='all')  %>%
              hot_col('ID',readOnly=TRUE) %>%
              hot_col('Experiment',readOnly=TRUE)

            output$tableSelection  <- renderRHandsontable({rdf})

        }

    }


    if (operation == 'average') {

        showModal(modalDialog(

            tags$h3('Please select the curves:'),

            tags$script(HTML("
              $(document).on('keypress', function(e) {
                if(e.which == 13) {
                  e.preventDefault();
                  $('#submitAverage').click();
                }
              });
            ")),

            rHandsontableOutput('tableSelection'),

            textInput('outputNameAverage','Output name','Average'),

            footer=tagList(
              actionButton('submitAverage', 'Submit'),
              modalButton('Cancel')
            )
        ))

    }

    if (operation == 'align_association') {

        showModal(modalDialog(

            tags$h3('Please select the curves:'),

            tags$script(HTML("
              $(document).on('keypress', function(e) {
                if(e.which == 13) {
                  e.preventDefault();
                  $('#submitAlign').click();
                }
              });
            ")),

            rHandsontableOutput('tableSelection'),

            checkboxInput('inPlaceAlignment','Perform in-place alignment',TRUE),
            conditionalPanel(
                condition = "input.inPlaceAlignment",
                checkboxInput('createNewSensorNames','Create new sensor names',FALSE),
            ),
            checkboxInput('keepRegeneration','Keep the regeneration data',FALSE),
            checkboxInput('keepLoading','Keep the loading data',FALSE),
            checkboxInput('keepBaseline','Keep the baseline data',FALSE),
            checkboxInput('keepActivation','Keep the activation data',FALSE),
            checkboxInput('keepQuenching','Keep the quenching data',FALSE),
            checkboxInput('keepCustom','Keep the custom step data',FALSE),

            footer=tagList(
              actionButton('submitAlign', 'Submit'),
              modalButton('Cancel')
            )
        ))

    }

    if (operation == 'correct_dissociation') {

        showModal(modalDialog(

            tags$h3('Please select the curves:'),

            tags$script(HTML("
              $(document).on('keypress', function(e) {
                if(e.which == 13) {
                  e.preventDefault();
                  $('#submitCorrectDis').click();
                }
              });
            ")),

            rHandsontableOutput('tableSelection'),

            checkboxInput('inPlaceCorrection','Perform in-place inter-step correction',TRUE),
            conditionalPanel(
                condition = "input.inPlaceCorrection",
                checkboxInput('createNewSensorNames','Create new sensor names',FALSE),
            ),

            footer=tagList(
              actionButton('submitCorrectDis', 'Submit'),
              modalButton('Cancel')
            )
        ))

    }

    if (operation == 'merge_steps') {

        showModal(modalDialog(

            tags$h3('Please choose the mode for merging steps:'),

            tags$script(HTML("
              $(document).on('keypress', function(e) {
                if(e.which == 13) {
                  e.preventDefault();
                  $('#submitMerge').click();
                }
              });
            ")),

            selectInput("merging_mode", NULL,choices = c(
                'By step name'  = 'step_name',
                'By step index' = 'step_index'
            )),

            footer=tagList(
              actionButton('submitMerge', 'Submit'),
              modalButton('Cancel')
            )
        ))

    }

})

observeEvent(input$submitTimeCutOff,{

    removeModal()

    reactives$traces_loaded <- FALSE

    exp <- pyKinetics$experiments[[input$selectedExperiment]]

    exp$cut_off_time(input$inputTimeCutOff)

    Sys.sleep(0.6)

    # Include the time cutoff step in the logbook
    append_record_to_logbook(paste0('Time cutoff set to ',input$inputTimeCutOff,' seconds.'),include_time = TRUE,add_empty_line = TRUE)

    reactives$traces_loaded <- TRUE

})

observeEvent(input$submitExpSubtraction,{

    removeModal()
    reactives$traces_loaded <- FALSE

    exp <- pyKinetics$experiments[[input$selectedExperiment]]

    other_exp <- pyKinetics$experiments[[input$baselineExperiment]]

    result <- tryCatch(
        {

            exp$subtract_experiment(other_exp,inplace = input$expSubtractionIsInPlace)

        }, error = function(e) {

            if (inherits(e, "python.builtin.RuntimeError")) {
                err <- py_last_error()
                popUpWarning(
                    paste0("⚠ Processing error: ", err$value)
                )
                return('Error')
            } else {
                stop(e) # rethrow non-Python errors
            }
        }
    )

    if (!is.null(result)) return(NULL)

    append_record_to_logbook('Experiment subtraction performed',include_time = TRUE,add_empty_line = TRUE)
    append_record_to_logbook(paste0('Experiment ',input$baselineExperiment,' subtracted from ',input$selectedExperiment,'.'),include_time = FALSE, add_empty_line = FALSE)

    reactives$traces_loaded <- TRUE

    Sys.sleep(0.6)

    render_legend_df()
    render_ligand_info_df()

})

observeEvent(input$submitMerge,{

    removeModal()

    merging_mode <- input$merging_mode

    exp <- pyKinetics$experiments[[input$selectedExperiment]]

    df_steps <- exp$df_steps

    step_names <- unique(df_steps$Name)

    if (merging_mode == "step_name") {

        showModal(modalDialog(

            tags$h3('Please choose the step name.'),
            tags$h4('All steps with the specified name will be merged with the step that immediately follows them.
            To view the outcome of this operation, please check the "Steps" tab.
            This is useful for correcting mislabeled steps—for example, if a "Baseline" step appears after a "Dissociation" step but
            actually represents a continuation of the dissociation phase.'),

            tags$script(HTML("
              $(document).on('keypress', function(e) {
                if(e.which == 13) {
                  e.preventDefault();
                  $('#submitMerge_2').click();
                }
              });
            ")),

            selectInput("step_name", NULL,choices = step_names),

            footer=tagList(
              actionButton('submitMerge_2', 'Submit'),
              modalButton('Cancel')
            )
        ))

    }

   if (merging_mode == "step_index") {

        showModal(modalDialog(

            tags$h3('Please choose the index of the reference step and the index of the step to combine.'),
            tags$h4('As a result, the x and y values of both steps will be combined into a single step.
            The analyte concentration and column location will be inherited from the reference step.
            Note: The steps to be merged must have consecutive indexes.'),

            tags$script(HTML("
              $(document).on('keypress', function(e) {
                if(e.which == 13) {
                  e.preventDefault();
                  $('#submitMerge_3').click();
                }
              });
            ")),

            numericInput("step_index_ref", "Reference step index", value = 1, min = 1, max = nrow(df_steps)),
            numericInput("step_index_to_merge", "Step index to merge", value = 2, min = 1, max = nrow(df_steps)),

            footer=tagList(
              actionButton('submitMerge_3', 'Submit'),
              modalButton('Cancel')
            )
        ))
    }
})

observeEvent(input$submitMerge_3,{

    removeModal()

    reactives$traces_loaded <- FALSE

    exp <- pyKinetics$experiments[[input$selectedExperiment]]

    exp$merge_consecutive_steps(
        idx_ref = input$step_index_ref,
        idx_to_merge = input$step_index_to_merge
    )

    reactives$traces_loaded <- TRUE

    append_record_to_logbook('Merging of steps by index performed.',include_time = TRUE,add_empty_line = TRUE)
    # append the selecte experiment, reference index and index to merge to the logbook
    append_record_to_logbook(
        paste0(
            'Experiment: ',input$selectedExperiment,
            ', Reference step index: ',input$step_index_ref,
            ', Step index to merge: ',input$step_index_to_merge
        )
    )

})


observeEvent(input$submitMerge_2,{

    removeModal()

    reactives$traces_loaded <- FALSE

    exp <- pyKinetics$experiments[[input$selectedExperiment]]

    exp$merge_consecutive_steps_by_name(step_name=input$step_name)

    reactives$traces_loaded <- TRUE

    append_record_to_logbook('Merging of steps by name performed.',include_time = TRUE,add_empty_line = TRUE)
    append_record_to_logbook(
        paste0(
            'Experiment: ',input$selectedExperiment,
            ', Step name: ',input$step_name
        )
    )

})

observeEvent(input$submitCorrectDis,{

    removeModal()

    tableCorrectDis <- hot_to_r(input$tableSelection)
    tableCorrectDis <- tableCorrectDis[tableCorrectDis$Select,]

    reactives$traces_loaded <- FALSE

    exps_to_analyse <- list()
    sensor_to_correct <- list()

    if (input$selectedExperiment == 'All') {

        unq_exps <- unique(tableCorrectDis$Experiment)

        for (exp in unq_exps) {

            temp_df <- tableCorrectDis[tableCorrectDis$Experiment == exp, ]
            samples <- c(temp_df$ID)

            exps_to_analyse[[length(exps_to_analyse)+1]] <- exp
            sensor_to_correct[[length(sensor_to_correct)+1]] <- samples

        }

    } else {

        exps_to_analyse[[1]] <- input$selectedExperiment
        sensor_to_correct[[1]]  <- c(tableCorrectDis$ID)

    }

    for (exp in exps_to_analyse) {

        samples <- sensor_to_correct[[which(exps_to_analyse == exp)]]

        py_exp <- pyKinetics$experiments[[exp]]

        py_exp$align_dissociation(samples,input$inPlaceCorrection,input$createNewSensorNames)

    }

    Sys.sleep(0.6)

    render_legend_df()
    render_ligand_info_df()
    # Include the correction step in the logbook

    append_record_to_logbook('Inter-step correction performed (between dissociation and association).',
    include_time = TRUE,add_empty_line = TRUE)

    append_record_to_logbook(paste0('Samples ',paste(samples,collapse = ', '),' corrected.'))

    # Append in-place option to the logbook
    append_record_to_logbook(paste0('In-place correction:',input$inPlaceCorrection))

    reactives$traces_loaded <- TRUE

})

# Modal dialog to ask for the input and desired units
observeEvent(input$submitAlign,{

    removeModal()

    tableAlign     <- hot_to_r(input$tableSelection)
    tableAlign     <- tableAlign[tableAlign$Select,]

    reactives$traces_loaded <- FALSE

    sel_exp <- input$selectedExperiment

    exps_to_analyse <- list()
    sensor_to_align <- list()

    if (sel_exp == 'All') {

        unq_exps <- unique(tableAlign$Experiment)

        for (exp in unq_exps) {

            temp_df <- tableAlign[tableAlign$Experiment == exp, ]
            samples <- c(temp_df$ID)

            exps_to_analyse[[length(exps_to_analyse)+1]] <- exp
            sensor_to_align[[length(sensor_to_align)+1]] <- samples

        }

    } else {

        exps_to_analyse[[1]] <- sel_exp
        sensor_to_align[[1]] <- c(tableAlign$ID)

    }


    for (exp in exps_to_analyse) {

        samples <- sensor_to_align[[which(exps_to_analyse == exp)]]

        py_exp <- pyKinetics$experiments[[exp]]

        py_exp$align_association(samples,input$inPlaceAlignment,input$createNewSensorNames)

        # Find new generated names
        all_names <- unlist(py_exp$sensor_names)
        new_names <- all_names[!(all_names %in% samples)]

        # Include the previous names, if we did the alignment it place
        if (input$inPlaceAlignment) new_names <- c(new_names,samples)

        if (!input$keepRegeneration)  py_exp$discard_steps(new_names,list('KREGENERATION'))
        if (!input$keepLoading)       py_exp$discard_steps(new_names,list('LOADING'))
        if (!input$keepBaseline)      py_exp$discard_steps(new_names,list('BASELINE'))
        if (!input$keepActivation)    py_exp$discard_steps(new_names,list('ACTIVATION'))
        if (!input$keepQuenching)     py_exp$discard_steps(new_names,list('QUENCHING'))
        if (!input$keepCustom)        py_exp$discard_steps(new_names,list('CUSTOM'))

    }

    Sys.sleep(0.6)
    render_legend_df()
    render_ligand_info_df()

    # Include the alignment step in the logbook
    append_record_to_logbook('Alignment step performed.',include_time = TRUE,add_empty_line = TRUE)
    append_record_to_logbook(paste0('Samples ',paste(samples,collapse = ', '),' aligned.'))

    # Append the discarded steps to the logbook
    if (!input$keepRegeneration) append_record_to_logbook('Regeneration steps removed.')
    if (!input$keepLoading)      append_record_to_logbook('Loading steps removed.')
    if (!input$keepBaseline)     append_record_to_logbook('Baseline steps removed.')
    if (!input$keepActivation)   append_record_to_logbook('Activation steps removed.')
    if (!input$keepQuenching)    append_record_to_logbook('Quenching steps removed.')
    if (!input$keepCustom)       append_record_to_logbook('Custom steps removed.')

    # Append in-place option to the logbook
    append_record_to_logbook(paste0('In-place alignment:',input$inPlaceAlignment))

    reactives$traces_loaded <- TRUE

})


# Modal dialog to ask for the input and desired units
observeEvent(input$submitAverage,{

    removeModal()

    tableAverage     <- hot_to_r(input$tableSelection)
    samples          <- c(tableAverage$ID[tableAverage$Select])

    if (length(samples) < 2) return(NULL)

    reactives$traces_loaded <- FALSE

    exp <- pyKinetics$experiments[[input$selectedExperiment]]

    exp$average(samples,input$outputNameAverage)

    labels <- unlist(pyKinetics$get_experiment_properties('sensor_names'))
    ids    <- unlist(pyKinetics$get_experiment_properties('sensor_names_unique'))

    output$legendInfo <- render_RHandsontable(get_plotting_df(ids,labels))

    Sys.sleep(0.6)
    render_legend_df()
    render_ligand_info_df()

    # Include the average step in the logbook
    append_record_to_logbook('Average step performed.',include_time = TRUE,add_empty_line = TRUE)
    append_record_to_logbook(paste0('Samples ',paste(samples,collapse = ', '),' averaged.'))

    reactives$traces_loaded <- TRUE
})

# Modal dialog to ask for the input and desired units
observeEvent(input$submitSubtraction1,{

    exp   <- pyKinetics$experiments[[input$selectedExperiment]]
    names <- exp$sensor_names
    names <- names[names != input$inputBaseline]

    removeModal()

    showModal(modalDialog(

        tags$h3('Please select the sample(s):'),

        tags$script(HTML("
          $(document).on('keypress', function(e) {
            if(e.which == 13) {
              e.preventDefault();
              $('#submitSubtraction2').click();
            }
          });
        ")),

        rHandsontableOutput('tableSubtraction'),

        checkboxInput('inPlaceSubtraction','Perform in-place subtraction',TRUE),

        footer=tagList(
          actionButton('submitSubtraction2', 'Submit'),
          modalButton('Cancel')
        )
    ))

})

observeEvent(input$submitSubtraction2,{

    removeModal()

    reactives$traces_loaded <- FALSE

    tableSubtraction <- hot_to_r(input$tableSubtraction)
    samples          <- c(tableSubtraction$ID[tableSubtraction$Select])

    exp <- pyKinetics$experiments[[input$selectedExperiment]]

    exp$subtraction(as.list(samples),input$inputBaseline,input$inPlaceSubtraction)

    labels <- unlist(pyKinetics$get_experiment_properties('sensor_names'))
    ids    <- unlist(pyKinetics$get_experiment_properties('sensor_names_unique'))

    output$legendInfo <- render_RHandsontable(get_plotting_df(ids,labels))

    Sys.sleep(0.6)

    render_legend_df()
    render_ligand_info_df()

    # Include the subtraction step in the logbook
    append_record_to_logbook('Subtraction step performed.',include_time = TRUE,add_empty_line = TRUE)
    append_record_to_logbook(paste0('Samples ',paste(samples,collapse = ', '),' subtracted from ',input$inputBaseline))

    # In-place subtraction
    append_record_to_logbook(paste0('In-place subtraction:',input$inPlaceSubtraction))

    reactives$traces_loaded <- TRUE

})

