box::use(
  .. / .. / helpers[
    add_all_option,
    find_probable_baseline,
    pandas_to_r
  ],
  .. / .. / tables[
    get_plotting_df,
    get_rtable_processing
  ],
  reticulate[
    py_last_error
  ],
  rhandsontable[
    hot_to_r,
    renderRHandsontable,
    rhandsontable,
    rHandsontableOutput
  ],
  shiny[
    actionButton,
    checkboxInput,
    column,
    conditionalPanel,
    fluidRow,
    HTML,
    icon,
    modalButton,
    modalDialog,
    moduleServer,
    NS,
    numericInput,
    observeEvent,
    p,
    reactiveVal,
    removeModal,
    req,
    selectInput,
    showModal,
    sliderInput,
    span,
    tagList,
    tags,
    textInput,
    updateSelectInput
  ],
  shinydashboard[
    box
  ]
)

#' @export
processingUI <- function(id) {
  ns <- NS(id)

  tagList(
    box(
      title = "2. Processing", width = 12, solidHeader = T, status = "primary",
      fluidRow(
        column(6, p(
          HTML("<b>Experiment name</b>"),
          span(shiny::icon("info-circle"), id = "info_uu-selectedExperiment"),
          selectInput(ns("selectedExperiment"), NULL,
            c("None" = "none"),
            selectize = FALSE
          ),
          tippy::tippy_this(
            elementId = "info_uu-selectedExperiment",
            tooltip = "Select the experiment to process.
                            ", placement = "right"
          )
        )),
        column(6, p(
          HTML("<b>Operation</b>"),
          span(shiny::icon("info-circle"), id = "info_uu-inputOperation"),
          selectInput(ns("operation"), NULL,
            c(
              "Subtract baseline" = "subtract",
              "Align association phase" = "align_association",
              "Inter-step correction (dissociation)" = "correct_dissociation",
              "Average" = "average",
              "Merge steps" = "merge_steps"
            ),
            selectize = FALSE
          ),
          #' Smooth'                  = 'smooth')),
          tippy::tippy_this(
            elementId = "info_uu-inputOperation",
            tooltip = "Select and press 'Apply'.
                            ", placement = "right"
          )
        ))
      ),
      fluidRow(
        column(4, p(
          HTML("<b><br></b>"),
          actionButton(
            inputId = ns("triggerProcessing"), label = "Apply",
            icon("forward"),
            style = "color: #fff; background-color: #337ab7;
                    border-color: #2e6da4"
          )
        ))
      )
    )
  )
}

#' @export
processingServer <- function(id, state, pyKinetics, legend_df, logbook) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(state$selected_experiment,
      {
        req(state$traces_loaded)
        req(state$selected_experiment)

        updateSelectInput(
          session,
          "selectedExperiment",
          choices = add_all_option(pyKinetics$experiment_names),
          selected = state$selected_experiment
        )
      },
      ignoreInit = TRUE
    )

    observeEvent(input$selectedExperiment,
      {
        req(state$traces_loaded)
        req(input$selectedExperiment)

        sel_exp <- input$selectedExperiment

        # IF we have a BLI experiments, leave all the options

        if (sel_exp == "All") {
          types <- unlist(pyKinetics$get_experiment_properties("type"))

          unq_types <- unique(types)

          if (length(unq_types) > 1) {
            pop_up_warning("Processing multiple experiments of different types is not supported. Please select a single experiment.")
            req(FALSE)
          } else {
            py_type <- types[1]
          }
        } else {
          experiment <- pyKinetics$experiments[[sel_exp]]
          py_type <- experiment$type
        }

        if (py_type %in% c("BLI_experiment", "Gator_experiment")) {
          choices <- c(
            "Align association phase" = "align_association",
            "Inter-step correction (dissociation)" = "correct_dissociation"
          )

          if (sel_exp != "All") {
            choices <- c(choices,
              "Subtract baseline"       = "subtract",
              "Average"                 = "average"
            )
          }

          # If we have two experiments, add an option to subtract a whole experiment
          n_exps <- length(pyKinetics$experiment_names)

          if (n_exps > 1 && sel_exp != "All") {
            choices <- c(choices, "Subtract experiment" = "subtract_experiment")
          }

          if (py_type == "BLI_experiment" && sel_exp != "All") {
            choices <- c(choices, "Merge steps" = "merge_steps")
          }
        } else {
          choices <- c("Set time cutoff" = "set_time_cutoff")
        }

        updateSelectInput(session, "operation", choices = choices)
      },
      ignoreInit = TRUE
    )

    # Modal dialog to ask for the input and desired units
    observeEvent(input$triggerProcessing, {
      req(state$traces_loaded)

      sel_exp <- input$selectedExperiment

      if (sel_exp != "All") {
        experiment <- pyKinetics$experiments[[sel_exp]]
      }

      operation <- input$operation

      if (operation == "set_time_cutoff") {
        showModal(modalDialog(
          tags$h3("Please select the time threshold:"),
          numericInput(ns("inputTimeCutOff"), NULL, value = 1, min = 0, max = 1e6),
          tags$script(HTML(paste0("
                    $(document).on('keypress', function(e) {
                        if (e.which == 13) {
                            e.preventDefault();
                            $('#", ns("submitTimeCutOff"), "').click();
                        }
                    });
                    "))),
          footer = tagList(
            actionButton(ns("submitTimeCutOff"), "Submit"),
            modalButton("Cancel")
          )
        ))
      } else {
        # Operations for surface-based experiments require the sensor names
        if (sel_exp != "All") {
          sensor_names <- experiment$sensor_names
        }
      }

      if (operation == "subtract_experiment") {
        # Guess which sensor to subtract using the one with the lowest max signal
        all_experiments <- pyKinetics$experiment_names
        other_experiments <- all_experiments[all_experiments != sel_exp]

        showModal(modalDialog(
          tags$h3("Please select the experiment to use as reference:"),
          selectInput(ns("baselineExperiment"), NULL, other_experiments, selectize = FALSE),
          checkboxInput(ns("expSubtractionIsInPlace"), "Apply in-place subtraction", TRUE),
          tags$script(HTML(paste0("
                    $(document).on('keypress', function(e) {
                        if (e.which == 13) {
                            e.preventDefault();
                            $('#", ns("submitExpSubtraction"), "').click();
                        }
                    });
                    "))),
          footer = tagList(
            actionButton(ns("submitExpSubtraction"), "Submit"),
            modalButton("Cancel")
          )
        ))
      }

      if (operation == "subtract") {
        # Guess which sensor to subtract using the one with the lowest max signal
        ys <- experiment$ys
        idx <- find_probable_baseline(ys)

        sel_sensor_name <- sensor_names[idx]
        sel_sensor_names <- c(sel_sensor_name, sensor_names[sensor_names != sel_sensor_name])

        showModal(modalDialog(
          tags$h3("Please select the baseline:"),
          selectInput(ns("inputBaseline"), NULL, sel_sensor_names, selectize = FALSE),
          tags$script(HTML(paste0("
                    $(document).on('keypress', function(e) {
                        if (e.which == 13) {
                            e.preventDefault();
                            $('#", ns("submitSubtraction1"), "').click();
                        }
                    });
                    "))),
          footer = tagList(
            actionButton(ns("submitSubtraction1"), "Submit"),
            modalButton("Cancel")
          )
        ))
      }

      if (operation %in% c("average", "align_association", "correct_dissociation")) {
        if (sel_exp != "All") {
          rdf <- get_rtable_processing(sensor_names)
          print(rdf)
          output$tableSelection <- renderRHandsontable({
            rdf
          })
        } else {
          dfs <- list()

          for (i in 1:length(pyKinetics$experiment_names)) {
            exp_name <- pyKinetics$experiment_names[i]
            experiment <- pyKinetics$experiments[[exp_name]]
            sensor_names <- experiment$sensor_names

            df_temp <- get_sensor_df(sensor_names, exp_name)
            dfs[[i]] <- df_temp
          }

          df <- do.call(rbind, dfs)

          rdf <- rhandsontable(df) %>%
            hot_col("Select") %>%
            hot_table(stretchH = "all") %>%
            hot_col("ID", readOnly = TRUE) %>%
            hot_col("Experiment", readOnly = TRUE)

          output$tableSelection <- renderRHandsontable({
            rdf
          })
        }
      }


      if (operation == "average") {
        showModal(modalDialog(
          tags$h3("Please select the curves:"),
          tags$script(HTML(paste0("
                    $(document).on('keypress', function(e) {
                        if (e.which == 13) {
                            e.preventDefault();
                            $('#", ns("submitAverage"), "').click();
                        }
                    });
                    "))),
          rHandsontableOutput(ns("tableSelection")),
          textInput(ns("outputNameAverage"), "Output name", "Average"),
          footer = tagList(
            actionButton(ns("submitAverage"), "Submit"),
            modalButton("Cancel")
          )
        ))
      }

      if (operation == "align_association") {
        showModal(modalDialog(
          tags$h3("Please select the curves:"),
          tags$script(HTML(paste0("
                    $(document).on('keypress', function(e) {
                        if (e.which == 13) {
                            e.preventDefault();
                            $('#", ns("submitAlign"), "').click();
                        }
                    });
                    "))),
          rHandsontableOutput(ns("tableSelection")),
          checkboxInput(ns("inPlaceAlignment"), "Perform in-place alignment", TRUE),
          conditionalPanel(
            condition = paste0("input['", ns("inPlaceAlignment"), "']"),
            checkboxInput(ns("createNewSensorNames"), "Create new sensor names", FALSE),
          ),
          checkboxInput(ns("keepRegeneration"), "Keep the regeneration data", FALSE),
          checkboxInput(ns("keepLoading"), "Keep the loading data", FALSE),
          checkboxInput(ns("keepBaseline"), "Keep the baseline data", FALSE),
          checkboxInput(ns("keepActivation"), "Keep the activation data", FALSE),
          checkboxInput(ns("keepQuenching"), "Keep the quenching data", FALSE),
          checkboxInput(ns("keepCustom"), "Keep the custom step data", FALSE),
          tags$script(HTML(paste0("
                    $(document).on('keypress', function(e) {
                        if (e.which == 13) {
                            e.preventDefault();
                            $('#", ns("submitAlign"), "').click();
                        }
                    });
                    "))),
          footer = tagList(
            actionButton(ns("submitAlign"), "Submit"),
            modalButton("Cancel")
          )
        ))
      }

      if (operation == "correct_dissociation") {
        showModal(modalDialog(
          tags$h3("Please select the curves:"),
          tags$script(HTML(paste0("
                    $(document).on('keypress', function(e) {
                        if (e.which == 13) {
                            e.preventDefault();
                            $('#", ns("submitCorrectDis"), "').click();
                        }
                    });
                    "))),
          rHandsontableOutput(ns("tableSelection")),
          tags$h4(""),
          sliderInput(ns("nPointsCorrectDis"), "Number of points to use for correction:", min = 1, max = 20, value = 1),
          checkboxInput(ns("inPlaceCorrection"), "Perform in-place inter-step correction", TRUE),
          conditionalPanel(
            condition = paste0("input['", ns("inPlaceCorrection"), "']"),
            checkboxInput(ns("createNewSensorNames"), "Create new sensor names", FALSE),
          ),
          footer = tagList(
            actionButton(ns("submitCorrectDis"), "Submit"),
            modalButton("Cancel")
          )
        ))
      }

      if (operation == "merge_steps") {
        showModal(modalDialog(
          tags$h3("Please choose the mode for merging steps:"),
          tags$script(HTML(paste0("
                    $(document).on('keypress', function(e) {
                        if (e.which == 13) {
                            e.preventDefault();
                            $('#", ns("submitMerge"), "').click();
                        }
                    });
                    "))),
          selectInput(ns("merging_mode"), NULL, choices = c(
            "By step name"  = "step_name",
            "By step index" = "step_index"
          )),
          footer = tagList(
            actionButton(ns("submitMerge"), "Submit"),
            modalButton("Cancel")
          )
        ))
      }
    })

    observeEvent(input$submitTimeCutOff, {
      removeModal()

      state$traces_loaded <- FALSE

      exp <- pyKinetics$experiments[[input$selectedExperiment]]

      exp$cut_off_time(input$inputTimeCutOff)

      Sys.sleep(0.6)

      # Include the time cutoff step in the logbook
      logbook$append(paste0("Time cutoff set to ", input$inputTimeCutOff, " seconds."), include_time = TRUE, add_empty_line = TRUE)

      state$traces_loaded <- TRUE
    })

    observeEvent(input$submitExpSubtraction, {
      removeModal()
      state$traces_loaded <- FALSE

      exp <- pyKinetics$experiments[[input$selectedExperiment]]

      other_exp <- pyKinetics$experiments[[input$baselineExperiment]]

      result <- tryCatch(
        {
          exp$subtract_experiment(other_exp, inplace = input$expSubtractionIsInPlace)
        },
        error = function(e) {
          if (inherits(e, "python.builtin.RuntimeError")) {
            err <- py_last_error()
            pop_up_warning(
              paste0("⚠ Processing error: ", err$value)
            )
            return("Error")
          } else {
            stop(e) # rethrow non-Python errors
          }
        }
      )

      if (!is.null(result)) {
        req(FALSE)
      }

      logbook$append("Experiment subtraction performed", include_time = TRUE, add_empty_line = TRUE)
      logbook$append(paste0("Experiment ", input$baselineExperiment, " subtracted from ", input$selectedExperiment, "."), include_time = FALSE, add_empty_line = FALSE)

      state$traces_loaded <- TRUE

      Sys.sleep(0.6)

      state$legend_version <- state$legend_version + 1
      state$ligand_info_version <- state$ligand_info_version + 1
    })

    observeEvent(input$submitMerge, {
      removeModal()

      merging_mode <- input$merging_mode

      exp <- pyKinetics$experiments[[input$selectedExperiment]]

      df_steps <- exp$df_steps
      df_steps <- pandas_to_r(df_steps)

      step_names <- unique(df_steps$Name)

      if (merging_mode == "step_name") {
        showModal(modalDialog(
          tags$h3("Please choose the step name."),
          tags$h4('All steps with the specified name will be merged with the step that immediately follows them.
                    To view the outcome of this operation, please check the "Steps" tab.
                    This is useful for correcting mislabeled steps—for example, if a "Baseline" step appears after a "Dissociation" step but
                    actually represents a continuation of the dissociation phase.'),
          tags$script(HTML(paste0("
                    $(document).on('keypress', function(e) {
                        if (e.which == 13) {
                            e.preventDefault();
                            $('#", ns("submitMerge_2"), "').click();
                        }
                    });
                    "))),
          selectInput(ns("step_name"), NULL, choices = step_names),
          footer = tagList(
            actionButton(ns("submitMerge_2"), "Submit"),
            modalButton("Cancel")
          )
        ))
      }

      if (merging_mode == "step_index") {
        showModal(modalDialog(
          tags$h3("Please choose the index of the reference step and the index of the step to combine."),
          tags$h4("As a result, the x and y values of both steps will be combined into a single step.
                    The analyte concentration and column location will be inherited from the reference step.
                    Note: The steps to be merged must have consecutive indexes."),
          tags$script(HTML(paste0("
                    $(document).on('keypress', function(e) {
                        if (e.which == 13) {
                            e.preventDefault();
                            $('#", ns("submitMerge_3"), "').click();
                        }
                    });
                    "))),
          numericInput("step_index_ref", "Reference step index", value = 1, min = 1, max = nrow(df_steps)),
          numericInput("step_index_to_merge", "Step index to merge", value = 2, min = 1, max = nrow(df_steps)),
          footer = tagList(
            actionButton(ns("submitMerge_3"), "Submit"),
            modalButton("Cancel")
          )
        ))
      }
    })

    observeEvent(input$submitMerge_3, {
      removeModal()

      state$traces_loaded <- FALSE

      exp <- pyKinetics$experiments[[input$selectedExperiment]]

      exp$merge_consecutive_steps(
        idx_ref = input$step_index_ref,
        idx_to_merge = input$step_index_to_merge
      )

      state$traces_loaded <- TRUE

      logbook$append("Merging of steps by index performed.", include_time = TRUE, add_empty_line = TRUE)
      # append the selecte experiment, reference index and index to merge to the logbook
      logbook$append(
        paste0(
          "Experiment: ", input$selectedExperiment,
          ", Reference step index: ", input$step_index_ref,
          ", Step index to merge: ", input$step_index_to_merge
        )
      )
    })


    observeEvent(input$submitMerge_2, {
      removeModal()

      state$traces_loaded <- FALSE

      exp <- pyKinetics$experiments[[input$selectedExperiment]]

      exp$merge_consecutive_steps_by_name(step_name = input$step_name)

      state$traces_loaded <- TRUE

      logbook$append("Merging of steps by name performed.", include_time = TRUE, add_empty_line = TRUE)
      logbook$append(
        paste0(
          "Experiment: ", input$selectedExperiment,
          ", Step name: ", input$step_name
        )
      )
    })

    observeEvent(input$submitCorrectDis, {
      removeModal()

      tableCorrectDis <- hot_to_r(input$tableSelection)
      tableCorrectDis <- tableCorrectDis[tableCorrectDis$Select, ]

      state$traces_loaded <- FALSE

      exps_to_analyse <- list()
      sensor_to_correct <- list()

      if (input$selectedExperiment == "All") {
        unq_exps <- unique(tableCorrectDis$Experiment)

        for (exp in unq_exps) {
          temp_df <- tableCorrectDis[tableCorrectDis$Experiment == exp, ]
          samples <- c(temp_df$ID)

          exps_to_analyse[[length(exps_to_analyse) + 1]] <- exp
          sensor_to_correct[[length(sensor_to_correct) + 1]] <- samples
        }
      } else {
        exps_to_analyse[[1]] <- input$selectedExperiment
        sensor_to_correct[[1]] <- c(tableCorrectDis$ID)
      }

      for (exp in exps_to_analyse) {
        samples <- sensor_to_correct[[which(exps_to_analyse == exp)]]

        py_exp <- pyKinetics$experiments[[exp]]

        py_exp$align_dissociation(
          samples,
          input$inPlaceCorrection,
          input$createNewSensorNames,
          input$nPointsCorrectDis
        )
      }

      Sys.sleep(0.6)

      state$legend_version <- state$legend_version + 1
      state$ligand_info_version <- state$ligand_info_version + 1
      # Include the correction step in the logbook

      logbook$append("Inter-step correction performed (between dissociation and association).",
        include_time = TRUE, add_empty_line = TRUE
      )

      logbook$append(paste0("Samples ", paste(samples, collapse = ", "), " corrected."))

      # Append in-place option to the logbook
      logbook$append(paste0("In-place correction:", input$inPlaceCorrection))

      state$traces_loaded <- TRUE
    })

    # Modal dialog to ask for the input and desired units
    observeEvent(input$submitAlign, {
      removeModal()

      tableAlign <- hot_to_r(input$tableSelection)
      tableAlign <- tableAlign[tableAlign$Select, ]

      state$traces_loaded <- FALSE

      sel_exp <- input$selectedExperiment

      exps_to_analyse <- list()
      sensor_to_align <- list()

      if (sel_exp == "All") {
        unq_exps <- unique(tableAlign$Experiment)

        for (exp in unq_exps) {
          temp_df <- tableAlign[tableAlign$Experiment == exp, ]
          samples <- c(temp_df$ID)

          exps_to_analyse[[length(exps_to_analyse) + 1]] <- exp
          sensor_to_align[[length(sensor_to_align) + 1]] <- samples
        }
      } else {
        exps_to_analyse[[1]] <- sel_exp
        sensor_to_align[[1]] <- c(tableAlign$ID)
      }


      for (exp in exps_to_analyse) {
        samples <- sensor_to_align[[which(exps_to_analyse == exp)]]

        py_exp <- pyKinetics$experiments[[exp]]

        py_exp$align_association(
          samples,
          input$inPlaceAlignment,
          input$createNewSensorNames
        )

        # Find new generated names
        all_names <- unlist(py_exp$sensor_names)
        new_names <- all_names[!(all_names %in% samples)]

        # Include the previous names, if we did the alignment it place
        if (input$inPlaceAlignment && !input$createNewSensorNames) {
          new_names <- c(new_names, samples)
        }

        if (!input$keepRegeneration) py_exp$discard_steps(new_names, list("KREGENERATION"))
        if (!input$keepLoading) py_exp$discard_steps(new_names, list("LOADING"))
        if (!input$keepBaseline) py_exp$discard_steps(new_names, list("BASELINE"))
        if (!input$keepActivation) py_exp$discard_steps(new_names, list("ACTIVATION"))
        if (!input$keepQuenching) py_exp$discard_steps(new_names, list("QUENCHING"))
        if (!input$keepCustom) py_exp$discard_steps(new_names, list("CUSTOM"))
      }

      Sys.sleep(0.6)
      state$legend_version <- state$legend_version + 1
      state$ligand_info_version <- state$ligand_info_version + 1

      # Include the alignment step in the logbook
      logbook$append("Alignment step performed.", include_time = TRUE, add_empty_line = TRUE)
      logbook$append(paste0("Samples ", paste(samples, collapse = ", "), " aligned."))

      # Append the discarded steps to the logbook
      if (!input$keepRegeneration) logbook$append("Regeneration steps removed.")
      if (!input$keepLoading) logbook$append("Loading steps removed.")
      if (!input$keepBaseline) logbook$append("Baseline steps removed.")
      if (!input$keepActivation) logbook$append("Activation steps removed.")
      if (!input$keepQuenching) logbook$append("Quenching steps removed.")
      if (!input$keepCustom) logbook$append("Custom steps removed.")

      # Append in-place option to the logbook
      logbook$append(paste0("In-place alignment:", input$inPlaceAlignment))

      state$traces_loaded <- TRUE
    })


    # Modal dialog to ask for the input and desired units
    observeEvent(input$submitAverage, {
      removeModal()

      tableAverage <- hot_to_r(input$tableSelection)
      samples <- c(tableAverage$ID[tableAverage$Select])

      if (length(samples) < 2) {
        req(FALSE)
      }

      state$traces_loaded <- FALSE

      exp <- pyKinetics$experiments[[input$selectedExperiment]]

      exp$average(samples, input$outputNameAverage)

      labels <- unlist(pyKinetics$get_experiment_properties("sensor_names"))
      ids <- unlist(pyKinetics$get_experiment_properties("sensor_names_unique"))

      legend_df(get_plotting_df(ids, labels))

      Sys.sleep(0.6)
      state$legend_version <- state$legend_version + 1
      state$ligand_info_version <- state$ligand_info_version + 1

      # Include the average step in the logbook
      logbook$append("Average step performed.", include_time = TRUE, add_empty_line = TRUE)
      logbook$append(paste0("Samples ", paste(samples, collapse = ", "), " averaged."))

      state$traces_loaded <- TRUE
    })

    # Modal dialog to ask for the input and desired units
    observeEvent(input$submitSubtraction1, {
      exp <- pyKinetics$experiments[[input$selectedExperiment]]
      names <- exp$sensor_names
      names <- names[names != input$inputBaseline]

      removeModal()

      showModal(modalDialog(
        tags$h3("Please select the sample(s):"),
        tags$script(HTML(paste0("
                $(document).on('keypress', function(e) {
                    if (e.which == 13) {
                        e.preventDefault();
                        $('#", ns("submitSubtraction2"), "').click();
                    }
                });
                "))),
        rHandsontableOutput(ns("tableSubtraction")),
        checkboxInput(ns("inPlaceSubtraction"), "Perform in-place subtraction", TRUE),
        footer = tagList(
          actionButton(ns("submitSubtraction2"), "Submit"),
          modalButton("Cancel")
        )
      ))
    })

    observeEvent(input$submitSubtraction2, {
      removeModal()

      state$traces_loaded <- FALSE

      tableSubtraction <- hot_to_r(input$tableSubtraction)
      samples <- c(tableSubtraction$ID[tableSubtraction$Select])

      exp <- pyKinetics$experiments[[input$selectedExperiment]]

      exp$subtraction(as.list(samples), input$inputBaseline, input$inPlaceSubtraction)

      labels <- unlist(pyKinetics$get_experiment_properties("sensor_names"))
      ids <- unlist(pyKinetics$get_experiment_properties("sensor_names_unique"))

      legend_df(get_plotting_df(ids, labels))

      Sys.sleep(0.6)

      state$legend_version <- state$legend_version + 1
      state$ligand_info_version <- state$ligand_info_version + 1

      # Include the subtraction step in the logbook
      logbook$append("Subtraction step performed.", include_time = TRUE, add_empty_line = TRUE)
      logbook$append(paste0("Samples ", paste(samples, collapse = ", "), " subtracted from ", input$inputBaseline))

      # In-place subtraction
      logbook$append(paste0("In-place subtraction:", input$inPlaceSubtraction))

      state$traces_loaded <- TRUE
    })

    output$tableSubtraction <- renderRHandsontable({
      req(state$traces_loaded)
      req(input$inputBaseline)

      names <- pyKinetics$experiments[[input$selectedExperiment]]$sensor_names

      rdf <- get_rtable_processing(names[names != input$inputBaseline])

      rdf
    })
  })
}
