box::use(
  .. / .. / busy_indicator[
    withBusyIndicatorUI
  ],
  .. / .. / dialogs[
    pop_up_warning
  ],
  reticulate[
    py_last_error
  ],
  shiny[
    actionButton,
    column,
    fileInput,
    fluidRow,
    HTML,
    icon,
    modalButton,
    modalDialog,
    moduleServer,
    NS,
    observeEvent,
    p,
    req,
    selectInput,
    showModal,
    sliderInput,
    span,
    tagList,
    textInput,
    updateSelectInput,
    tags,
    removeModal
  ],
  shinydashboard[
    box
  ],
  shinyjs[
    hidden
  ],
  tippy[
    tippy_this
  ]
)

#' @export
loadDataUI <- function(id) {
  ns <- NS(id)

  tagList(
    box(
      title = "1. Input", width = 12, solidHeader = T, status = "primary",
      fluidRow(
        column(9, p(
          HTML("<b>Kinetic data files</b>"),
          span(icon("info-circle"), id = "info_uu1-1"),
          fileInput(ns("kineticFiles"), NULL, multiple = TRUE),
          tippy_this(
            elementId = "info_uu1-1",
            tooltip = "For Octet BLI, load (simultaneously) all the sensor files (.frd)  and the Method file (.fmf).",
            placement = "right"
          )
        )),
        column(3, p(
          HTML("<b></b>"),
          withBusyIndicatorUI(hidden(actionButton(ns("Go"), "", class = "btn-primary")))
        ))
      ),
      fluidRow(
        column(4, p(
          HTML("<b></b>"),
          actionButton(
            inputId = ns("loadExampleData"),
            label = "Load example data!",
            icon("caret-right"),
            style = "color: #FFFFFF; background-color: #00829c;
                    border-color: #00829c"
          )
        ))
      ),
      fluidRow(
        column(8, p(
          HTML("<b>Select experiment to delete</b>"),
          selectInput(ns("experiment2delete"), NULL, c("None" = "None"), selectize = FALSE)
        )),
        column(
          4, p(HTML("<b><br></b>")),
          actionButton(
            inputId = ns("triggerDeletion"), label = "",
            icon("trash-can"),
            style = "color: #0E090D; background-color: #DABFDF;
                    border-color: #6A4D71"
          )
        )
      )
    )
  )
}

#' @export
loadDataServer <- function(id, state, pyKinetics, pykingenie, logbook) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$loadExampleData, {
      state$traces_loaded <- FALSE
      state$surface_based_binding <- TRUE
      # Clear all experiments
      pyKinetics$delete_experiment(pyKinetics$experiment_names)

      # list files in test_bli_folder
      files <- list.files("./www/test_bli_folder", full.names = TRUE)
      names <- list.files("./www/test_bli_folder")

      newExperimentName <- "Example"

      bli_experiment <- pykingenie$OctetExperiment(newExperimentName)

      bli_experiment$read_sample_plate_info(files, names)
      bli_experiment$read_sensor_data(files, names)

      pyKinetics$add_experiment(bli_experiment, newExperimentName)

      state$legend_version <- state$legend_version + 1
      state$selected_experiment <- newExperimentName

      updateSelectInput(session, "experiment2delete", choices = c("ALL", pyKinetics$experiment_names))
      state$ligand_info_version <- state$ligand_info_version + 1

      state$traces_loaded <- TRUE
      state$show_ligand_info <- TRUE
    })

    observeEvent(input$kineticFiles, {
      req(input$kineticFiles)

      files <- input$kineticFiles$datapath

      frd_files <- files[grepl(".frd", files)]
      csv_files <- files[grepl(".csv", files)]

      n_frd <- length(frd_files)
      n_csv <- length(csv_files)

      if (n_frd == 0 && n_csv == 0) {
        pop_up_warning("No .frd or .csv files found.")
        req(FALSE)
      }

      if (n_frd > 0) {
        newExperimentName <- pykingenie$guess_experiment_name(frd_files[1])
      } else {
        newExperimentName <- "Experiment"
      }

      showModal(modalDialog(
        tags$h3("Set the experiment name:"),
        textInput(ns("newExperimentName"), NULL, newExperimentName),
        footer = tagList(
          actionButton(ns("submitLoadExperiment"), "Submit"),
          modalButton("Cancel")
        )
      ))
    })

    observeEvent(input$submitLoadExperiment, {
      removeModal()
      traces_loaded <- FALSE

      files <- input$kineticFiles$datapath
      names <- input$kineticFiles$name

      sorted_indices <- order(names)

      files <- files[sorted_indices]
      names <- names[sorted_indices]

      # Check if the type of experiment
      exp_type <- pykingenie$guess_experiment_type(files)

      # Check if it is surface-based
      exp_is_surface_based <- exp_type == "surface"

      # Check if we have data loaded and if it is of the same type
      if (state$traces_loaded) {
        if (exp_is_surface_based && !state$surface_based_binding) {
          pop_up_warning("Oops! We can not load the file,
                    because the type of experiment (surface-based) is different from the one loaded.")
          req(FALSE)
        }

        if (!exp_is_surface_based && state$surface_based_binding) {
          pop_up_warning("Oops! We can not load the file,
                    because the type of experiment (solution-based) is different from the one loaded.")
          req(FALSE)
        }
      }

      state$surface_based_binding <- exp_is_surface_based

      state$traces_loaded <- FALSE

      frd_files <- files[grepl(".frd", files)]
      n_frds <- length(frd_files)

      if (n_frds > 0) {
        bli_experiment <- pykingenie$OctetExperiment(input$newExperimentName)
        bli_experiment$read_sample_plate_info(files, names)

        result <- tryCatch(
          {
            bli_experiment$read_sensor_data(files, names)
          },
          error = function(e) {
            if (inherits(e, "python.builtin.ValueError")) {
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

        pyKinetics$add_experiment(bli_experiment, input$newExperimentName)

        if (bli_experiment$traces_loaded) {
          state$legend_version <- state$legend_version + 1
          state$selected_experiment <- input$newExperimentName
        }

        traces_loaded <- any(unlist(pyKinetics$get_experiment_properties("traces_loaded")))
      } else if (any(grepl("simulation_KinGenie_", names))) {
        nFiles <- length(files)

        for (i in 1:nFiles) {
          if (nFiles > 1) {
            newExperimentName <- paste0(input$newExperimentName, i)
          } else {
            newExperimentName <- input$newExperimentName
          }

          if (exp_is_surface_based) {
            kingenie_simulation <- pykingenie$KinGenieCsv(newExperimentName)
          } else {
            kingenie_simulation <- pykingenie$KinGenieCsvSolution(newExperimentName)
          }

          kingenie_simulation$read_csv(files[i])
          pyKinetics$add_experiment(kingenie_simulation, newExperimentName)
        }

        state$legend_version <- state$legend_version + 1
        state$selected_experiment <- input$newExperimentName

        traces_loaded <- TRUE
      } else if (any(grepl(".ini", files))) {
        gt <- pykingenie$GatorExperiment(input$newExperimentName)
        gt$read_all_gator_data(files, names)
        pyKinetics$add_experiment(gt, input$newExperimentName)

        Sys.sleep(0.5)
        state$legend_version <- state$legend_version + 1
        state$selected_experiment <- input$newExperimentName

        traces_loaded <- any(unlist(pyKinetics$get_experiment_properties("traces_loaded")))
      } else if (all(grepl(".csv", files))) {
        showModal(modalDialog(
          tags$h3("Please select if the data corresponds to surface or in-solution binding:"),
          selectInput(ns("bindingType"), NULL, choices = c(
            "Surface-based binding"  = "surface",
            "In-solution binding"    = "solution"
          )),
          tags$script(HTML("
                    $(document).on('keypress', function(e) {
                        if(e.which == 13) {
                        e.preventDefault();
                        $('#" + ns("importCSVs") + "').click();
                        }
                    });
                    ")),
          footer = tagList(
            actionButton(ns("importCSVs"), "Submit"),
            modalButton("Cancel")
          )
        ))

        req(FALSE)
      }

      Sys.sleep(0.6)

      state$traces_loaded <- traces_loaded
      state$ligand_info_version <- state$ligand_info_version + 1
      state$show_ligand_info <- TRUE

      if (!traces_loaded) {
        pop_up_warning("No traces were loaded. Please check the input files.")
        req(FALSE)
      }

      logbook$append(paste0("Loading experiment: ", input$newExperimentName), include_time = TRUE, add_empty_line = TRUE)

      # Append files to logbook
      logbook$append(paste0("Files: ", paste(names, collapse = ", ")), include_time = FALSE)

      updateSelectInput(session, "experiment2delete", choices = c("ALL", pyKinetics$experiment_names))
    })

    observeEvent(input$importCSVs, {
      removeModal()

      exp_is_surface_based <- input$bindingType == "surface"

      state$surface_based_binding <- exp_is_surface_based

      files <- input$kineticFiles$datapath
      names <- input$kineticFiles$name

      sorted_indices <- order(names)

      files <- files[sorted_indices]
      names <- names[sorted_indices]

      nFiles <- length(files)

      for (i in 1:nFiles) {
        if (nFiles > 1) {
          newExperimentName <- paste0(input$newExperimentName, i)
        } else {
          newExperimentName <- input$newExperimentName
        }

        if (exp_is_surface_based) {
          kingenie_csv <- pykingenie$KinGenieCsv(newExperimentName)
        } else {
          kingenie_csv <- pykingenie$KinGenieCsvSolution(newExperimentName)
        }

        kingenie_csv$read_csv(files[i])
        pyKinetics$add_experiment(kingenie_csv, newExperimentName)
      }

      pyKinetics$collapse_solution_experiments(input$newExperimentName)

      state$legend_version <- state$legend_version + 1
      state$selected_experiment <- input$newExperimentName

      traces_loaded <- TRUE
      state$show_ligand_info <- TRUE
      Sys.sleep(0.6)

      state$traces_loaded <- traces_loaded
      state$ligand_info_version <- state$ligand_info_version + 1

      if (!traces_loaded) {
        pop_up_warning("No traces were loaded. Please check the input files.")
        req(FALSE)
      }

      logbook$append(paste0("Loading experiment: ", input$newExperimentName), include_time = TRUE, add_empty_line = TRUE)

      # Append files to logbook
      logbook$append(paste0("Files: ", paste(names, collapse = ", ")), include_time = FALSE)

      updateSelectInput(session, "experiment2delete", choices = c("ALL", pyKinetics$experiment_names))
    })

    observeEvent(input$triggerDeletion, {
      req(state$traces_loaded)

      experiment2delete <- input$experiment2delete

      if (experiment2delete == "None") {
        req(FALSE)
      }

      state$traces_loaded <- FALSE

      if (experiment2delete == "ALL") {
        pyKinetics$delete_experiment(pyKinetics$experiment_names)

        updateSelectInput(session, "experiment2delete", NULL, "None")
        state$selected_experiment <- NULL

        logbook$append("All experiments were deleted.", include_time = TRUE, add_empty_line = TRUE)
      } else {
        pyKinetics$delete_experiment(experiment2delete)

        logbook$append(paste0("Experiment ", experiment2delete, " was deleted."), include_time = TRUE, add_empty_line = TRUE)

        current_experiments <- pyKinetics$experiment_names

        if (length(current_experiments) == 0) {
          updateSelectInput(session, "experiment2delete", NULL, "None")
          state$selected_experiment <- NULL
        } else {
          updateSelectInput(session, "experiment2delete", NULL, c("ALL", pyKinetics$experiment_names))

          state$legend_version <- state$legend_version + 1
          state$ligand_info_version <- state$ligand_info_version + 1

          state$traces_loaded <- TRUE
          state$show_ligand_info <- TRUE
        }
      }
    })
  })
}
