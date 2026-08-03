box::use(
  .. / .. / busy_indicator[
    withBusyIndicatorServer,
    withBusyIndicatorUI
  ],
  .. / .. / dialogs[
    pop_up_info,
    pop_up_success,
    pop_up_warning
  ],
  .. / .. / helpers[
    df_to_lines
  ],
  reticulate[
    py_last_error
  ],
  shiny[
    actionButton,
    checkboxInput,
    column,
    conditionalPanel,
    fluidRow,
    HTML,
    icon,
    modalDialog,
    moduleServer,
    NS,
    observeEvent,
    p,
    removeModal,
    renderUI,
    req,
    selectInput,
    showModal,
    span,
    tagList,
    tags,
    uiOutput,
    updateSelectInput
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
fitControlsUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        12,
        uiOutput(ns("fitting_box"))
      )
    )
  )
}

#' @export
fitControlsServer <- function(id, state, dataset, pyKinetics, logbook) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$fitting_box <- renderUI({
      req(state$traces_loaded)

      if (state$surface_based_binding) {
        box(
          title = "2. Fitting", width = 12, solidHeader = T, status = "primary",
          fluidRow(
            column(3, p(
              HTML('<p style="margin-bottom:0px;"><br></p>'),
              actionButton(
                inputId = ns("triggerCreateDataset"), label = "a. Create dataset",
                style = "color: #fff; background-color: #337ab7;
                            border-color: #2e6da4"
              )
            )),
            column(4, p(
              HTML("<b>b. Region</b>"),
              span(icon("info-circle"), id = "info_uu-fitRegion"),
              selectInput(ns("fittingRegion"), NULL,
                c(
                  "Association and dissociation" = "association_dissociation",
                  "Steady-state" = "steady_state",
                  "Association" = "association",
                  "Dissociation" = "dissociation"
                ),
                selectize = FALSE
              ),
              tippy_this(
                elementId = "info_uu-fitRegion",
                tooltip = "Select association and/or dissociation to fit K_d and k_off.
                                Select Steady-state to fit K_d.'.
                                ", placement = "right"
              )
            )),
            column(4, p(
              HTML("<b>c. Model</b>"),
              span(icon("info-circle"), id = "info_uu-fitModel"),
              selectInput(
                ns("fittingModel"), NULL,
                c(
                  "One-to-one" = "one_to_one",
                  "One-to-one (MTL)" = "one_to_one_mtl"
                )
              ),
              tippy_this(
                elementId = "info_uu-fitModel",
                tooltip = "Select and press 'Fit!'.
                                ", placement = "right"
              )
            ))
          ),
          fluidRow(
            conditionalPanel(
              "input.fittingModel == 'two_to_one'",
              ns = ns,
              column(
                2,
                p(
                  HTML("<b>Fit sigma</b>"),
                  checkboxInput(ns("fitSigmaTwoToOne"), NULL, FALSE)
                )
              )
            ),
            column(2, p(
              HTML('<p style="margin-bottom:0px;"><br></p>'),
              actionButton(
                inputId = ns("triggerFitting"), label = "d. Fit!",
                icon("meteor"),
                style = "color: #fff; background-color: #337ab7;
                            border-color: #2e6da4"
              )
            )),

            # Little hack to use the withBusyIndicatorUI function (loading spinner)
            column(1, p(
              HTML("<b><br></b>"),
              withBusyIndicatorUI(
                hidden(actionButton(ns("fitSurfaceBusyIndicator"), "", class = "btn-primary"))
              )
            )),

            # Button to calculate the asymmetric error
            column(2, p(
              HTML('<p style="margin-bottom:0px;"><br></p>'),
              actionButton(
                inputId = ns("triggerAsymError"), label = "e. Confidence interval",
                icon("calculator"),
                style = "color: #fff; background-color: #337ab7;
                            border-color: #2e6da4"
              )
            ))
          )
        )
      } else {
        box(
          title = "2. Fitting", width = 12, solidHeader = T, status = "primary",
          fluidRow(
            column(3, p(
              HTML('<p style="margin-bottom:0px;"><br></p>'),
              actionButton(
                inputId = ns("triggerCreateDatasetSolution"), label = "a. Create dataset",
                style = "color: #fff; background-color: #337ab7;
                            border-color: #2e6da4"
              )
            )),

            # Model selection - between single and double exponential
            column(3, p(
              HTML("b. Model"),
              selectInput(
                inputId = ns("modelSelectionSolution"),
                label = NULL,
                choices = c(
                  "Single exponential" = "single",
                  "Double exponential" = "double",
                  "One binding site" = "one_binding_site",
                  "One binding site (induced fit)" = "one_binding_site_if"
                ),
                selected = "single",
                selectize = FALSE,
              )
            )),
            column(3, p(
              HTML('<p style="margin-bottom:0px;"><br></p>'),
              actionButton(
                inputId = ns("triggerFitSolution"), label = "c. Fit",
                style = "color: #fff; background-color: #337ab7;
                            border-color: #2e6da4"
              )
            ))
          )
        )
      }
    })

    observeEvent(input$fittingRegion, {
      if (input$fittingRegion == "association_dissociation") {
        updateSelectInput(session, "fittingModel",
          choices = c(
            "One-to-one" = "one_to_one",
            "One-to-one (MTL)" = "one_to_one_mtl",
            "One-to-one (induced fit)" = "one_to_one_if",
            "Two-to-one" = "two_to_one"
          )
        )
      } else {
        updateSelectInput(session, "fittingModel",
          choices = c(
            "One-to-one" = "one_to_one",
            "Two-to-one" = "two_to_one"
          )
        )
      }
    })

    observeEvent(input$triggerCreateDataset, {
      req(state$traces_loaded)
      req(dataset())

      df <- dataset()

      # First filter, verify that we do not have any negative concentration in rows where the 'Select' column is TRUE
      sel_concs <- df[, 2][df$Select] # The second column is the ligand concentration, both for surface-based and solution-based binding

      if (any(sel_concs < 0)) {
        pop_up_warning("Negative ligand concentrations found in selected rows.
                Please correct them before creating the dataset.")
        req(FALSE)
      }

      for (name in c(
        "kinetics_table_shown",
        "ss_plot_shown",
        "ss_table_shown",
        "kinetics_ci_95_table_shown",
        "diagnostics_plot_shown",
        "fit_dataset_loaded",
        "ss_fit_done",
        "kinetics_fit_done"
      )) {
        state[[name]] <- FALSE
      }

      pyKinetics$init_fittings()

      logbook$append("Dataset creation triggered", include_time = TRUE, add_empty_line = TRUE)
      logbook$append(df_to_lines(df))

      logbook_messages <- pyKinetics$generate_fittings(df)

      for (msg in logbook_messages) {
        logbook$append(msg, add_empty_line = TRUE)
      }

      state$is_single_cycle <- any(pyKinetics$get_experiment_properties("is_single_cycle", fittings = TRUE))

      if (state$is_single_cycle) {
        # remove  the option to fit the steady-state region
        updateSelectInput(session, "fittingRegion",
          choices = c(
            "Association and dissociation" = "association_dissociation",
            "Dissociation" = "dissociation"
          )
        )
      }

      state$fit_dataset_loaded <- TRUE
      state$show_ligand_info <- FALSE

      pop_up_info(
        "Dataset(s) created.
                Traces with the same sample ID will share thermodynamic parameters
                (e.g., <i>K</i><sub>d</sub> and <i>k</i><sub>off</sub>).
                Traces with the same sample ID and Smax ID will share the <i>Smax</i> parameter,
                if the option linked by Smax is set to TRUE."
      )
    })

    # Show or hide the confidence interval button
    observeEvent(list(input$fittingRegion, input$fittingModel), {
      req(state$fit_dataset_loaded)

      if (input$fittingRegion == "steady_state") {
        state$ss_plot_shown <- TRUE
      }

      # Show the button if the fitting region is association and dissociation
      # and the selected model is one to one
      if (
        input$fittingRegion == "association_dissociation" &&
          input$fittingModel == "one_to_one"
      ) {
        shinyjs::show("triggerAsymError")
      } else {
        shinyjs::hide("triggerAsymError")
      }
    })

    observeEvent(input$triggerFitting, {
      req(state$fit_dataset_loaded)

      state$kinetics_fit_done <- FALSE
      state$ss_fit_done <- FALSE

      state$fit_dataset_loaded <- FALSE

      result <- tryCatch(
        {
          pyKinetics$submit_steady_state_fitting()
        },
        error = function(e) {
          if (inherits(e, "python.builtin.RuntimeError")) {
            err <- py_last_error()
            pop_up_warning(
              paste0("⚠ Fitting error: ", err$value)
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

      fittingRegion <- input$fittingRegion

      if (fittingRegion == "steady_state") {
        fittingModel <- input$fittingModel

        # Re fit if the user-selected the two_to_one model, because the steady-state fitting is done with the one-to-one model
        if (fittingModel == "two_to_one") {
          pyKinetics$submit_steady_state_fitting(
            fittingModel,
            fit_sigma = input$fitSigmaTwoToOne
          )
        }


        pop_up_success(
          "Steady state fitting done.
                    Remember that the model assumes
                    equal binding capacity for all the samples with the same combination of
                    sample ID,
                    and Smax ID."
        )

        Sys.sleep(0.5)

        logbook$append("Steady state fitting done.", add_empty_line = TRUE)

        state$ss_fit_done <- TRUE
        state$fit_dataset_loaded <- TRUE

        state$kinetics_fit_done <- FALSE

        state$ss_plot_shown <- TRUE
        state$ss_table_shown <- TRUE
      } else {
        showModal(modalDialog(
          tags$h4(
            "Can you assume equal binding capacity for all the samples with the same combination
                    of sample ID, experiment and Smax ID?
                    If yes, set the option below to TRUE.
                    Hint: If the loading levels are different, you should set this option to FALSE.
                    Only if the same sensor was used for all the ligand concentrations in a series,
                    you should set this option to TRUE."
          ),
          checkboxInput(ns("linkedRmax"), "Smax Linked By Sensor", FALSE),
          footer = tagList(
            actionButton(ns("submitKineticsFitting"), "Submit"),
            actionButton(ns("cancelKineticsFitting"), "Cancel")
          )
        ))
      }
    })

    observeEvent(input$submitKineticsFitting, {
      removeModal()

      withBusyIndicatorServer(ns("fitSurfaceBusyIndicator"), {
        state$kinetics_fit_done <- FALSE

        fittingModel <- input$fittingModel
        fittingRegion <- input$fittingRegion
        shared_smax <- input$linkedRmax

        result <- tryCatch(
          {
            pyKinetics$submit_kinetics_fitting(
              fitting_model = fittingModel,
              fitting_region = fittingRegion,
              shared_smax = shared_smax,
              fit_sigma = input$fitSigmaTwoToOne
            )
          },
          error = function(e) {
            if (inherits(e, "python.builtin.RuntimeError")) {
              err <- py_last_error()
              pop_up_warning(
                paste0("⚠ Fitting error: ", err$value)
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

        pop_up_success(
          "Fitting done.
                    The fitted parameters are shown in the 'Fitted params (Kinetics)' table.
                    For the fitting algorithm to work, we use boundaries.
                    Check them in the 'Fitted params (boundaries)' table. If a fitted parameter is too close to the
                    lower or upper boundary, it will be highlighted in red."
        )

        logbook$append("Kinetics fitting done.", add_empty_line = TRUE)
        logbook$append(paste0("Region:", fittingRegion), add_empty_line = FALSE)
        logbook$append(paste0("Model:", fittingModel), add_empty_line = FALSE)

        logbook$append(paste0("Linked Rmax:", input$linkedRmax))

        state$kinetics_fit_done <- TRUE
        state$fit_dataset_loaded <- TRUE

        state$kinetics_table_shown <- TRUE

        c1 <- fittingModel == "one_to_one"
        c2 <- grepl("asso", fittingRegion)

        state$diagnostics_plot_shown <- c1 & c2
      })
    })

    observeEvent(input$cancelKineticsFitting, {
      removeModal()
      for (name in c(
        "kinetics_table_shown",
        "ss_plot_shown",
        "ss_table_shown",
        "kinetics_ci_95_table_shown",
        "diagnostics_plot_shown",
        "ss_fit_done",
        "kinetics_fit_done"
      )) {
        state[[name]] <- FALSE
      }
      pop_up_info("Kinetics fitting was cancelled.")
    })

    observeEvent(input$triggerAsymError, {
      req(state$fit_dataset_loaded)
      req(state$kinetics_fit_done)

      pop_up_info(
        "Calculating the asymmetric error for K_d and k_off.
                Please wait some minutes for the calculations to finish."
      )

      pyKinetics$calculate_asymmetric_error()

      # Create the tab with the asymmetric error

      state$kinetics_ci_95_table_shown <- TRUE

      pop_up_success("Asymmetric error calculated for K_d and k_off.")

      logbook$append("Asymmetric error calculated for K_d and k_off.", add_empty_line = TRUE)
    })

    observeEvent(input$triggerCreateDatasetSolution, {
      req(state$traces_loaded)
      req(dataset())
      req(!state$surface_based_binding)


      state$kinetics_fit_done <- FALSE
      state$fit_dataset_loaded <- FALSE

      # Set the checkboxInput showLigandInfo to FALSE
      state$show_ligand_info <- FALSE

      pyKinetics$init_fittings()
      df <- dataset()


      logbook$append("Dataset creation triggered", include_time = TRUE, add_empty_line = TRUE)
      logbook$append(df_to_lines(df))

      logbook_messages <- pyKinetics$generate_fittings_solution(df)

      for (msg in logbook_messages) {
        logbook$append(msg, add_empty_line = TRUE)
      }

      state$fit_dataset_loaded <- TRUE
    })

    observeEvent(input$triggerFitSolution, {
      req(state$fit_dataset_loaded)
      req(input$modelSelectionSolution)

      state$kinetics_fit_done <- FALSE

      model <- input$modelSelectionSolution

      logbook$append("Fitting triggered", include_time = TRUE, add_empty_line = TRUE)
      # append model selection to the logbook
      logbook$append(paste0("Model selection: ", model), add_empty_line = FALSE)

      # If the model is induced fit - we need to ask which signals to fit
      if (model == "one_binding_site_if") {
        # Create a modal dialog to ask the user if they want to fit the signal of
        # the free protein, the free ligand, the intermediate complex, or the bound complex
        # and to select if the signal of the intermediate complex equals the signal of the trapped complex
        # Each of them is a checkboxInput

        showModal(modalDialog(
          tags$h3("Please select the species producing the signal:"),
          tags$script(HTML(paste0("
                    $(document).on('keypress', function(e) {
                        if(e.which == 13) {
                        e.preventDefault();
                        $('#", ns("submitIF"), "').click();
                        }
                    });
                    "))),
          checkboxInput(ns("fit_signal_E"), "Free protein", FALSE),
          checkboxInput(ns("fit_signal_S"), "Free ligand", FALSE),
          checkboxInput(ns("fit_signal_ES"), "Complex", TRUE),
          tags$h4("Please select if the signal produced by the induced complex equals the signal produced by the intermediate complex:"),
          checkboxInput(ns("ESint_equals_ES"), "Signal(ES) = Signal(ES_intermediate)", TRUE),
          footer = tagList(
            actionButton(ns("submitIF"), "Submit"),
            modalButton("Cancel")
          )
        ))

        req(FALSE)
      }

      # If the model is simple one binding site - we need to ask which signals to fit and if we need to fit t0
      if (model == "one_binding_site") {
        # Create a modal dialog to ask the user if they want to fit the signal of
        # the free protein, the free ligand, or the bound complex
        # and to select if they want to fit t0
        # Each of them is a checkboxInput

        showModal(modalDialog(
          tags$h3("Please select the species producing the signal:"),
          tags$script(HTML(paste0("
                    $(document).on('keypress', function(e) {
                        if(e.which == 13) {
                        e.preventDefault();
                        $('#", ns("submitOBS"), "').click();
                        }
                    });
                    "))),
          checkboxInput(ns("fit_signal_E"), "Free protein", FALSE),
          checkboxInput(ns("fit_signal_S"), "Free ligand", FALSE),
          checkboxInput(ns("fit_signal_ES"), "Complex", TRUE),
          checkboxInput(ns("fit_t0"), "Fit t0", TRUE),
          footer = tagList(
            actionButton(ns("submitOBS"), "Submit"),
            modalButton("Cancel")
          )
        ))

        req(FALSE)
      }

      pyKinetics$submit_fitting_solution(model)

      pop_up_success("The fitting has been completed.")

      state$solution_model <- model

      state$k_obs_plot_shown <- model == "double"

      state$kinetics_table_shown_sol <- TRUE
      state$kinetics_fit_done <- TRUE
    })


    observeEvent(input$submitIF, {
      # close the modal dialog
      removeModal()

      pop_up_info("Fitting with induced fit model is in progress.
            A new popup will appear when the fitting is done. Please wait...")

      pyKinetics$submit_fitting_solution(
        "one_binding_site_if",
        fit_signal_E = input$fit_signal_E,
        fit_signal_S = input$fit_signal_S,
        fit_signal_ES = input$fit_signal_ES,
        ESint_equals_ES = input$ESint_equals_ES
      )

      pop_up_success("The fitting has been completed.")

      state$kinetics_table_shown_sol <- TRUE
      state$kinetics_fit_done <- TRUE
    })

    observeEvent(input$submitOBS, {
      # close the modal dialog
      removeModal()

      pop_up_info("Fitting with one binding site model is in progress.
            A new popup will appear when the fitting is done. Please wait...")

      pyKinetics$submit_fitting_solution(
        "one_binding_site",
        fit_signal_E = input$fit_signal_E,
        fit_signal_S = input$fit_signal_S,
        fit_signal_ES = input$fit_signal_ES,
        fixed_t0 = !input$fit_t0
      )

      pop_up_success("The fitting has been completed.")

      state$kinetics_table_shown_sol <- TRUE
      state$kinetics_fit_done <- TRUE
    })
  })
}
