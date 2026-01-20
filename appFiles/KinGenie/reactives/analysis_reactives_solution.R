removeDynamicTabsSolution <- function() {

    # Remove all the tabs that we dynamically create
    if (reactives$k_obs_plot_shown) {

        removeTab(
          inputId='tabBoxFitInSolution',
          target = 'k_obs',
          session = session
        )

        reactives$k_obs_plot_shown <- FALSE

    }

    if (reactives$kinetics_table_shown_sol) {

        removeTab(
          inputId='tabBoxFitInSolution',
          target = 'Fitted params (Kinetics)',
          session = session
        )

        reactives$kinetics_table_shown_sol <- FALSE

    }

}

observeEvent(input$triggerCreateDatasetSolution, {

    req(reactives$traces_loaded)
    req(input$legendInfo)
    req(!reactives$surface_based_binding)

    removeDynamicTabsSolution()

    reactives$kinetics_fit_done  <- FALSE
    reactives$fit_dataset_loaded <- FALSE

    # Set the checkboxInput showLigandInfo to FALSE
    if (input$showLigandInfo) {
        updateCheckboxInput(session, "showLigandInfo", value = FALSE)
        Sys.sleep(1)
    }

    pyKinetics$init_fittings()
    df <- hot_to_r(input$ligandInfo)

    append_record_to_logbook('Dataset creation triggered',include_time = TRUE,add_empty_line = TRUE)
    append_record_to_logbook(df_to_lines(df))

    logbook_messages <- pyKinetics$generate_fittings_solution(df)

    for (msg in logbook_messages) {

        append_record_to_logbook(msg,add_empty_line = TRUE)

    }

    reactives$fit_dataset_loaded <- TRUE

})

observeEvent(input$triggerFitSolution, {

    req(reactives$fit_dataset_loaded)
    req(input$modelSelectionSolution)

    removeDynamicTabsSolution()

    reactives$kinetics_fit_done <- FALSE

    model <- input$modelSelectionSolution

    append_record_to_logbook('Fitting triggered',include_time = TRUE,add_empty_line = TRUE)
    # append model selection to the logbook
    append_record_to_logbook(paste0('Model selection: ', model), add_empty_line = FALSE)

    # If the model is induced fit - we need to ask which signals to fit
    if (model == "one_binding_site_if") {

        # Create a modal dialog to ask the user if they want to fit the signal of
        # the free protein, the free ligand, the intermediate complex, or the bound complex
        # and to select if the signal of the intermediate complex equals the signal of the trapped complex
        # Each of them is a checkboxInput

        showModal(modalDialog(

            tags$h3('Please select the species producing the signal:'),

            tags$script(HTML("
              $(document).on('keypress', function(e) {
                if(e.which == 13) {
                  e.preventDefault();
                  $('#submitIF').click();
                }
              });
            ")),

            checkboxInput('fit_signal_E','Free protein',FALSE),
            checkboxInput('fit_signal_S','Free ligand',FALSE),
            checkboxInput('fit_signal_ES','Complex',TRUE),

            tags$h4('Please select if the signal produced by the induced complex equals the signal produced by the intermediate complex:'),

            checkboxInput('ESint_equals_ES','Signal(ES) = Signal(ES_intermediate)',TRUE),

            footer=tagList(
              actionButton('submitIF', 'Submit'),
              modalButton('Cancel')
            )
        ))

        return(NULL)
    }

    # If the model is simple one binding site - we need to ask which signals to fit and if we need to fit t0
    if (model == "one_binding_site") {

        # Create a modal dialog to ask the user if they want to fit the signal of
        # the free protein, the free ligand, or the bound complex
        # and to select if they want to fit t0
        # Each of them is a checkboxInput

        showModal(modalDialog(

            tags$h3('Please select the species producing the signal:'),

            tags$script(HTML("
              $(document).on('keypress', function(e) {
                if(e.which == 13) {
                  e.preventDefault();
                  $('#submitOBS').click();
                }
              });
            ")),

            checkboxInput('fit_signal_E','Free protein',FALSE),
            checkboxInput('fit_signal_S','Free ligand',FALSE),
            checkboxInput('fit_signal_ES','Complex',TRUE),
            checkboxInput('fit_t0','Fit t0',TRUE),

            footer=tagList(
              actionButton('submitOBS', 'Submit'),
              modalButton('Cancel')
            )
        ))

        return(NULL)
    }

    pyKinetics$submit_fitting_solution(model)

    popUpSuccess("The fitting has been completed.")

    reactives$solution_model <- model

    # IF model is double exponential, include the plot for Kobs
    # tabPanel(HTML("K<sub>obs</sub>"),plotlyOutput("kobs_plot"))
     if (!reactives$k_obs_plot_shown && model == "double") {

        tab <- tabPanel('k_obs',plotlyOutput("kobs_plot"))

        appendTab(
            inputId='tabBoxFitInSolution',
            tab=tab,
            session = session
        )
        reactives$k_obs_plot_shown <- TRUE
     }

    # Add a Table with the fitted params
    tab <- tabPanel("Fitted params (Kinetics)",tableOutput("fittingInfoKinetics"))

    # Append the Tab with the kinetic fitted parameters
    appendTab(
      inputId='tabBoxFitInSolution',
      tab=tab,
      select = FALSE,
      menuName = NULL,
      session = session
    )

    reactives$kinetics_table_shown_sol <- TRUE
    reactives$kinetics_fit_done <- TRUE

})


observeEvent(input$submitIF,{

    #close the modal dialog
    removeModal()

    popUpInfo("Fitting with induced fit model is in progress.
    A new popup will appear when the fitting is done. Please wait...")

    pyKinetics$submit_fitting_solution(
        'one_binding_site_if',
        fit_signal_E = input$fit_signal_E,
        fit_signal_S = input$fit_signal_S,
        fit_signal_ES = input$fit_signal_ES,
        ESint_equals_ES = input$ESint_equals_ES
    )

    popUpSuccess("The fitting has been completed.")

    # Add a Table with the fitted params
    tab <- tabPanel("Fitted params (Kinetics)",tableOutput("fittingInfoKinetics"))

    # Append the Tab with the kinetic fitted parameters
    appendTab(
      inputId='tabBoxFitInSolution',
      tab=tab,
      select = FALSE,
      menuName = NULL,
      session = session
    )

    reactives$kinetics_table_shown_sol <- TRUE
    reactives$kinetics_fit_done <- TRUE

})

observeEvent(input$submitOBS,{

    #close the modal dialog
    removeModal()

    popUpInfo("Fitting with one binding site model is in progress.
    A new popup will appear when the fitting is done. Please wait...")

    pyKinetics$submit_fitting_solution(
        'one_binding_site',
        fit_signal_E = input$fit_signal_E,
        fit_signal_S = input$fit_signal_S,
        fit_signal_ES = input$fit_signal_ES,
        fixed_t0 = !input$fit_t0
    )

    popUpSuccess("The fitting has been completed.")

    # Add a Table with the fitted params
    tab <- tabPanel("Fitted params (Kinetics)",tableOutput("fittingInfoKinetics"))

    # Append the Tab with the kinetic fitted parameters
    appendTab(
      inputId='tabBoxFitInSolution',
      tab=tab,
      select = FALSE,
      menuName = NULL,
      session = session
    )

    reactives$kinetics_table_shown_sol <- TRUE
    reactives$kinetics_fit_done <- TRUE

})

