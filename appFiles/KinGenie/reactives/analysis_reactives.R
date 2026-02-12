removeDynamicTabs <- function() {

    # Remove all the tabs that we dynamically create
    if (reactives$kinetics_table_shown) {

        removeTab(
          inputId='tabBoxFit',
          target = 'Fitted params (Kinetics)',
          session = session
        )

        reactives$kinetics_table_shown <- FALSE

    }

    if (reactives$ss_plot_shown) {

        removeTab(
          inputId='tabBoxFit',
          target = "Steady-state",
          session = session
        )

        reactives$ss_plot_shown <- FALSE

    }

    if (reactives$ss_table_shown) {

        removeTab(
          inputId='tabBoxFit',
          target = 'Fitted params (Steady-state)',
          session = session
        )

        reactives$ss_table_shown <- FALSE

    }

    if (reactives$kinetics_ci95_table_shown) {

        removeTab(
          inputId='tabBoxFit',
          target = 'Asymmetric error',
          session = session
        )

        reactives$kinetics_ci95_table_shown <- FALSE

    }

    if (reactives$diagnostic_plots_done) {

        removeTab(
          inputId='tabBoxFit',
          target = 'Observed constants',
          session = session
        )

        reactives$diagnostic_plots_done <- FALSE

    }

    return(NULL)

}

observeEvent(input$quickSelect,{

    df <- hot_to_r(input$ligandInfo)

    # Find the different sample IDs
    sample_ids <- unique(df$SampleID)

    # Option for surface-based binding
    if (reactives$surface_based_binding) {

        # Find different Smax IDs
        smax_ids <- c('any',unique(df$Smax_ID))

        # Show the modal dialog with the sample IDs
        showModal(modalDialog(
            tags$h4("Select the sample IDs"),
            selectInput("sample_ids_fit", "Sample IDs", choices = sample_ids, multiple = TRUE),
            tags$h4(),
            tags$h4("Select the Smax IDs"),
            selectInput("smax_ids_fit", "Smax IDs", choices = smax_ids, multiple = TRUE),

            footer = tagList(
                actionButton("quickSelectOK", "OK"),
                modalButton('Cancel')
            )
        ))

    } else {

        # Show the modal dialog with the sample IDs
        showModal(modalDialog(
            tags$h4("Select the sample IDs"),
            selectInput("sample_ids_fit", "Sample IDs", choices = sample_ids, multiple = TRUE),

            footer = tagList(
                actionButton("quickSelectOK", "OK"),
                modalButton('Cancel')
            )
        ))

    }

})


observeEvent(input$quickSelectOK, {

    removeModal()

    df <- hot_to_r(input$ligandInfo)

    # Get the selected sample IDs
    selected_sample_ids <- input$sample_ids_fit

    # Check all if it is NULL
    if (is.null(selected_sample_ids)) {
        selected_sample_ids <- unique(df$SampleID)
    }

    # Set the 'Select' column to TRUE for the selected sample IDs
    df$Select <- ifelse(df$SampleID %in% selected_sample_ids, TRUE, FALSE)

    if (reactives$surface_based_binding) {

        selected_smax_ids <- input$smax_ids_fit

        if (!is.null(selected_smax_ids) & !('any' %in% selected_smax_ids)) {
            df$Select <- df$Select & ifelse(df$Smax_ID %in% selected_smax_ids, TRUE, FALSE)
        }

        # Update the ligandInfo table with the filtered dataframe
        output$ligandInfo <- renderRHandsontable({render_combined_ligand_conc_df(df)})

    } else {

       output$ligandInfo <- renderRHandsontable({render_combined_ligand_conc_df_solution(df)})

    }



    # Update show dataset info
    updateCheckboxInput(session, "showLigandInfo", value = TRUE)

})

observeEvent(input$triggerCreateDataset, {

    req(reactives$traces_loaded)
    removeDynamicTabs()

    df <- hot_to_r(input$ligandInfo)

    # First filter, verify that we do not have any negative concentration in rows where the 'Select' column is TRUE
    sel_concs <- df[,2][df$Select] # The second column is the ligand concentration, both for surface-based and solution-based binding

    if (any(sel_concs < 0)) {
        popUpWarning("Negative ligand concentrations found in selected rows.
        Please correct them before creating the dataset.")
        return(NULL)
    }

    # Set the checkboxInput showLigandInfo to FALSE
    if (input$showLigandInfo) {
        updateCheckboxInput(session, "showLigandInfo", value = FALSE)
        Sys.sleep(1)
    }

    reactives$fit_dataset_loaded <- FALSE

    reactives$kinetics_fit_done  <- FALSE
    reactives$ss_fit_done        <- FALSE

    pyKinetics$init_fittings()

    append_record_to_logbook('Dataset creation triggered',include_time = TRUE,add_empty_line = TRUE)
    append_record_to_logbook(df_to_lines(df))

    logbook_messages <- pyKinetics$generate_fittings(df)

    for (msg in logbook_messages) {

        append_record_to_logbook(msg,add_empty_line = TRUE)

    }

    reactives$is_single_cycle <- any(pyKinetics$get_experiment_properties('is_single_cycle',fittings=TRUE))

    if (!reactives$is_single_cycle) {

        # include the tabPanel("Steady-state",plotlyOutput("steady_state"))
        tab <- tabPanel("Steady-state",plotlyOutput("steady_state"))
        # Append the Tab with the steady-state plot
        insertTab(
          inputId='tabBoxFit',
          tab=tab,
          session = session,
          target = "Assoc. and diss. traces",
          position = "after"
        )
        reactives$ss_plot_shown <- TRUE

    } else {

        # remove  the option to fit the steady-state region
        updateSelectInput(session, "fittingRegion",
            choices = c(
                'Association and dissociation' = 'association_dissociation',
                'Dissociation'                   = 'dissociation'
            )
        )

    }

    reactives$fit_dataset_loaded    <- TRUE

    updateCheckboxInput(session, "showLigandInfo", value = FALSE)

    popUpInfo("Dataset(s) created.
    Traces with the same sample ID will share thermodynamic parameters
    (e.g., <i>K</i><sub>d</sub> and <i>k</i><sub>off</sub>).
    Traces with the same sample ID and Smax ID will share the <i>Smax</i> parameter,
    if the option linked by Smax is set to TRUE.")

})

observeEvent(input$fittingRegion, {

    if (input$fittingRegion == 'association_dissociation') {

        updateSelectInput(session, "fittingModel",
        choices = c('One-to-one' = 'one_to_one',
        'One-to-one (MTL)' = 'one_to_one_mtl',
        'One-to-one (induced fit)' = 'one_to_one_if'
        ))

    } else {

        updateSelectInput(session, "fittingModel",
        choices = c('One-to-one' = 'one_to_one'))

    }

})

observeEvent(input$triggerFitting, {

    req(reactives$fit_dataset_loaded)

    removeDynamicTabs()

    reactives$kinetics_fit_done <- FALSE
    reactives$ss_fit_done       <- FALSE

    reactives$fit_dataset_loaded <- FALSE

    result <- tryCatch(
        {

            pyKinetics$submit_steady_state_fitting()

        }, error = function(e) {
            if (inherits(e, "python.builtin.RuntimeError")) {
                err <- py_last_error()
                popUpWarning(
                    paste0("⚠ Fitting error: ", err$value)
                )
                return('Error')
            } else {
                stop(e) # rethrow non-Python errors
            }
        }
    )

    if (!is.null(result)) return(NULL)

    fittingRegion <- input$fittingRegion

    if (fittingRegion == 'steady_state') {

        popUpSuccess("Steady state fitting done.
        Remember that the model assumes
        equal binding capacity for all the samples with the same combination of sample ID,
        and Smax ID.")

        Sys.sleep(0.5)

        append_record_to_logbook('Steady state fitting done.',add_empty_line = TRUE)

        reactives$ss_fit_done         <- TRUE
        reactives$fit_dataset_loaded  <- TRUE

        reactives$kinetics_fit_done   <- FALSE

        # include the tabPanel("Steady-state",plotlyOutput("steady_state"))
        tab <- tabPanel("Steady-state",plotlyOutput("steady_state"))
        # Append the Tab with the steady-state plot
        insertTab(
          inputId='tabBoxFit',
          tab=tab,
          session = session,
          target = "Assoc. and diss. traces",
          position = "after"
        )
        reactives$ss_plot_shown <- TRUE

        reactives$ss_table_shown <- TRUE
        tab <- tabPanel("Fitted params (Steady-state)",tableOutput("fittingInfoSS"))

        # Append the Tab with the steady-state fitted parameters
        insertTab(
          inputId='tabBoxFit',
          tab=tab,
          session = session,
          target = "Steady-state",
          position = "after"
        )

    } else {

        showModal(modalDialog(

            tags$h4(
            'Can you assume equal binding capacity for all the samples with the same combination
            of sample ID, experiment and Smax ID?
            If yes, set the option below to TRUE.
            Hint: If the loading levels are different, you should set this option to FALSE.
            Only if the same sensor was used for all the ligand concentrations in a series,
            you should set this option to TRUE.'),

            checkboxInput('linkedRmax','Smax Linked By Sensor',FALSE),

            footer=tagList(
              actionButton('submitKineticsFitting', 'Submit'),
              actionButton('cancelKineticsFitting', 'Cancel')
            )
        ))

    }

})

observeEvent(input$submitKineticsFitting, {

    removeModal()

    reactives$kinetics_fit_done <- FALSE

    fittingModel  <- input$fittingModel
    fittingRegion <- input$fittingRegion
    shared_smax   <- input$linkedRmax

    result <- tryCatch(
        {

            pyKinetics$submit_kinetics_fitting(fittingModel,fittingRegion,shared_smax)

        }, error = function(e) {
            if (inherits(e, "python.builtin.RuntimeError")) {
                err <- py_last_error()
                popUpWarning(
                    paste0("⚠ Fitting error: ", err$value)
                )
                return('Error')
            } else {
                stop(e) # rethrow non-Python errors
            }
        }
    )

    if (!is.null(result)) return(NULL)

    popUpSuccess("Fitting done.
    The fitted parameters are shown in the 'Fitted params (Kinetics)' table.
    For the fitting algorithm to work, we use boundaries.
    Check them in the 'Fitted params (boundaries)' table. If a fitted parameter is too close to the
    lower or upper boundary, it will be highlighted in red.")

    append_record_to_logbook('Kinetics fitting done.',add_empty_line = TRUE)
    append_record_to_logbook(paste0('Region:', fittingRegion),add_empty_line = FALSE)
    append_record_to_logbook(paste0('Model:', fittingModel),add_empty_line = FALSE)

    append_record_to_logbook(paste0('Linked Rmax:',input$linkedRmax))

    reactives$kinetics_fit_done   <- TRUE
    reactives$fit_dataset_loaded  <- TRUE

    reactives$kinetics_table_shown <- TRUE
    tab <- tabPanel("Fitted params (Kinetics)",tableOutput("fittingInfoKinetics"))

    # Append the Tab with the kinetic fitted parameters
    insertTab(
      inputId='tabBoxFit',
      tab=tab,
      session = session,
      target = "Assoc. and diss. traces",
      position = "after"
    )

    c1 <- fittingModel == 'one_to_one'
    c2 <- grepl('asso',fittingRegion)
    #c3 <- !reactives$is_single_cycle

    if (c1 & c2) {

        tab <- tabPanel("Observed constants",plotlyOutput("diagnostic_plot"))
        # Append the Tab with the observed constants
        appendTab(
          inputId  = 'tabBoxFit',
          tab      = tab,
          select   = FALSE,
          menuName = NULL,
          session  = session
        )

        reactives$diagnostic_plots_done <- TRUE
    }

})

observeEvent(input$cancelKineticsFitting, {
    removeModal()
    # Add any additional actions to be taken when Cancel is clicked
    reactives$kinetics_fit_done   <- FALSE
    reactives$fit_dataset_loaded  <- TRUE
    popUpInfo("Kinetics fitting was cancelled.")
})

# Show or hide the confidence interval button
observeEvent(list(input$fittingRegion,input$fittingModel), {

    req(reactives$fit_dataset_loaded)

    # Show the button if the fitting region is association and dissociation
    # and the selected model is one to one
    if (input$fittingRegion == 'association_dissociation' &&
        input$fittingModel  == 'one_to_one') {

        shinyjs::show("triggerAsymError")

    } else {

        shinyjs::hide("triggerAsymError")

    }

})

# Calculate the asymmetric error
observeEvent(input$triggerAsymError, {

    req(reactives$fit_dataset_loaded)
    req(reactives$kinetics_fit_done)

    popUpInfo("Calculating the asymmetric error for K_d and k_off.
    Please wait some minutes for the calculations to finish.")

    pyKinetics$calculate_asymmetric_error()

    # Create the tab with the asymmetric error

    reactives$kinetics_ci95_table_shown <- TRUE
    tab <- tabPanel("Asymmetric error",tableOutput("fittingInfoKineticsCI95"))

    # Append the Tab with the asymmetric error
    insertTab(
      inputId='tabBoxFit',
      tab=tab,
      session = session,
      target = "Fitted params (Kinetics)",
      position = "after"
    )

    popUpSuccess("Asymmetric error calculated for K_d and k_off.")

    append_record_to_logbook('Asymmetric error calculated for K_d and k_off.',add_empty_line = TRUE)

})


