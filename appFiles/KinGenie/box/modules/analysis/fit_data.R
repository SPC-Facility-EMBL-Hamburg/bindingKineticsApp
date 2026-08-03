box::use(
  .. / .. / helpers[
    pandas_to_r
  ],
  .. / .. / tables[
    render_ligand_conc_df_solution,
    render_ligand_conc_df_surface
  ],
  rhandsontable[
    hot_to_r,
    renderRHandsontable,
    rHandsontableOutput
  ],
  shiny[
    actionButton,
    checkboxInput,
    column,
    conditionalPanel,
    fluidRow,
    HTML,
    modalButton,
    modalDialog,
    moduleServer,
    NS,
    observeEvent,
    p,
    reactiveVal,
    removeModal,
    renderText,
    renderUI,
    req,
    selectInput,
    showModal,
    tagList,
    tags,
    updateCheckboxInput
  ],
  shinydashboard[
    box
  ]
)

#' @export
fitDataUI <- function(id) {
  ns <- NS(id)

  tagList(
    box(
      title = "1. Dataset info", width = 12, solidHeader = T, status = "primary",
      fluidRow(
        # Add checkbox to show or hide the ligandInfo dataframe
        column(
          4,
          checkboxInput(ns("showLigandInfo"), "Show dataset info", value = TRUE)
        ),

        # Add button for quick selection of the ligand info
        column(4, p(
          HTML("<br>"),
          actionButton(ns("quickSelect"), "Quick select")
        )),
        conditionalPanel(
          "input.showLigandInfo",
          ns = ns,
          column(12, rHandsontableOutput(ns("ligandInfo")))
        )
      )
    )
  )
}

#' @export
fitDataServer <- function(id, state, pyKinetics) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ligandInfo <- reactiveVal()

    observeEvent(state$ligand_info_version,
      {
        req(state$traces_loaded)

        if (state$surface_based_binding) {
          pyKinetics$merge_ligand_conc_df()

          py_df <- pyKinetics$combined_ligand_conc_df

          df <- pandas_to_r(py_df)

          ligandInfo(df)
        } else {
          pyKinetics$merge_conc_df_solution()

          py_df <- pyKinetics$combined_conc_df
          df <- pandas_to_r(py_df)

          ligandInfo(df)
        }
      },
      ignoreInit = TRUE
    )

    output$ligandInfo <- renderRHandsontable({
      req(ligandInfo())

      if (state$surface_based_binding) {
        render_ligand_conc_df_surface(ligandInfo())
      } else {
        render_ligand_conc_df_solution(ligandInfo())
      }
    })

    observeEvent(input$ligandInfo,
      {
        ligandInfo(hot_to_r(input$ligandInfo))
      },
      ignoreInit = TRUE
    )

    observeEvent(input$quickSelect, {
      df <- hot_to_r(input$ligandInfo)

      # Find the different sample IDs
      sample_ids <- unique(df$SampleID)

      # Option for surface-based binding
      if (state$surface_based_binding) {
        # Find different Smax IDs
        smax_ids <- c("any", unique(df$Smax_ID))

        # Show the modal dialog with the sample IDs
        showModal(modalDialog(
          tags$h4("Select the sample IDs"),
          selectInput(ns("sample_ids_fit"), "Sample IDs", choices = sample_ids, multiple = TRUE),
          tags$h4(),
          tags$h4("Select the Smax IDs"),
          selectInput(ns("smax_ids_fit"), "Smax IDs", choices = smax_ids, multiple = TRUE),
          footer = tagList(
            actionButton(ns("quickSelectOK"), "OK"),
            modalButton("Cancel")
          )
        ))
      } else {
        # Show the modal dialog with the sample IDs
        showModal(modalDialog(
          tags$h4("Select the sample IDs"),
          selectInput(ns("sample_ids_fit"), "Sample IDs",
            choices = sample_ids, multiple = TRUE
          ),
          footer = tagList(
            actionButton(ns("quickSelectOK"), "OK"),
            modalButton("Cancel")
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

      if (state$surface_based_binding) {
        selected_smax_ids <- input$smax_ids_fit

        if (!is.null(selected_smax_ids) & !("any" %in% selected_smax_ids)) {
          df$Select <- df$Select & ifelse(df$Smax_ID %in% selected_smax_ids, TRUE, FALSE)
        }
      }

      ligandInfo(df)

      # Update show dataset info
      updateCheckboxInput(session, "showLigandInfo", value = TRUE)
    })

    observeEvent(state$show_ligand_info,
      {
        updateCheckboxInput(session, "showLigandInfo", value = state$show_ligand_info)
      },
      ignoreInit = TRUE
    )

    list(
      df = ligandInfo
    )
  })
}
