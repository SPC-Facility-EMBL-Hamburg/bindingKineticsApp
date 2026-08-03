box::use(
  colourpicker[
    colourInput
  ],
  shiny[
    actionButton,
    column,
    fluidRow,
    h3,
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
    tagList
  ]
)

#' @export
appearanceConfigUI <- function(id, width = 5) {
  ns <- NS(id)

  tagList(
    column(
      width = width,
      p(
        HTML('<p style="margin-bottom:0px;"></p>'),
        actionButton(
          inputId = ns("configure_appearance"),
          label = "Font, markers, colors",
          icon("marker"),
          style = "color: #fff; background-color: #337ab7;
              border-color: #2e6da4"
        )
      )
    )
  )
}

#' @export
appearanceConfigServer <- function(
    id,
    state,
    plot_config,
    show_colour_controls = TRUE,
    legend_df_reac = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Appearance modal: marker size, font size, line width, max points
    observeEvent(input$configure_appearance,
      {
        req(state$traces_loaded)

        if (show_colour_controls) {
          req(legend_df_reac())
          legend <- legend_df_reac()
          ids <- as.character(legend$Internal_ID)
        }

        margin_bottom <- ifelse(show_colour_controls, "110px", "10px")

        showModal(modalDialog(
          h3("Configure appearance settings"),
          fluidRow(
            style = paste0("margin-bottom: ", margin_bottom, ";"),
            column(
              width = 6,
              p(
                HTML("<b>Marker size</b>"),
                sliderInput(
                  ns("appearance_marker_size"), NULL,
                  value = plot_config()$marker_size, min = 1, max = 20, step = 0.5
                )
              )
            ),
            column(
              width = 6,
              p(
                HTML("<b>Line width</b>"),
                sliderInput(
                  ns("appearance_line_width"), NULL,
                  value = plot_config()$line_width, min = 0.5, max = 10, step = 0.5
                )
              )
            ),
            column(
              width = 12,
              p(
                HTML("<b>Font size</b>"),
                sliderInput(
                  ns("appearance_font_size"), NULL,
                  value = plot_config()$axis_size, min = 8, max = 34, step = 1
                )
              )
            ),
            if (show_colour_controls) {
              column(
                12,
                column(
                  width = 6,
                  p(
                    HTML("<b>Set colour</b>"),
                    selectInput(
                      inputId = ns("mol2changeColor"),
                      label = NULL,
                      choices = ids,
                      selected = state$selected_trace,
                      selectize = FALSE
                    )
                  )
                ),
                column(
                  width = 6,
                  p(
                    HTML('<p style="margin-bottom:0px;"><br></p>'),
                    colourInput(ns("colorForLegend"), NULL,
                      value = state$selected_color
                    )
                  )
                )
              )
            }
          ),
          easyClose = TRUE,
          footer = tagList(
            modalButton("Close")
          )
        ))
      },
      ignoreInit = TRUE
    )

    # Observers to sync appearance inputs into reactives
    observeEvent(input$appearance_marker_size, {
      val <- as.numeric(input$appearance_marker_size)
      if (is.na(val)) {
        req(FALSE)
      }
      val <- max(0.1, min(100, val))
      cfg <- plot_config()
      cfg$marker_size <- val
      plot_config(cfg)
    })

    observeEvent(input$appearance_font_size, {
      val <- as.numeric(input$appearance_font_size)
      if (is.na(val)) {
        req(FALSE)
      }
      val <- max(8, min(34, round(val)))
      cfg <- plot_config()
      cfg$axis_size <- val
      plot_config(cfg)
    })


    observeEvent(input$appearance_line_width, {
      val <- as.numeric(input$appearance_line_width)
      if (is.na(val)) {
        req(FALSE)
      }
      val <- max(0.1, min(100, val))
      cfg <- plot_config()
      cfg$line_width <- val
      plot_config(cfg)
    })

    observeEvent(input$colorForLegend, {
      df <- legend_df_reac()

      idx <- which(df$Internal_ID == input$mol2changeColor)
      df$Color[idx] <- input$colorForLegend

      legend_df_reac(df)

      state$selected_color <- input$colorForLegend
      state$selected_trace <- input$mol2changeColor
    })
  })
}
