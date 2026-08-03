box::use(
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
    numericInput,
    observeEvent,
    p,
    req,
    selectInput,
    showModal,
    tagList
  ]
)

#' @export
exportConfigUI <- function(id) {
  ns <- NS(id)
  tagList(
    column(
      width = 3,
      p(
        HTML('<p style="margin-bottom:0px;"></p>'),
        actionButton(
          inputId = ns("configure_export"),
          label = "Export",
          icon("maximize"),
          style = "color: #fff; background-color: #337ab7;
              border-color: #2e6da4"
        )
      )
    )
  )
}

#' @export
exportConfigServer <- function(id, plot_config) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$configure_export,
      {
        showModal(modalDialog(
          h3("Configure export settings"),
          fluidRow(
            column(
              width = 12,
              p(
                HTML("<b>File type</b>"),
                selectInput(
                  inputId = ns("export_file_type_input"),
                  label = NULL,
                  choices = c("PNG" = "png", "SVG" = "svg", "JPEG" = "jpeg"),
                  selected = plot_config()$type,
                )
              )
            ),
            column(
              width = 12,
              p(
                HTML("<b>Plot width (px)</b>"),
                numericInput(
                  inputId = ns("export_plot_width_input"),
                  label = NULL,
                  value = plot_config()$width,
                  min = 100,
                  max = 10000,
                  step = 50
                )
              )
            ),
            column(
              width = 12,
              p(
                HTML("<b>Plot height (px)</b>"),
                numericInput(
                  inputId = ns("export_plot_height_input"),
                  label = NULL,
                  value = plot_config()$height,
                  min = 100,
                  max = 10000,
                  step = 50
                )
              )
            )
          ),
          easyClose = TRUE,
          footer = tagList(
            modalButton("Close")
          )
        ))
      },
      ignoreInit = TRUE
    )

    observeEvent(input$export_file_type_input,
      {
        req(input$export_file_type_input)
        cfg <- plot_config()
        cfg$type <- input$export_file_type_input
        plot_config(cfg)
      },
      ignoreInit = TRUE
    )


    observeEvent(input$export_plot_width_input,
      {
        val <- as.numeric(input$export_plot_width_input)
        if (is.na(val)) {
          req(FALSE)
        }
        val <- min(10000, max(100, round(val)))
        cfg <- plot_config()
        cfg$width <- val
        plot_config(cfg)
      },
      ignoreInit = TRUE
    )

    observeEvent(input$export_plot_height_input,
      {
        val <- as.numeric(input$export_plot_height_input)
        if (is.na(val)) {
          req(FALSE)
        }
        val <- min(10000, max(100, round(val)))
        cfg <- plot_config()
        cfg$height <- val
        plot_config(cfg)
      },
      ignoreInit = TRUE
    )
  })
}
