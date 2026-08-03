box::use(
  shiny[
    actionButton,
    checkboxInput,
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
    showModal,
    sliderInput,
    tagList,
    textInput
  ]
)

#' @export
axisConfigUI <- function(id) {
  ns <- NS(id)
  tagList(
    column(
      width = 3,
      p(
        HTML('<p style="margin-bottom:0px;"></p>'),
        actionButton(
          inputId = ns("configure_axis"),
          label = "Axis",
          icon("grip"),
          style = "color: #fff; background-color: #337ab7;
              border-color: #2e6da4"
        )
      )
    )
  )
}

#' @export
axisConfigServer <- function(id, plot_config) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$configure_axis,
      {
        showModal(modalDialog(
          h3("Configure axis settings"),
          fluidRow(
            column(
              width = 12,
              p(
                HTML("<b>X-axis label</b>"),
                textInput(
                  inputId = ns("x_axis_label_input"),
                  label = NULL,
                  value = plot_config()$x_label
                )
              )
            ),
            column(
              width = 12,
              p(
                HTML("<b>Y-axis label</b>"),
                textInput(
                  inputId = ns("y_axis_label_input"),
                  label = NULL,
                  value = plot_config()$y_label
                )
              )
            ),
            column(
              width = 12,
              p(
                HTML("<b>Number of ticks</b>"),
                column(
                  6,
                  sliderInput(
                    inputId = ns("axis_n_xticks_input"),
                    label = "X-axis",
                    value = plot_config()$n_xticks,
                    min = 3,
                    max = 10,
                    step = 1
                  )
                ),
                column(
                  6,
                  sliderInput(
                    inputId = ns("axis_n_yticks_input"),
                    label = "Y-axis",
                    value = plot_config()$n_yticks,
                    min = 3,
                    max = 10,
                    step = 1
                  )
                )
              )
            ),
            column(
              width = 6,
              p(
                HTML("<b>Tick length (px)</b>"),
                sliderInput(
                  inputId = ns("axis_tick_length_input"),
                  label = NULL,
                  value = plot_config()$tick_length,
                  min = 4,
                  max = 16,
                  step = 1
                )
              )
            ),
            column(
              width = 6,
              p(
                HTML("<b>Tick width (px)</b>"),
                sliderInput(
                  inputId = ns("axis_tick_width_input"),
                  label = NULL,
                  value = plot_config()$tick_width,
                  min = 0,
                  max = 10,
                  step = 1
                )
              )
            ),
            column(
              width = 12,
              p(
                HTML("<b>Show grid</b>"),
                column(
                  6,
                  checkboxInput(
                    inputId = ns("x_axis_show_grid_input"),
                    label = "x-axis",
                    value = isTRUE(plot_config()$show_grid_x)
                  )
                ),
                column(
                  6,
                  checkboxInput(
                    inputId = ns("y_axis_show_grid_input"),
                    label = "y-axis",
                    value = isTRUE(plot_config()$show_grid_y)
                  )
                )
              )
            )
          ),
          size = "m",
          easyClose = TRUE,
          footer = tagList(
            modalButton("Close")
          )
        ))
      },
      ignoreInit = TRUE
    )

    # Observers to update reactives when the modal inputs change
    observeEvent(input$axis_n_xticks_input, {
      val <- as.integer(input$axis_n_xticks_input)
      if (is.na(val) || val < 1) {
        req(FALSE)
      }
      val <- min(10, max(3, val))
      cfg <- plot_config()
      cfg$n_xticks <- val
      plot_config(cfg)
    })

    observeEvent(input$axis_n_yticks_input, {
      val <- as.integer(input$axis_n_yticks_input)
      if (is.na(val) || val < 1) {
        req(FALSE)
      }
      val <- min(10, max(3, val))
      cfg <- plot_config()
      cfg$n_yticks <- val
      plot_config(cfg)
    })

    observeEvent(input$axis_tick_length_input, {
      val <- as.numeric(input$axis_tick_length_input)
      if (is.na(val)) {
        req(FALSE)
      }
      val <- min(16, max(4, val))
      cfg <- plot_config()
      cfg$tick_length <- val
      plot_config(cfg)
    })

    observeEvent(input$axis_tick_width_input, {
      val <- as.numeric(input$axis_tick_width_input)
      if (is.na(val)) {
        req(FALSE)
      }
      val <- min(10, max(0, val))
      cfg <- plot_config()
      cfg$tick_width <- val
      plot_config(cfg)
    })

    observeEvent(input$y_axis_show_grid_input, {
      cfg <- plot_config()
      cfg$show_grid_y <- isTRUE(input$y_axis_show_grid_input)
      plot_config(cfg)
    })

    observeEvent(input$x_axis_show_grid_input, {
      cfg <- plot_config()
      cfg$show_grid_x <- isTRUE(input$x_axis_show_grid_input)
      plot_config(cfg)
    })

    observeEvent(input$x_axis_label_input, {
      cfg <- plot_config()
      cfg$x_label <- input$x_axis_label_input
      plot_config(cfg)
    })

    observeEvent(input$y_axis_label_input, {
      cfg <- plot_config()
      cfg$y_label <- input$y_axis_label_input
      plot_config(cfg)
    })
  })
}
