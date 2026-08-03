box::use(
  shiny[
    actionButton,
    checkboxInput,
    column,
    conditionalPanel,
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
    tagList
  ]
)

#' @export
visualizationConfigUI <- function(id) {
  ns <- NS(id)
  tagList(
    column(
      width = 3,
      p(
        HTML('<p style="margin-bottom:0px;"></p>'),
        actionButton(
          inputId = ns("configure_visualization"),
          label = "Visualization",
          icon("wrench"),
          style = "color: #fff; background-color: #337ab7;border-color: #2e6da4"
        )
      )
    )
  )
}

#' @export
visualizationConfigServer <- function(id, state, plot_config) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns


    observeEvent(input$configure_visualization,
      {
        req(state$traces_loaded)

        showModal(modalDialog(
          h3("Configure visualization settings"),
          fluidRow(
            column(
              width = 12,
              p(
                HTML("<b>Max points per subplot</b>"),
                sliderInput(
                  ns("visualization_max_points"), NULL,
                  value = plot_config()$max_points, min = 500, max = 10000, step = 750
                )
              )
            ),
            column(
              width = 6,
              p(
                HTML("<b>Smooth curves</b>"),
                checkboxInput(
                  ns("visualization_smooth_curves"), NULL,
                  value = isTRUE(plot_config()$smooth_curves)
                )
              )
            ),
            conditionalPanel(
              condition = "input.visualization_smooth_curves",
              ns = ns,
              column(
                width = 6,
                p(
                  HTML("<b>Rolling window size (s)</b>"),
                  sliderInput(
                    ns("visualization_rolling_window_size"), NULL,
                    value = plot_config()$rolling_window_size,
                    min = 0.5,
                    max = 5,
                    step = 0.5
                  )
                )
              )
            ),
            column(
              width = 6,
              p(
                HTML("<b>Split by smax</b>"),
                checkboxInput(
                  ns("visualization_split_by_smax"), NULL,
                  value = isTRUE(plot_config()$split_by_smax)
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

    # Observers to sync visualization inputs into reactives
    observeEvent(input$visualization_max_points, {
      val <- as.integer(input$visualization_max_points)
      if (is.na(val) || val < 1) {
        req(FALSE)
      }
      val <- min(10000, max(500, val))
      cfg <- plot_config()
      cfg$max_points <- val
      plot_config(cfg)
    })

    observeEvent(input$visualization_smooth_curves, {
      cfg <- plot_config()
      cfg$smooth_curves <- isTRUE(input$visualization_smooth_curves)
      plot_config(cfg)
    })

    observeEvent(input$visualization_rolling_window_size, {
      val <- as.numeric(input$visualization_rolling_window_size)
      if (is.na(val) || val <= 0) {
        req(FALSE)
      }
      val <- min(5, max(0.5, val))
      cfg <- plot_config()
      cfg$rolling_window_size <- val
      plot_config(cfg)
    })

    observeEvent(input$visualization_split_by_smax, {
      cfg <- plot_config()
      cfg$split_by_smax <- isTRUE(input$visualization_split_by_smax)
      plot_config(cfg)
    })
  })
}
