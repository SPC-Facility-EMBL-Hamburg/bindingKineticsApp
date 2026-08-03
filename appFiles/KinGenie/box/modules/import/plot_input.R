box::use(
  .. / .. / helpers[
    pandas_to_r
  ],
  .. / .. / plot_functions[
    plot_plate_info,
    plot_traces_all,
    plotlyOutput600px
  ],
  DT[
    datatable,
    DTOutput,
    renderDT
  ],
  plotly[
    plotlyOutput,
    renderPlotly
  ],
  shiny[
    column,
    moduleServer,
    NS,
    observeEvent,
    renderUI,
    req,
    tabPanel,
    tagList,
    uiOutput
  ],
  shinydashboard[
    tabBox
  ]
)

#' @export
plotInputUI <- function(id) {
  ns <- NS(id)

  tagList(
    uiOutput(ns("plot_panel"))
  )
}

#' @export
plotInputServer <- function(id, state, pyKinetics,
                            legend_df_reac = NULL,
                            plot_config_reac = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$plot_panel <- renderUI({
      req(state$traces_loaded)

      if (state$surface_based_binding) {
        tabs <- list(
          tabPanel(
            "Traces",
            plotlyOutput600px(ns("traces"))
          ),
          tabPanel(
            "Assoc. and Diss. traces",
            plotlyOutput600px(ns("tracesAssDiss"))
          ),
          tabPanel(
            "Steps",
            DTOutput(ns("stepsInfo"))
          )
        )

        for (exp in pyKinetics$experiments) {
          if (isTRUE(exp$sample_plate_loaded)) {
            tabs <- c(
              tabs,
              list(
                tabPanel(
                  paste0("Sample plate: ", exp$name),
                  plotlyOutput600px(ns(paste0("samplePlatePlot_", exp$name)))
                )
              )
            )
          }
        }

        do.call(
          tabBox,
          c(
            list(
              title = "",
              width = 12,
              id = ns("tabBoxImport")
            ),
            tabs
          )
        )
      } else {
        tabBox(
          title = "",
          width = 12,
          id = ns("tabBoxImportInSolution"),
          tabPanel(
            "Traces",
            plotlyOutput600px(ns("tracesInSolution"))
          )
        )
      }
    })

    observeEvent(state$traces_loaded,
      {
        req(state$surface_based_binding)

        lapply(pyKinetics$experiments, function(exp) {
          if (!isTRUE(exp$sample_plate_loaded)) {
            req(FALSE)
          }

          local({
            pySingleExp <- exp
            plot_output <- paste0("samplePlatePlot_", pySingleExp$name)

            output[[plot_output]] <- renderPlotly({
              req(plot_config_reac())

              plot_plate_info(
                pySingleExp$sample_row,
                pySingleExp$sample_column,
                pySingleExp$sample_type,
                pySingleExp$sample_id,
                pySingleExp$sample_conc_labeled,
                pySingleExp$name,
                plot_config = plot_config_reac()
              )
            })
          })
        })
      },
      ignoreInit = FALSE
    )

    output$traces <- renderPlotly({
      req(state$traces_loaded)
      req(state$surface_based_binding)
      req(legend_df_reac())

      legend_df <- legend_df_reac()

      xs_all <- lapply(pyKinetics$experiments, \(x) x$xs)
      ys_all <- lapply(pyKinetics$experiments, \(x) x$ys)

      plot_traces_all(
        xs_all,
        ys_all,
        legend_df$Legend,
        legend_df$Color,
        legend_df$Show,
        plot_config = plot_config_reac()
      )
    })

    output$stepsInfo <- renderDT({
      req(state$traces_loaded)
      req(state$surface_based_binding)

      n <- length(pyKinetics$experiment_names)

      dfs <- lapply(pyKinetics$experiments, function(exp) {
        py_df <- exp$df_steps
        df <- pandas_to_r(py_df)

        if (n > 1) df["Name"] <- exp$name

        df
      })

      df <- do.call(rbind, dfs)

      datatable(df, options = list(scrollY = "400px", paging = FALSE), rownames = FALSE)
    })

    output$tracesAssDiss <- renderPlotly({
      req(state$traces_loaded)
      req(legend_df_reac())
      req(plot_config_reac())
      req(state$surface_based_binding)

      legend_df <- legend_df_reac()

      xs_all <- list()
      ys_all <- list()

      cnt <- 1
      for (exp in pyKinetics$experiments) {
        py_df <- exp$df_steps
        df <- pandas_to_r(py_df)

        selected <- grepl("ASS", df$Type)

        xs <- exp$xs
        ys <- exp$ys

        xs_subset <- lapply(xs, function(x) x[selected])
        ys_subset <- lapply(ys, function(y) y[selected])

        xs_all[[cnt]] <- xs_subset
        ys_all[[cnt]] <- ys_subset

        cnt <- cnt + 1
      }

      fig <- plot_traces_all(
        xs_all,
        ys_all,
        legend_df$Legend,
        legend_df$Color,
        legend_df$Show,
        plot_config = plot_config_reac()
      )

      fig
    })

    output$tracesInSolution <- renderPlotly({
      req(state$traces_loaded)
      req(legend_df_reac())
      req(!state$surface_based_binding)

      legend_df <- legend_df_reac()

      xs_all <- lapply(pyKinetics$experiments, function(x) x$xs)
      ys_all <- lapply(pyKinetics$experiments, function(x) x$ys)

      fig <- plot_traces_all(
        xs_all,
        ys_all,
        legend_df$Legend,
        legend_df$Color,
        legend_df$Show,
        plot_config = plot_config_reac()
      )

      fig
    })
  })
}
