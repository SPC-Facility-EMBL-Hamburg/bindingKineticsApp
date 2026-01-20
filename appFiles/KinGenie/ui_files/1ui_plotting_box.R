box(title = "3. Plotting", width = 12, solidHeader = T, status = "primary",

    column(6,

        fluidRow(
          column(width = 3,
                 p(HTML("<b><br>Choose ID</b>")),
                 selectInput("mol2changeColor", label = NULL,c("X"))),

          column(width = 3,
                 p(HTML("<b><br>Set colour</b>")),
                 colourpicker::colourInput("colorForLegend",label=NULL, value = "#E41A1C")),

          column(3, p(HTML("<b><br>Show plot export options</b>")),
                      checkboxInput('showAdvancedPlottingOptions',NULL,TRUE))
        ),

        conditionalPanel(

            "input.showAdvancedPlottingOptions",

        fluidRow(

            column(3, p(HTML("<b>File type</b>"),
                        selectInput("plot_type", NULL,
                                    c("PNG"    = "png",
                                      "SVG"    = "svg",
                                      "JPEG"    = "jpeg")))),

            column(3, p(HTML("<b>Width</b>"),
                        span(shiny::icon("info-circle"), id = "info_uuPlotWidth"),
                        numericInput('plot_width',NULL, 12,min = 1, max = 100),
                        tippy::tippy_this(elementId = "info_uuPlotWidth",
                                          tooltip = "Units are pixels * 50",placement = "right"))),

            column(3, p(HTML("<b>Height</b>"),
                        span(shiny::icon("info-circle"), id = "info_uuPlotHeight"),
                        numericInput('plot_height',NULL, 11,min = 1, max = 100),
                        tippy::tippy_this(elementId = "info_uuPlotHeight",
                                          tooltip = "Units are pixels * 50",
                                          placement = "right"))),

            column(3, p(HTML("<b>Text size</b>"),
                        numericInput('plot_axis_size',NULL, 18,min = 4, max = 40)))

        ),

        fluidRow(

               column(3, p(HTML("<b>Show X-grid</b>")),checkboxInput('plot_show_grid_x',NULL,FALSE)),

               column(3, p(HTML("<b>Show Y-grid</b>")),checkboxInput('plot_show_grid_y',NULL,FALSE)),

               column(3, p(HTML("<b>Marker size</b>"),sliderInput('plot_marker_size',NULL, 3,min = 2, max = 8,step = 0.5))),

               column(3, p(HTML("<b>Line width</b>"),sliderInput('plot_line_width',NULL, 2,min = 0.5, max = 6,step = 0.5)))

        ))

    ),

    column(6,

        fluidRow(column(width = 12,rHandsontableOutput('legendInfo')))

    )

)