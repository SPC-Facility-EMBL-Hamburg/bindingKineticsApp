box(title = "Plotting", width = 12, solidHeader = T, status = "primary",

    fluidRow(

        column(2, p(HTML("<b>File type</b>"),
                    selectInput("plot_type_sim", NULL,
                                c("PNG"    = "png",
                                  "SVG"    = "svg",
                                  "JPEG"    = "jpeg")))),

        column(2, p(HTML("<b>Width</b>"),
                    span(shiny::icon("info-circle"), id = "info_uuPlotWidth_sim"),
                    numericInput('plot_width_sim',NULL, 12,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uuPlotWidth_sim",
                                      tooltip = "Units are pixels * 50",placement = "right"))),

        column(2, p(HTML("<b>Height</b>"),
                    span(shiny::icon("info-circle"), id = "info_uuPlotHeight_sim"),
                    numericInput('plot_height_sim',NULL, 11,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uuPlotHeight_sim",
                                      tooltip = "Units are pixels * 50",
                                      placement = "right"))),

        column(2, p(HTML("<b>Text size</b>"),
                    numericInput('plot_axis_size_sim',NULL, 18,min = 4, max = 40))),

        column(2, p(HTML("<b>Show X-grid</b>")),checkboxInput('plot_show_grid_x_sim',NULL,FALSE)),

        column(2, p(HTML("<b>Show Y-grid</b>")),checkboxInput('plot_show_grid_y_sim',NULL,FALSE))

    ),

    fluidRow(

           column(2, p(HTML("<b>Marker size</b>"),sliderInput('plot_marker_size_sim',NULL, 6,min = 2, max = 12,step = 0.5))),

    )

)