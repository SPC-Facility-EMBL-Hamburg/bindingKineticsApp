box(title = "3. Plotting", width = 12, solidHeader = T, status = "primary",

    fluidRow(

        column(2, p(HTML("<b>File type</b>"),
                    selectInput("plot_type_fit", NULL,
                                c("PNG"    = "png",
                                  "SVG"    = "svg",
                                  "JPEG"    = "jpeg")))),

        column(2, p(HTML("<b>Width</b>"),
                    span(shiny::icon("info-circle"), id = "info_uuPlotWidth"),
                    numericInput('plot_width_fit',NULL, 12,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uuPlotWidth",
                                      tooltip = "Units are pixels * 50",placement = "right"))),

        column(2, p(HTML("<b>Height</b>"),
                    span(shiny::icon("info-circle"), id = "info_uuPlotHeight"),
                    numericInput('plot_height_fit',NULL, 11,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uuPlotHeight",
                                      tooltip = "Units are pixels * 50",
                                      placement = "right"))),

        column(2, p(HTML("<b>Text size</b>"),
                    numericInput('plot_axis_size_fit',NULL, 18,min = 4, max = 40))),

        column(2, p(HTML("<b>Show X-grid</b>")),checkboxInput('plot_show_grid_x_fit',NULL,FALSE)),

        column(2, p(HTML("<b>Show Y-grid</b>")),checkboxInput('plot_show_grid_y_fit',NULL,FALSE))

    ),

    fluidRow(

        column(2, p(HTML("<b>Marker size</b>"),sliderInput('plot_marker_size_fit',NULL, 6,min = 2, max = 12,step = 0.5))),
        column(2, p(HTML("<b>Line width</b>"),sliderInput('plot_line_width_fit',NULL, 3,min = 0.5, max = 10,step = 0.5))),
        conditionalPanel('output.surface_based_binding',
            column(2, p(HTML("<b>Split by Smax ID</b>"),checkboxInput('split_replicates',NULL,TRUE)))
        ),
        column(2, p(HTML("<b>Max points per subplot</b>"),numericInput('max_points_per_plot',NULL, 1500,min = 500, max = 1e4,step = 500))),
        column(2, p(HTML("<b>Smooth curves</b>"),
            span(shiny::icon("info-circle"), id = "info_uuSmoothCurves"),
            checkboxInput('smooth_curves_fit',NULL,FALSE)),
            tippy::tippy_this(elementId = "info_uuSmoothCurves",
                              tooltip = "Smooth curves using a median filter.
                              This affects only the visualization. The raw data is used for fitting.
                              ",placement = "right")),
        column(2, p(HTML("<b>Smooth window size (sec)</b>"),numericInput('median_window_size',NULL, 0.2,min = 0, max = 10,step = 0.1)))

    )

)