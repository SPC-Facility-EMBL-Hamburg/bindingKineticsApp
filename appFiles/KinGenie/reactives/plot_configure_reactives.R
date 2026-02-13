# Export configuration modal: file type, plot width, plot height
observeEvent(list(input$configure_export,input$configure_export_fit), {
    showModal(modalDialog(
        h3("Configure export settings"),
        fluidRow(

            column(
                width = 12,

                p(HTML("<b>File type</b>"),

                    selectInput(
                        inputId = "export_file_type_input",
                        label = NULL,
                        choices = c("PNG" = "png", "SVG" = "svg", "JPEG" = "jpeg"),
                        selected = reactives$plot_config$type,
                    )
                )
            ),

            column(
                width = 12,

                p(HTML("<b>Plot width (px)</b>"),

                    numericInput(
                        inputId = "export_plot_width_input",
                        label = NULL,
                        value =  reactives$plot_config$width,
                        min = 100,
                        max = 10000,
                        step = 50
                    )
                )
            ),

            column(
                width = 12,

                p(HTML("<b>Plot height (px)</b>"),

                    numericInput(
                        inputId = "export_plot_height_input",
                        label = NULL,
                        value = reactives$plot_config$height,
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
}, ignoreInit = TRUE)

observeEvent(input$export_file_type_input,{
    req(input$export_file_type_input)
    reactives$plot_config$type <- input$export_file_type_input
})


observeEvent(input$export_plot_width_input, {
    val <- as.numeric(input$export_plot_width_input)
    if (is.na(val)) return()
    val <- min(10000, max(100, round(val)))
    reactives$plot_config$width <- val
})

observeEvent(input$export_plot_height_input, {
    val <- as.numeric(input$export_plot_height_input)
    if (is.na(val)) return()
    val <- min(10000, max(100, round(val)))
    reactives$plot_config$height <- val
})

# Axis configuration modal: show controls for x/y ticks, tick length/width, and grid toggles
observeEvent(list(input$configure_axis,input$configure_axis_fit), {
    showModal(modalDialog(
        h3("Configure axis settings"),

        fluidRow(

            column(
                width = 12,

                p(HTML("<b>X-axis label</b>"),

                    textInput(
                        inputId = "x_axis_label_input",
                        label = NULL,
                        value = reactives$plot_config$x_label
                    )
                )
            ),

            column(
                width = 12,

                p(HTML("<b>Y-axis label</b>"),

                    textInput(
                        inputId = "y_axis_label_input",
                        label = NULL,
                        value = reactives$plot_config$y_label
                    )
                )
            ),

            column(
                width = 12,

                p(HTML("<b>Number of x-ticks</b>"),

                    sliderInput(
                        inputId = "axis_n_xticks_input",
                        label = NULL,
                        value = reactives$plot_config$n_xticks,
                        min = 3,
                        max = 10,
                        step = 1
                    )
                )
            ),

            column(
                width = 12,

                p(HTML("<b>Number of y-ticks</b>"),

                    sliderInput(
                        inputId = "axis_n_yticks_input",
                        label = NULL,
                        value = reactives$plot_config$n_yticks,
                        min = 3,
                        max = 10,
                        step = 1
                    )
                )
            ),

            column(
                width = 12,

                p(HTML("<b>Tick length (px)</b>"),

                    sliderInput(
                        inputId = "axis_tick_length_input",
                        label = NULL,
                        value = reactives$plot_config$tick_length,
                        min = 4,
                        max = 16,
                        step = 1
                    )
                )
            ),

            column(
                width = 12,

                p(HTML("<b>Tick width (px)</b>"),

                    sliderInput(
                        inputId = "axis_tick_width_input",
                        label = NULL,
                        value = reactives$plot_config$tick_width,
                        min = 0,
                        max = 10,
                        step = 1
                    )
                )
            ),

            column(
                width = 12,

                p(HTML("<b>Show grid</b>"),

                    column(6,

                        checkboxInput(
                            inputId = "x_axis_show_grid_input",
                            label = "x-axis",
                            value = isTRUE(reactives$plot_config$show_grid_x)
                        )

                    ),

                    column(6,

                        checkboxInput(
                            inputId = "y_axis_show_grid_input",
                            label = "y-axis",
                            value = isTRUE(reactives$plot_config$show_grid_y)
                        )

                    )

                )
            )

        ),
        size = 'm',
        easyClose = TRUE,
        footer = tagList(
            modalButton("Close")
        )
    ))
}, ignoreInit = TRUE)

# Observers to update reactives when the modal inputs change
observeEvent(input$axis_n_xticks_input, {
    val <- as.integer(input$axis_n_xticks_input)
    if (is.na(val) || val < 1) return()
    val <- min(10, max(3, val))
    reactives$plot_config$n_xticks <- val
})

observeEvent(input$axis_n_yticks_input, {
    val <- as.integer(input$axis_n_yticks_input)
    if (is.na(val) || val < 1) return()
    val <- min(10, max(3, val))
    reactives$plot_config$n_yticks <- val
})

observeEvent(input$axis_tick_length_input, {
    val <- as.numeric(input$axis_tick_length_input)
    if (is.na(val)) return()
    val <- min(16, max(4, val))
    reactives$plot_config$tick_length <- val
})

observeEvent(input$axis_tick_width_input, {
    val <- as.numeric(input$axis_tick_width_input)
    if (is.na(val)) return()
    val <- min(10, max(0, val))
    reactives$plot_config$tick_width <- val
})

observeEvent(input$y_axis_show_grid_input, {
    reactives$plot_config$show_grid_x <- isTRUE(input$y_axis_show_grid_input)
})

observeEvent(input$x_axis_show_grid_input, {
    reactives$plot_config$show_grid_y <- isTRUE(input$x_axis_show_grid_input)
})

observeEvent(input$x_axis_label_input, {
    reactives$plot_config$x_label <- input$x_axis_label_input
})

observeEvent(input$y_axis_label_input, {
    reactives$plot_config$y_label <- input$y_axis_label_input
})

# Appearance modal: marker size, font size, line width, max points
observeEvent(list(input$configure_appearance), {

    req(reactives$traces_loaded)

    ids <- as.character(hot_to_r(input$legendInfo)$Internal_ID)

    showModal(modalDialog(
        h3("Configure appearance settings"),

        fluidRow(

            column(
                width = 12,

                p(HTML("<b>Marker size</b>"),

                    sliderInput(
                        'appearance_marker_size', NULL,
                        value = reactives$plot_config$marker_size, min = 1, max = 20, step = 0.5
                    )
                )
            ),

            column(
                width = 12,

                p(HTML("<b>Line width</b>"),

                    sliderInput(
                        'appearance_line_width', NULL,
                        value = reactives$plot_config$line_width, min = 0.5, max = 10, step = 0.5
                    )
                )
            ),

            column(
                width = 12,

                p(HTML("<b>Font size</b>"),

                    sliderInput(
                        'appearance_font_size', NULL,
                        value = reactives$plot_config$axis_size, min = 8, max = 34, step = 1
                    )
                )
            ),

            column(
                12,
                column(
                    width = 6,
                    p(HTML('<b>Set colour</b>'),
                        selectInput(inputId="mol2changeColor",
                          label=NULL,
                          choices=ids,
                          selected = reactives$selected_trace,
                          selectize = FALSE
                        )
                    )
                ),

                column(
                    width = 6,
                    p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                        colourpicker::colourInput("colorForLegend", NULL, value = reactives$selected_color)
                    )
                )
            ),

        ),
        easyClose = TRUE,
        footer = tagList(
            modalButton('Close')
        )
    ))
}, ignoreInit = TRUE)

# Observers to sync appearance inputs into reactives
observeEvent(input$appearance_marker_size, {
    val <- as.numeric(input$appearance_marker_size)
    if (is.na(val)) return()
    val <- max(0.1, min(100, val))
    reactives$plot_config$marker_size <- val
})

observeEvent(input$appearance_font_size, {
    val <- as.numeric(input$appearance_font_size)
    if (is.na(val)) return()
    val <- max(8, min(34, round(val)))
    reactives$plot_config$axis_size <- val
})


observeEvent(input$appearance_line_width, {
    val <- as.numeric(input$appearance_line_width)
    if (is.na(val)) return()
    val <- max(0.1, min(100, val))
    reactives$plot_config$line_width <- val
})


