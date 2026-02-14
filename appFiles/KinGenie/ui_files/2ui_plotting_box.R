box(title = "3. Plotting", width = 8, solidHeader = T, status = "primary",

    fluidRow(

        column(
            width=3,
            p(HTML('<p style="margin-bottom:0px;"></p>'),
            actionButton(
                inputId = "configure_appearance_fit",
                label = "Font, markers, lines",
                icon("marker"),
                style="color: #fff; background-color: #337ab7;border-color: #2e6da4")
            )
        ),

        column(
            width=3,
            p(HTML('<p style="margin-bottom:0px;"></p>'),
                actionButton(
                    inputId = "configure_axis_fit",
                    label = "Axis",
                    icon("grip"),
                    style="color: #fff; background-color: #337ab7;border-color: #2e6da4"
                )
            )
        ),

        column(
            width=3,
            p(HTML('<p style="margin-bottom:0px;"></p>'),
               actionButton(
                  inputId = "configure_export_fit",
                  label = "Export",
                    icon("maximize"),
                    style="color: #fff; background-color: #337ab7;border-color: #2e6da4"
               )
            )
        ),

        column(
            width=3,
            p(HTML('<p style="margin-bottom:0px;"></p>'),
               actionButton(
                  inputId = "configure_visualization",
                  label = "Visualization",
                    icon("wrench"),
                    style="color: #fff; background-color: #337ab7;border-color: #2e6da4"
               )
            )
        )

    )

)