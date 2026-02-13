box(title = "3. Plotting", width = 12, solidHeader = T, status = "primary",

    column(6,

        column(
          width=5,
          p(HTML('<p style="margin-bottom:0px;"></p>'),
            actionButton(
              inputId = "configure_appearance",
              label = "Font, markers, colors",
              icon("marker"),
              style="color: #fff; background-color: #337ab7;
              border-color: #2e6da4")
          )
        ),

        column(
          width=3,
          p(HTML('<p style="margin-bottom:0px;"></p>'),
            actionButton(
              inputId = "configure_axis",
              label = "Axis",
              icon("grip"),
              style="color: #fff; background-color: #337ab7;
              border-color: #2e6da4")
          )
        ),

        column(
          width=3,
          p(HTML('<p style="margin-bottom:0px;"></p>'),
            actionButton(
              inputId = "configure_export",
              label = "Export",
              icon("maximize"),
              style="color: #fff; background-color: #337ab7;
              border-color: #2e6da4")
          )
        )

    ),

    column(6,

        fluidRow(column(width = 12,rHandsontableOutput('legendInfo')))

    )

)