box(title = "1. Dataset info", width = 12, solidHeader = T, status = "primary",
    fluidRow(

        # Add checkbox to show or hide the ligandInfo dataframe
        column(4,
            checkboxInput("showLigandInfo", "Show dataset info", value = TRUE)
        ),

        # Add button for quick selection of the ligand info
        column(4,p(HTML("<br>"),
            actionButton("quickSelect", "Quick select")
        )),

        conditionalPanel('input.showLigandInfo',

            column(12, rHandsontableOutput("ligandInfo"))

        )
    )
)
