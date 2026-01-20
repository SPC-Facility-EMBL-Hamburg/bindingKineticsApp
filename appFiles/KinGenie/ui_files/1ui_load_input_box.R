box(title = "1. Input", width = 12, solidHeader = T, status = "primary",

    fluidRow(

        column(9, p(HTML("<b>Kinetic data files</b>"),
            span(shiny::icon("info-circle"), id = "info_uu1-1"),
            fileInput("kineticFiles", NULL,multiple = TRUE),
            tippy::tippy_this(
            elementId = "info_uu1-1",
            tooltip = "For Octet BLI, load (simultaneously) all the sensor files (.frd)  and the Method file (.fmf).",
            placement = "right"))),

        column(3, p(HTML("<b></b>"),
            withBusyIndicatorUI(shinyjs::hidden(actionButton("Go","",class = "btn-primary")))))

    ),

    fluidRow(

        column(4, p(HTML("<b></b>"),
            actionButton(inputId = "loadExampleData",
            label = "Load example data!",
            icon("caret-right"),
            style="color: #FFFFFF; background-color: #00829c;
            border-color: #00829c")))

    ),

    fluidRow(

        column(8, p(HTML("<b>Select experiment to delete</b>"),
        selectInput("experiment2delete", NULL,c(	'None' = 'None'),selectize = FALSE))),

        column(4, p(HTML("<b><br></b>")),
            actionButton(inputId = "triggerDeletion",label = "",
             icon("trash-can"),
             style="color: #0E090D; background-color: #DABFDF;
             border-color: #6A4D71"))

    )

)