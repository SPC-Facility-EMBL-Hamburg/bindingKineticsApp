box(title = "1. Model", width = 3, solidHeader = T, status = "primary",

    fluidRow(

        column(6, p(HTML("<b>Type</b>"),
            span(shiny::icon("info-circle"), id = "info_uu_sim_2-0"),
            selectInput("model_type_sim", NULL,
            c( "Surface binding"  = "surface",
                "In-solution"     = "solution")),
            tippy::tippy_this(elementId = "info_uu_sim_2-0",
            tooltip = "Choose if the target protein is bound to a surface (e.g., BLI) or in solution (e.g., NMR).
            In the case of surface binding, the ligand concentration remains constant.
            In the case of in-solution binding, ligand depletion takes place.",
            placement = "right"))),

        column(6, p(HTML("<b>Model</b>"),
            span(shiny::icon("info-circle"), id = "info_uu_sim_2-1"),

                selectInput("model_selected_sim", NULL,
                    c( "1:1"  = "one_site",
                    "1:1 (induced fit)"  = "one_site_induced_fit",
                    "1:1 (conformational selection)"  = "one_site_conformational_selection"#,
                    #"1:1 (two ligands)"  = "heteregeneous_ligand"
                    ), selectize = FALSE), # Set selectize to FALSE to avoid text overflowing outside the box

            tippy::tippy_this(elementId = "info_uu_sim_2-1",
            tooltip = "Select the model for the simulation.
            More information in the User Guide.",placement = "right")))
    ),

    fluidRow(

      column(7, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
        actionButton(
        inputId = "btn_cal_simulation",label = "Run simulation",
        icon("meteor"),
        style="color: #fff; background-color: #337ab7;
        border-color: #2e6da4"))),

      column(5, p(HTML("<b>Time step (s)</b>"),
        span(shiny::icon("info-circle"), id = "info_uu_sim_2-2"),
        numericInput('time_step',NULL, 0.5,min = 0, max = 10),
        tippy::tippy_this(elementId = "info_uu_sim_2-2",
        tooltip = "Time step to simulate the association and dissociation curves.",
        placement = "right"))),

      column(7, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
        downloadButton(
        "btn_export_simulation",
        label = "Export sim. results",
        style = "color: #fff; background-color: #6c757d; border-color: #545b62;"
        )))

    )

)
