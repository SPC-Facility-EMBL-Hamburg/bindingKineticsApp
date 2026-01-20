box(title = "2. Fitting", width = 12, solidHeader = T, status = "primary",
    fluidRow(

        column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
            actionButton(
            inputId = "triggerCreateDataset",label = "a. Create dataset",
            style="color: #fff; background-color: #337ab7;
            border-color: #2e6da4"))
        ),

        column(4, p(HTML("<b>b. Region</b>"),
            span(shiny::icon("info-circle"), id = "info_uu-fitRegion"),
            selectInput("fittingRegion", NULL,
                c('Association and dissociation' = 'association_dissociation',
                'Steady-state'                 = 'steady_state',
                'Association'                  = 'association',
                'Dissociation'                 = 'dissociation')),
            tippy::tippy_this(
                elementId = "info_uu-fitRegion",
                tooltip = "Select association and/or dissociation to fit K_d and k_off.
                Select Steady-state to fit K_d.'.
                ",placement = "right"))
        ),

        column(4, p(HTML("<b>c. Model</b>"),
            span(shiny::icon("info-circle"), id = "info_uu-fitModel"),
            selectInput("fittingModel", NULL,
                c('One-to-one'       = 'one_to_one',
                'One-to-one (MTL)' = 'one_to_one_mtl'
                )
            ),
            tippy::tippy_this(
                elementId = "info_uu-fitModel",
                tooltip = "Select and press 'Fit!'.
                ",placement = "right"))
        )
    ),

    fluidRow(

        column(2, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
            actionButton(
            inputId = "triggerFitting",label = "d. Fit!",
            icon("meteor"),
            style="color: #fff; background-color: #337ab7;
            border-color: #2e6da4"))
        ),

        # Button to calculate the asymmetric error
        column(2, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
            actionButton(
            inputId = "triggerAsymError",label = "e. Confidence interval",
            icon("calculator"),
            style="color: #fff; background-color: #337ab7;
            border-color: #2e6da4")))


))
