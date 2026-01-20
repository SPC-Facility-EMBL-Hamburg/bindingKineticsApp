box(title = "2. Fitting", width = 12, solidHeader = T, status = "primary",
    fluidRow(

        column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
            actionButton(
            inputId = "triggerCreateDatasetSolution",label = "a. Create dataset",
            style="color: #fff; background-color: #337ab7;
            border-color: #2e6da4"))
        ),

        # Model selection - between single and double exponential
        column(3, p(HTML('b. Model'),
            selectInput(
            inputId = "modelSelectionSolution",
            label = NULL,
            choices = c(
            "Single exponential" = "single",
            "Double exponential" = "double",
            "One binding site"   = "one_binding_site",
            "One binding site (induced fit)" = "one_binding_site_if"
            ),
            selected = "single",
            selectize = FALSE,
        ))),

        column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
            actionButton(
            inputId = "triggerFitSolution",label = "c. Fit",
            style="color: #fff; background-color: #337ab7;
            border-color: #2e6da4"))
        )

    )

)
