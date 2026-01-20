box(title = "2. Processing", width = 12, solidHeader = T, status = "primary",
    fluidRow(

      column(6, p(HTML("<b>Experiment name</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-selectedExperiment"),
                  selectInput("selectedExperiment", NULL,
                              c('None'       = 'none'),selectize = FALSE),
                  tippy::tippy_this(
                    elementId = "info_uu-selectedExperiment",
                    tooltip = "Select the experiment to process.
                    ",placement = "right"))),

      column(6, p(HTML("<b>Operation</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-inputOperation"),
                  selectInput("operation", NULL,
                              c('Subtract baseline'       = 'subtract',
                                'Align association phase' = 'align_association',
                                'Inter-step correction (dissociation)'= 'correct_dissociation',
                                'Average'                 = 'average',
                                'Merge steps'             = 'merge_steps'),
                                selectize = FALSE),
                                #'Smooth'                  = 'smooth')),
                  tippy::tippy_this(
                    elementId = "info_uu-inputOperation",
                    tooltip = "Select and press 'Apply'.
                    ",placement = "right")))
    ),

    fluidRow(

      column(4, p(HTML("<b><br></b>"),
                  actionButton(
                    inputId = "triggerProcessing",label = "Apply",
                    icon("forward"),
                    style="color: #fff; background-color: #337ab7;
               border-color: #2e6da4")))

))
