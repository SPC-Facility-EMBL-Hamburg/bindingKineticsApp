box(title = "2. Experimental parameters", width = 9, solidHeader = T, status = "primary",

    fluidRow(

        conditionalPanel(
            condition = "input.model_type_sim == 'surface'",
            column(4, p(HTML("<b>[Analyte]  (μM)</b>"),
            numericInput('init_lig_sim',NULL, 0.1,min = 0, max = 1e6)))
        ),

        conditionalPanel(
            condition = "input.model_type_sim != 'surface'",
            column(4, p(HTML("<b>[Ligand]  (μM)</b>"),
            numericInput('init_lig_sim',NULL, 0.1,min = 0, max = 1e6)))
        ),


        column(4, p(HTML("<b># Dilutions </b>"),
            numericInput('numb_dil_sim',NULL, 7,min = 0, max = 1e6))),

        column(4, p(HTML("<b>Dilution factor </b>"),
            numericInput('lig_dil_factor_sim',NULL, 2.2,min = 0, max = 1e6)))

    ),

    fluidRow(

        conditionalPanel(
            condition = "input.model_type_sim == 'surface' &&
            input.model_selected_sim != 'heterogeneous_analyte'",
            column(4, p(HTML("<b>Rmax (A.U.)</b>"),
            numericInput('protein_smax_sim',NULL, 5,min = 0, max = 1e6))),
        ),

        conditionalPanel(
            condition = "input.model_type_sim == 'surface' &&
            input.model_selected_sim == 'heterogeneous_analyte'",
            column(4, p(HTML("<b>Total Smax (A.U.)</b>"),
            numericInput('total_smax_sim',NULL, 5,min = 0, max = 1e6))),
        ),

        conditionalPanel(
            condition = "input.model_type_sim != 'surface'",
            column(4, p(HTML("<b>[Protein] (μM)</b>"),
            numericInput('protein_conc_sim',NULL, 0.1,min = 0, max = 1e6))),
        ),

        column(4, p(HTML("<b># Dilutions </b>"),
            numericInput('numb_dil_sim_prot',NULL, 0,min = 0, max = 1e6))),

        column(4, p(HTML("<b>Dilution factor </b>"),
            numericInput('prot_dil_factor_sim',NULL, 2,min = 0, max = 1e6)))

    ),

    fluidRow(

        conditionalPanel(

            condition = "input.model_type_sim == 'surface'",

            column(4, p(HTML("<b>Association time (sec)</b>"), numericInput('association_time',NULL, 300,min = 0, max = 5000))),
            column(4, p(HTML("<b>Dissociation time (sec)</b>"),numericInput('dissociation_time',NULL, 600,min = 0, max = 5000))),
            column(4, p(HTML("<b>Single-cycle</b>"),           checkboxInput('is_single_cycle_sim', NULL, value = FALSE)))
            #column(4, p(HTML("<b>Number of cycles</b>"),       numericInput('numb_cycles_sim',NULL, 1,min = 1, max = 10)))

        ),

        conditionalPanel(
            condition = "input.model_type_sim != 'surface'",

            column(4, p(HTML("<b>Total time (sec)</b>"), numericInput('total_time',NULL, 300,min = 0, max = 5000))),

        )


    )

)