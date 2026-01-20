box(title = "3. Binding parameters", width = 12, solidHeader = T, status = "primary",

    conditionalPanel("input.model_selected_sim == 'one_site' ||
                      input.model_selected_sim == 'one_site_mtl'",

        fluidRow(

            column(2, p(HTML("<b>K<sub>d</sub> [μM]</b>"),numericInput('kd_sim_1to1',NULL, 0.5,min = 0, max = 1e8))),
            column(2, p(HTML("<b>k<sub>off</sub> [1/s]</b>"),numericInput('koff_sim_1to1',NULL, 0.01,min = 0, max = 5000))),

            conditionalPanel("input.model_selected_sim == 'one_site_mtl'",

                column(2, p(HTML("<b>K<sub>tr</sub> [1/s]</b>"),numericInput('ktr_sim_1to1',NULL, 0.005,min = 0, max = 1e8)))

            )

        )

    ),

    conditionalPanel("input.model_selected_sim == 'heterogeneous_analyte'",

        fluidRow(

            column(2, p(HTML("<b>K<sub>d,1</sub> [μM]</b>"),numericInput('kd1_sim_hetero',NULL, 0.5,min = 0, max = 1e8))),
            column(2, p(HTML("<b>k<sub>off,1</sub> [1/s]</b>"),numericInput('koff1_sim_hetero',NULL, 0.01,min = 0, max = 5000))),
            column(2, p(HTML("<b>K<sub>d,2</sub> [μM]</b>"),numericInput('kd2_sim_hetero',NULL, 0.5,min = 0, max = 1e8))),
            column(2, p(HTML("<b>k<sub>off,2</sub> [1/s]</b>"),numericInput('koff2_sim_hetero',NULL, 0.01,min = 0, max = 5000)))

        ),

        fluidRow(

            column(2, p(HTML("<b>Population 1 [%]</b>"),numericInput('pop1_sim_hetero',NULL, 50,min = 0, max = 100))),
            column(2, p(HTML("<b>Population 2 [%]</b>"),textOutput('pop2_sim_hetero'))),
            column(2, p(HTML("<b>Population 1 Smax</b>"),numericInput('pop1_sim_smax',NULL, 2.5,min = 0, max = 100))),
            column(2, p(HTML("<b>Population 2 Smax</b>"),textOutput('pop2_sim_smax')))

        )

    ),

    conditionalPanel("input.model_selected_sim == 'one_site_induced_fit' ||
                          input.model_selected_sim == 'one_site_conformational_selection'",

        fluidRow(

            column(2, p(HTML("<b>k<sub>on</sub> [1/(μM*s)]</b>"),numericInput('kon_sim_1to1_adv',NULL, 0.5,min = 0, max = 1e8))),
            column(2, p(HTML("<b>k<sub>off</sub> [1/s]</b>"),numericInput('koff_sim_1to1_adv',NULL, 0.01,  min = 0, max = 5000))),
            column(2, p(HTML("<b>k<sub>c</sub>   [1/s]</b>"),numericInput('kc_sim_1to1_adv',NULL, 1,       min = 0, max = 1e8))),
            column(2, p(HTML("<b>k<sub>rev</sub> [1/s]</b>"),numericInput('krev_sim_1to1_adv',NULL, 10,    min = 0, max = 1e8)))

        )

    ),

    conditionalPanel("input.model_type_sim == 'solution'",

        fluidRow(

            # Remove the hidden statement to allow simulating a signal proportional to the ligand or protein alone concentration

            column(2, p(HTML("<b>Free protein signal</b>"),numericInput('signal_E',NULL, 0,min = 0, max = 1e8))),
            column(2, p(HTML("<b>Free ligand signal</b>"),numericInput('signal_S',NULL, 0, min = 0, max = 1e8))),

            conditionalPanel("input.model_selected_sim == 'one_site' ||
                              input.model_selected_sim == 'one_site_conformational_selection'",

                column(2, p(HTML("<b>Complex signal </b>"),numericInput('signal_ES_simple',NULL, 1, min = 0, max = 1e8)))

            ),

            conditionalPanel("input.model_selected_sim == 'one_site_induced_fit'",

                column(2, p(HTML("<b>Intermediate complex signal </b>"),numericInput('signal_ES_int',NULL, 1,       min = 0, max = 1e8))),
                column(2, p(HTML("<b>Induced complex signal</b>"),numericInput('signal_ES',NULL, 1,    min = 0, max = 1e8)))

            )

        )
    )
)



