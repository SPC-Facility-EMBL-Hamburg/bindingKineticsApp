box::use(
  shiny[
    column,
    fluidRow,
    HTML,
    moduleServer,
    NS,
    numericInput,
    observeEvent,
    p,
    renderText,
    renderUI,
    req,
    tagList,
    textOutput,
    uiOutput,
    updateNumericInput
  ],
  shinydashboard[
    box
  ]
)
#' @export
bindingParamsUI <- function(id) {
  ns <- NS(id)

  tagList(
    box(
      title = "3. Binding parameters",
      width = 12,
      solidHeader = TRUE,
      status = "primary",
      uiOutput(ns("binding_params"))
    )
  )
}

#' @export
bindingParamsServer <- function(id, sim_state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$binding_params <- renderUI({
      req(
        sim_state$model_selected_sim,
        sim_state$model_type_sim
      )

      model <- sim_state$model_selected_sim
      type <- sim_state$model_type_sim

      tagList(
        ## --------------------------------------------------
        ## One-site / MTL
        ## --------------------------------------------------
        if (model %in% c("one_site", "one_site_mtl")) {
          fluidRow(
            column(
              2,
              p(
                HTML("<b>K<sub>d</sub> [μM]</b>"),
                numericInput(ns("kd_sim_1to1"), NULL, 0.5, min = 0, max = 1e8)
              )
            ),
            column(
              2,
              p(
                HTML("<b>k<sub>off</sub> [1/s]</b>"),
                numericInput(ns("koff_sim_1to1"), NULL, 0.01, min = 0, max = 5000)
              )
            ),
            if (model == "one_site_mtl") {
              column(
                2,
                p(
                  HTML("<b>K<sub>tr</sub> [1/s]</b>"),
                  numericInput(ns("ktr_sim_1to1"), NULL, 0.005, min = 0, max = 1e8)
                )
              )
            }
          )
        },

        ## --------------------------------------------------
        ## Ligand has two sites
        ## --------------------------------------------------

        if (model == "ligand_has_two_sites") {
          fluidRow(
            column(
              2,
              p(
                HTML("<b>k<sub>on</sub> [1/(μM*s)]</b>"),
                numericInput(ns("kon_sim_ligand_two_sites"), NULL, 0.5)
              )
            ),
            column(
              2,
              p(
                HTML("<b>k<sub>off</sub> [1/s]</b>"),
                numericInput(ns("koff_sim_ligand_two_sites"), NULL, 0.01)
              )
            ),
            column(
              2,
              p(
                HTML("<b>Cooperativity factor</b>"),
                numericInput(ns("coop_sim_ligand_two_sites"), NULL, 1)
              )
            )
          )
        },

        ## --------------------------------------------------
        ## Heterogeneous analyte 
        ## --------------------------------------------------

        if (model == "heterogeneous_analyte") {
          tagList(
            fluidRow(
              column(2, p(
                HTML("<b>K<sub>d,1</sub> [μM]</b>"),
                numericInput(ns("kd1_sim_hetero"), NULL, 0.5)
              )),
              column(2, p(
                HTML("<b>k<sub>off,1</sub> [1/s]</b>"),
                numericInput(ns("koff1_sim_hetero"), NULL, 0.01)
              )),
              column(2, p(
                HTML("<b>K<sub>d,2</sub> [μM]</b>"),
                numericInput(ns("kd2_sim_hetero"), NULL, 0.5)
              )),
              column(2, p(
                HTML("<b>k<sub>off,2</sub> [1/s]</b>"),
                numericInput(ns("koff2_sim_hetero"), NULL, 0.01)
              ))
            ),
            fluidRow(
              column(2, p(
                HTML("<b>Population 1 [%]</b>"),
                numericInput(ns("pop1_sim_hetero"), NULL, 50)
              )),
              column(2, p(
                HTML("<b>Population 2 [%]</b>"),
                textOutput(ns("pop2_sim_hetero"))
              )),
              column(2, p(
                HTML("<b>Population 1 Smax</b>"),
                numericInput(ns("pop1_sim_smax"), NULL, 2.5)
              )),
              column(2, p(
                HTML("<b>Population 2 Smax</b>"),
                textOutput(ns("pop2_sim_smax"))
              ))
            )
          )
        },

        ## --------------------------------------------------
        ## Heterogeneous ligand 
        ## --------------------------------------------------

        if (model == "heterogeneous_ligand") {
          tagList(
            fluidRow(
              column(2, p(
                HTML("<b>K<sub>d,1</sub> [μM]</b>"),
                numericInput(ns("kd1_sim_hetero"), NULL, 0.5)
              )),
              column(2, p(
                HTML("<b>k<sub>off,1</sub> [1/s]</b>"),
                numericInput(ns("koff1_sim_hetero"), NULL, 0.01)
              )),
              column(2, p(
                HTML("<b>K<sub>d,2</sub> [μM]</b>"),
                numericInput(ns("kd2_sim_hetero"), NULL, 0.5)
              )),
              column(2, p(
                HTML("<b>k<sub>off,2</sub> [1/s]</b>"),
                numericInput(ns("koff2_sim_hetero"), NULL, 0.01)
              ))
            ),
            fluidRow(
              column(2, p(
                HTML("<b>Fraction site 1 [%]</b>"),
                numericInput(ns("pop1_sim_hetero"), NULL, 50)
              )),
              column(2, p(
                  HTML("<b>Fraction site 2 [%]</b>"),
                  textOutput(ns("pop2_sim_hetero"))
              ))
            )
          )
        },

        ## --------------------------------------------------
        ## Induced fit / Conformational selection
        ## --------------------------------------------------

        if (model %in% c(
          "one_site_induced_fit",
          "one_site_conformational_selection"
        )) {
          fluidRow(
            column(2, p(
              HTML("<b>k<sub>on</sub> [1/(μM*s)]</b>"),
              numericInput(ns("kon_sim_1to1_adv"), NULL, 0.5)
            )),
            column(2, p(
              HTML("<b>k<sub>off</sub> [1/s]</b>"),
              numericInput(ns("koff_sim_1to1_adv"), NULL, 0.01)
            )),
            column(2, p(
              HTML("<b>k<sub>c</sub> [1/s]</b>"),
              numericInput(ns("kc_sim_1to1_adv"), NULL, 1)
            )),
            column(2, p(
              HTML("<b>k<sub>rev</sub> [1/s]</b>"),
              numericInput(ns("krev_sim_1to1_adv"), NULL, 10)
            ))
          )
        },

        ## --------------------------------------------------
        ## Solution-only parameters
        ## --------------------------------------------------

        if (type == "solution") {
          fluidRow(
            column(2, p(
              HTML("<b>Free protein signal</b>"),
              numericInput(ns("signal_E"), NULL, 0)
            )),
            column(2, p(
              HTML("<b>Free ligand signal</b>"),
              numericInput(ns("signal_S"), NULL, 0)
            )),
            if (model %in% c(
              "one_site",
              "one_site_conformational_selection"
            )) {
              column(
                2,
                p(
                  HTML("<b>Complex signal</b>"),
                  numericInput(ns("signal_ES_simple"), NULL, 1)
                )
              )
            },
            if (model == "one_site_induced_fit") {
              tagList(
                column(
                  2,
                  p(
                    HTML("<b>Intermediate complex signal</b>"),
                    numericInput(ns("signal_ES_int"), NULL, 1)
                  )
                ),
                column(
                  2,
                  p(
                    HTML("<b>Induced complex signal</b>"),
                    numericInput(ns("signal_ES"), NULL, 1)
                  )
                )
              )
            }
          )
        }
      )
    })

    params <- c(
      "kd_sim_1to1",
      "koff_sim_1to1",
      "ktr_sim_1to1",
      "kon_sim_ligand_two_sites",
      "koff_sim_ligand_two_sites",
      "coop_sim_ligand_two_sites",
      "kd1_sim_hetero",
      "koff1_sim_hetero",
      "kd2_sim_hetero",
      "koff2_sim_hetero",
      "pop1_sim_hetero",
      "pop1_sim_smax",
      "kon_sim_1to1_adv",
      "koff_sim_1to1_adv",
      "kc_sim_1to1_adv",
      "krev_sim_1to1_adv",
      "signal_E",
      "signal_S",
      "signal_ES_simple",
      "signal_ES_int",
      "signal_ES"
    )

    lapply(params, function(x) {
      observeEvent(input[[x]],
        {
          sim_state[[x]] <- input[[x]]
        },
        ignoreNULL = FALSE
      )
    })

    output$pop2_sim_hetero <- renderText({
      paste0(100 - input$pop1_sim_hetero)
    })

    output$pop2_sim_smax <- renderText({
      paste0(sim_state$total_smax_sim - input$pop1_sim_smax)
    })

    # Update smax of each population to be between 0 and the global smax
    observeEvent(sim_state$total_smax_sim, {
      smax_global <- sim_state$total_smax_sim
      req(input$pop1_sim_smax)
      smax_1 <- input$pop1_sim_smax
      smax_2 <- smax_global - smax_1
      diff <- smax_global - smax_1 - smax_2

      updateNumericInput(session, "pop1_sim_smax", value = smax_1 + diff, max = smax_global)
    })
  })
}
