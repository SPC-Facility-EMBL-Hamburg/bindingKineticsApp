box::use(
  shiny[
    checkboxInput,
    column,
    fluidRow,
    HTML,
    moduleServer,
    NS,
    numericInput,
    observeEvent,
    p,
    renderUI,
    req,
    tagList,
    uiOutput
  ],
  shinydashboard[
    box
  ]
)
#' @export
experimentalParamsUI <- function(id) {
  ns <- NS(id)

  tagList(
    box(
      title = "2. Experimental parameters",
      width = 9,
      solidHeader = TRUE,
      status = "primary",
      uiOutput(ns("experimental_params"))
    )
  )
}

#' @export
experimentalParamsServer <- function(id, sim_state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$experimental_params <- renderUI({
      req(sim_state$model_type_sim)

      model_type <- sim_state$model_type_sim
      model_selected <- sim_state$model_selected_sim

      tagList(
        ## ---- First row -------------------------------------------------------
        fluidRow(
          if (model_type == "surface") {
            column(
              4,
              p(
                HTML("<b>[Analyte] (μM)</b>"),
                numericInput(ns("init_lig_sim_surface"), NULL, 0.1,
                  min = 0, max = 1e6
                )
              )
            )
          } else {
            column(
              4,
              p(
                HTML("<b>[Ligand] (μM)</b>"),
                numericInput(ns("init_lig_sim_solution"), NULL, 0.1,
                  min = 0, max = 1e6
                )
              )
            )
          },
          column(
            4,
            p(
              HTML("<b># Dilutions</b>"),
              numericInput(ns("numb_dil_sim"), NULL, 7,
                min = 0, max = 1e6
              )
            )
          ),
          column(
            4,
            p(
              HTML("<b>Dilution factor</b>"),
              numericInput(ns("lig_dil_factor_sim"), NULL, 2.2,
                min = 0, max = 1e6
              )
            )
          )
        ),

        ## ---- Second row ------------------------------------------------------

        fluidRow(
          if (model_type == "surface" &&
            !model_selected %in% c(
              "heterogeneous_analyte",
              "ligand_has_two_sites"
            )) {
            column(
              4,
              p(
                HTML("<b>Rmax (A.U.)</b>"),
                numericInput(ns("protein_smax_sim"), NULL, 2.5,
                  min = 0, max = 1e6
                )
              )
            )
          },
          if (model_type == "surface" &&
            model_selected == "ligand_has_two_sites") {
            tagList(
              column(
                4,
                p(
                  HTML("<b>PL Rmax (A.U.)</b>"),
                  numericInput(ns("pl_rmax_sim_two_sites"), NULL, 3.5,
                    min = 0, max = 1e6
                  )
                )
              ),
              column(
                4,
                p(
                  HTML("<b>PL<sub>2</sub> Rmax (A.U.)</b>"),
                  numericInput(ns("pl2_rmax_sim_two_sites"), NULL, 5,
                    min = 0, max = 1e6
                  )
                )
              )
            )
          },
          if (model_type == "surface" &&
            model_selected == "heterogeneous_analyte") {
            column(
              4,
              p(
                HTML("<b>Total Smax (A.U.)</b>"),
                numericInput(ns("total_smax_sim"), NULL, 5,
                  min = 0, max = 1e6
                )
              )
            )
          },
          if (model_type != "surface") {
            column(
              4,
              p(
                HTML("<b>[Protein] (μM)</b>"),
                numericInput(ns("protein_conc_sim"), NULL, 0.1,
                  min = 0, max = 1e6
                )
              )
            )
          },
          column(
            4,
            p(
              HTML("<b># Dilutions</b>"),
              numericInput(ns("numb_dil_sim_prot"), NULL, 0)
            )
          ),
          column(
            4,
            p(
              HTML("<b>Dilution factor</b>"),
              numericInput(ns("prot_dil_factor_sim"), NULL, 2)
            )
          )
        ),

        ## ---- Third row -------------------------------------------------------

        fluidRow(
          if (model_type == "surface") {
            tagList(
              column(
                4,
                p(
                  HTML("<b>Association time (sec)</b>"),
                  numericInput(ns("association_time"), NULL, 300,
                    min = 0, max = 5000
                  )
                )
              ),
              column(
                4,
                p(
                  HTML("<b>Dissociation time (sec)</b>"),
                  numericInput(ns("dissociation_time"), NULL, 600,
                    min = 0, max = 5000
                  )
                )
              ),
              column(
                4,
                p(
                  HTML("<b>Single-cycle</b>"),
                  checkboxInput(ns("is_single_cycle_sim"), NULL, FALSE)
                )
              )
            )
          } else {
            column(
              4,
              p(
                HTML("<b>Total time (sec)</b>"),
                numericInput(ns("total_time"), NULL, 300,
                  min = 0, max = 5000
                )
              )
            )
          }
        )
      )
    })

    params <- c(
      "init_lig_sim_surface",
      "init_lig_sim_solution",
      "protein_conc_sim",
      "protein_smax_sim",
      "pl_rmax_sim_two_sites",
      "pl2_rmax_sim_two_sites",
      "total_smax_sim",
      "numb_dil_sim",
      "lig_dil_factor_sim",
      "numb_dil_sim_prot",
      "prot_dil_factor_sim",
      "association_time",
      "dissociation_time",
      "is_single_cycle_sim",
      "total_time"
    )

    lapply(params, function(param) {
      observeEvent(input[[param]],
        {
          sim_state[[param]] <- input[[param]]
        },
        ignoreNULL = FALSE
      )
    })
  })
}
