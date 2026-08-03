box::use(
  shiny[
    column,
    downloadButton,
    downloadHandler,
    HTML,
    moduleServer,
    NS,
    p,
    req,
    tagList
  ],
  utils[
    write.csv
  ],
)

#' @export
simResultsExportUI <- function(id) {
  ns <- NS(id)

  tagList(
    column(7, p(
      HTML('<p style="margin-bottom:0px;"><br></p>'),
      downloadButton(
        ns("btn_export_simulation"),
        label = "Export sim. results",
        style = "color: #fff; background-color: #6c757d; border-color: #545b62;"
      )
    ))
  )
}

#' @export
simResultsExportServer <- function(id, sim_results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns


    output$btn_export_simulation <- downloadHandler(
      filename = function() {

          if (getOption("shiny.testmode", FALSE)) {
            return("simulation.csv")
          }

          paste0(
            "simulation_KinGenie_",
            Sys.Date(),
            ".csv"
          )

        },

      content = function(file) {
        req(sim_results$available)

        sig_digits <- 6

        sim_model <- sim_results[["model_type"]]

        dfs <- list()

        if (sim_model == "solution") {
          interaction_time <- sim_results$data[["interaction_time"]]
          signal_per_pc <- sim_results$data[["signal_per_pc"]]
          prot_concs <- sim_results$data[["prot_concs"]]
          lig_concs <- sim_results$data[["lig_concs"]]

          for (i in 1:length(prot_concs)) {
            for (j in 1:length(lig_concs)) {
              df <- data.frame(
                "Time" = signif(interaction_time, sig_digits),
                "Signal" = signif(signal_per_pc[[i]][[j]], sig_digits),
                "Protein_concentration_micromolar" = prot_concs[i],
                "Ligand_concentration_micromolar" = lig_concs[j]
              )

              dfs[[length(dfs) + 1]] <- df
            }
          }
        } else {
          is_single_cycle <- sim_results[["is_single_cycle"]]

          association_time_per_cycle <- sim_results$data[["association_time_per_cycle"]]
          dissociation_time_per_cycle <- sim_results$data[["dissociation_time_per_cycle"]]
          signal_per_cycle_assoc <- sim_results$data[["signal_per_cycle_assoc"]]
          signal_per_cycle_disso <- sim_results$data[["signal_per_cycle_disso"]]
          smax_all <- sim_results$data[["smax_all"]]
          lig_concs <- sim_results$data[["lig_concs"]]

          # Option 1 - multi-cycle kinetics, all ligand concentration measured in parallel
          if (!is_single_cycle) {
            for (i in 1:length(smax_all)) {
              for (j in 1:length(lig_concs)) {
                df1 <- data.frame(
                  "Time" = signif(association_time_per_cycle[[1]], sig_digits),
                  "Signal" = signif(signal_per_cycle_assoc[[1]][[i]][[j]], sig_digits),
                  "Smax" = smax_all[[i]],
                  "Analyte_concentration_micromolar_constant" = lig_concs[j]
                )

                df2 <- data.frame(
                  "Time" = signif(dissociation_time_per_cycle[[1]], sig_digits),
                  "Signal" = signif(signal_per_cycle_disso[[1]][[i]][[j]], sig_digits),
                  "Smax" = smax_all[[i]],
                  "Analyte_concentration_micromolar_constant" = 0
                )

                dfs[[length(dfs) + 1]] <- df1
                dfs[[length(dfs) + 1]] <- df2
              }
            }

            # Option 2 - single-cycle kinetics, all ligand concentration measured one after another
          } else {
            for (cycle in 1:length(association_time_per_cycle)) {
              for (i in 1:length(smax_all)) {
                df1 <- data.frame(
                  "Time" = signif(association_time_per_cycle[[cycle]], sig_digits),
                  "Signal" = signif(signal_per_cycle_assoc[[cycle]][[i]][[1]], sig_digits),
                  "Smax" = smax_all[[i]],
                  "Analyte_concentration_micromolar_constant" = lig_concs[cycle],
                  "Cycle" = cycle
                )

                df2 <- data.frame(
                  "Time" = signif(dissociation_time_per_cycle[[cycle]], sig_digits),
                  "Signal" = signif(signal_per_cycle_disso[[cycle]][[i]][[1]], sig_digits),
                  "Smax" = smax_all[[i]],
                  "Analyte_concentration_micromolar_constant" = 0,
                  "Cycle" = cycle
                )

                dfs[[length(dfs) + 1]] <- df1
                dfs[[length(dfs) + 1]] <- df2
              }
            }
          }
        }

        df <- do.call(rbind, dfs)
        write.csv(df, file = file, row.names = FALSE)
      }
    )
  })
}
