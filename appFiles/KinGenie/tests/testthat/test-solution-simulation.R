# AI-generated app integration tests: check solution input forwarding, result data and validation.
box::use(
  ../../box/modules/simulate/model_selection_and_run[modelSelectionAndRunServer]
)

for (solution_model in c("one_site", "one_site_induced_fit", "one_site_conformational_selection")) {
  local({
    # AI-generated: give each test a stable model value without patching dependencies.
    model <- solution_model
    testthat::test_that(paste(model, "forwards simulation inputs and supplies plot data"), {
      # AI-generated: run real panel callbacks with multiple concentrations and unequal species signals.
      state <- shiny::reactiveValues()
      results <- shiny::reactiveValues()
      shiny::testServer(modelSelectionAndRunServer,
        args = list(sim_state = state, sim_results = results, pykingenie = reticulate::import("pykingenie")), {
        session$setInputs(model_type_sim = "solution", model_selected_sim = model, time_step = 0.05)
        if (model == "one_site_induced_fit") {
          testthat::expect_equal(state$kc_sim_1to1_adv, 10)
          testthat::expect_equal(state$krev_sim_1to1_adv, 1)
        }
        state$total_time <- 3
        state$time_step <- 0.05
        state$init_lig_sim_solution <- 2.4
        state$protein_conc_sim <- 0.8
        state$lig_dil_factor_sim <- 2
        state$prot_dil_factor_sim <- 2
        state$numb_dil_sim <- 2
        state$numb_dil_sim_prot <- 1
        state$kd_sim_1to1 <- 0.4
        state$koff_sim_1to1 <- 0.3
        state$kon_sim_1to1_adv <- 1.3
        state$koff_sim_1to1_adv <- 0.3
        state$kc_sim_1to1_adv <- 0.7
        state$krev_sim_1to1_adv <- 1.1
        state$signal_E <- 0.2
        state$signal_S <- 0.4
        state$signal_ES_simple <- 1.7
        state$signal_ES_int <- 1.7
        state$signal_ES <- 2.3
        session$setInputs(btn_cal_simulation = 1)
        testthat::expect_true(results$available)
        data <- results$data
        testthat::expect_equal(data$prot_concs, c(0.8, 0.4))
        testthat::expect_equal(data$lig_concs, c(0.6, 1.2, 2.4))
        for (i in seq_along(data$prot_concs)) {
          for (j in seq_along(data$lig_concs)) {
            pc <- data$prot_concs[i]
            lc <- data$lig_concs[j]
            # Compare the app's forwarding with the package API, without retesting its equations.
            py <- reticulate::import("pykingenie")
            times <- reticulate::np_array(data$interaction_time)
            if (model == "one_site") {
              expected <- py$signal_ode_one_site_insolution(t = times,
                koff = 0.3, Kd = 0.4, a_total = pc, b_total = lc,
                signal_a = 0.2, signal_b = 0.4, signal_complex = 1.7)
            } else if (model == "one_site_induced_fit") {
              expected <- py$signal_ode_induced_fit_insolution(t = times, y = list(0, 0),
                k1 = 1.3, k_minus1 = 0.3, k2 = 0.7, k_minus2 = 1.1,
                E_tot = pc, S_tot = lc,
                signal_E = 0.2, signal_S = 0.4, signal_ES_int = 1.7, signal_ES = 2.3)
            } else {
              initial <- py$get_initial_concentration_conformational_selection(
                total_protein = pc, k1 = 0.7, k_minus1 = 1.1)
              expected <- py$signal_ode_conformational_selection_insolution(t = times,
                y = list(initial[2], 0), k1 = 0.7, k_minus1 = 1.1,
                k2 = 1.3, k_minus2 = 0.3, E_tot = pc, S_tot = lc,
                signal_E1 = 0.2, signal_E2 = 0.2, signal_S = 0.4, signal_E2S = 1.7)
            }
            testthat::expect_equal(as.numeric(data$signal_per_pc[[i]][[j]]), as.numeric(expected))
            if (model != "one_site") {
              rate_function <- if (model == "one_site_induced_fit") {
                py$get_kobs_induced_fit
              } else py$get_kobs_conformational_selection
              testthat::expect_equal(data$k_obs_dominant_per_pc[[i]][j],
                rate_function(tot_lig = lc, tot_prot = pc, kr = 1.1, ke = 0.7,
                  kon = 1.3, koff = 0.3, dominant = TRUE))
              testthat::expect_equal(data$k_obs_non_dominant_per_pc[[i]][j],
                rate_function(tot_lig = lc, tot_prot = pc, kr = 1.1, ke = 0.7,
                  kon = 1.3, koff = 0.3, dominant = FALSE))
            }
          }
        }
        # Reject invalid grids before calling Python, and allow the next valid run.
        for (field in c("time_step", "lig_dil_factor_sim", "numb_dil_sim_prot", "kd_sim_1to1")) {
          if (field == "kd_sim_1to1" && model != "one_site") next
          previous <- state[[field]]
          state[[field]] <- if (field == "numb_dil_sim_prot") 0.5 else 0
          session$setInputs(btn_cal_simulation = input$btn_cal_simulation + 1)
          testthat::expect_false(results$available)
          testthat::expect_null(results$data)
          state[[field]] <- previous
        }
        session$setInputs(btn_cal_simulation = input$btn_cal_simulation + 1)
        testthat::expect_true(results$available)
        if (model == "one_site_induced_fit") {
          state$signal_ES <- state$signal_ES_int
          session$setInputs(btn_cal_simulation = input$btn_cal_simulation + 1)
          testthat::expect_true(results$relaxation_plot_done)
        }
        state$init_lig_sim_solution <- 0
        session$setInputs(btn_cal_simulation = input$btn_cal_simulation + 1)
        testthat::expect_true(results$available)
        testthat::expect_true(all(is.finite(unlist(results$data$signal_per_pc))))
      })
    })
  })
}
