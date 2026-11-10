# AI-generated integration tests for the mass-transport simulation API.
box::use(
  ../../box/modules/simulate/model_selection_and_run[modelSelectionAndRunServer]
)

check_mass_transport_simulation <- function(single_cycle) {
  # AI-generated: run the real Shiny module with the installed Python solvers.
  pykingenie <- reticulate::import("pykingenie")
  sim_state <- shiny::reactiveValues()
  sim_results <- shiny::reactiveValues()

  shiny::testServer(
    modelSelectionAndRunServer,
    args = list(sim_state = sim_state, sim_results = sim_results, pykingenie = pykingenie),
    {
      # AI-generated: exercise initial and repeated-cycle calls without replacing dependencies.
      session$setInputs(
        model_type_sim = "surface",
        model_selected_sim = "one_site_mtl",
        time_step = 0.5
      )
      sim_state$is_single_cycle_sim <- single_cycle
      sim_state$init_lig_sim_surface <- 4
      sim_state$lig_dil_factor_sim <- 2
      sim_state$numb_dil_sim <- 2
      sim_state$protein_smax_sim <- 5
      sim_state$prot_dil_factor_sim <- 2
      sim_state$numb_dil_sim_prot <- 1
      sim_state$association_time <- 20
      sim_state$dissociation_time <- 30
      sim_state$kd_sim_1to1 <- 0.5
      sim_state$koff_sim_1to1 <- 0.01
      sim_state$ktr_sim_1to1 <- 0.25

      session$setInputs(btn_cal_simulation = 1)
      testthat::expect_true(sim_results$available)
      data <- sim_results$data
      n_cycles <- if (single_cycle) 3L else 1L
      n_traces <- if (single_cycle) 1L else 3L
      testthat::expect_length(data$signal_per_cycle_assoc, n_cycles)
      testthat::expect_equal(data$smax_all, c(5, 2.5))
      testthat::expect_equal(data$lig_concs, c(1, 2, 4))

      for (cycle in seq_len(n_cycles)) {
        ta <- data$association_time_per_cycle[[cycle]]
        td <- data$dissociation_time_per_cycle[[cycle]]
        for (capacity in seq_along(data$smax_all)) {
          assoc <- data$signal_per_cycle_assoc[[cycle]][[capacity]]
          disso <- data$signal_per_cycle_disso[[cycle]][[capacity]]
          testthat::expect_length(assoc, n_traces)
          for (trace in seq_len(n_traces)) {
            concentration <- data$lig_concs[if (single_cycle) cycle else trace]
            initial <- if (cycle == 1L) 0 else as.numeric(tail(
              data$signal_per_cycle_disso[[cycle - 1L]][[capacity]][[trace]], 1
            ))
            expected_assoc <- as.numeric(pykingenie$solve_ode_one_site_mass_transport_association(
              t = reticulate::np_array(ta), s1_0 = initial,
              analyte_conc = concentration, K_d = 0.5, k_off = 0.01,
              k_tr = 0.25, s_max = data$smax_all[capacity]
            ))
            expected_disso <- as.numeric(pykingenie$solve_ode_one_site_mass_transport_dissociation(
              t = reticulate::np_array(td), s1_0 = tail(expected_assoc, 1),
              K_d = 0.5, k_off = 0.01, k_tr = 0.25, s_max = data$smax_all[capacity]
            ))
            testthat::expect_length(assoc[[trace]], length(ta))
            testthat::expect_length(disso[[trace]], length(td))
            testthat::expect_equal(as.numeric(assoc[[trace]]), expected_assoc, tolerance = 1e-8)
            testthat::expect_equal(as.numeric(disso[[trace]]), expected_disso, tolerance = 1e-8)
            testthat::expect_equal(assoc[[trace]][1], initial, tolerance = 1e-10)
            testthat::expect_equal(disso[[trace]][1], as.numeric(tail(assoc[[trace]], 1)), tolerance = 1e-10)
            testthat::expect_true(all(is.finite(c(assoc[[trace]], disso[[trace]]))))
          }
        }
      }
    }
  )
}

testthat::test_that("mass-transport simulation supports multiple concentrations and capacities", {
  # AI-generated regression for the multi-cycle simulation panel.
  check_mass_transport_simulation(FALSE)
})

testthat::test_that("mass-transport single-cycle simulation carries signal between injections", {
  # AI-generated regression for the later association/dissociation call sites.
  check_mass_transport_simulation(TRUE)
})
