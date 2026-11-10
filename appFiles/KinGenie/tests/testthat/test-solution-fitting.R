# AI-generated integration tests for solution fitting against pyKinGenie, without patched dependencies.
box::use(
  ../../box/modules/analysis/fit_controls[fitControlsServer],
  ../../box/modules/global_logbook[logbookServer]
)

solutionFittingTestServer <- function(id, state, pyKinetics) {
  # AI-generated: compose actual fitting and logging modules in a real Shiny session.
  shiny::moduleServer(id, function(input, output, session) {
    # AI-generated: supply the normal application module dependencies.
    logbook <- logbookServer("logbook")
    fitControlsServer("fit", state, shiny::reactive(data.frame()), pyKinetics, logbook)
  })
}

for (solution_model in c("one_binding_site", "one_binding_site_if")) {
  local({
    # AI-generated: keep each regression's model fixed.
    model <- solution_model
    testthat::test_that(paste(model, "forwards signal choices and clears exponential result state"), {
      # AI-generated: check UI parameter forwarding with six traces covering early and late dynamics.
      py <- reticulate::import("pykingenie")
      analyzer <- reticulate::import("pykingenie.main")$KineticsAnalyzer()
      times <- c(0, exp(seq(log(0.001), log(8), length.out = 60)))
      proteins <- c(0.4, 0.4, 0.4, 0.8, 0.8, 0.8)
      ligands <- rep(c(0.3, 0.9, 2.7), 2)
      signals <- lapply(seq_along(proteins), function(i) {
        # AI-generated: create input traces for the app workflow using the package as supplied.
        if (model == "one_binding_site") {
          py$signal_ode_one_site_insolution(reticulate::np_array(times), 1, 1.3,
            proteins[i], ligands[i], signal_a = 0.2, signal_b = 0.4, signal_complex = 1.7)
        } else {
          py$signal_ode_induced_fit_insolution(reticulate::np_array(times), list(0, 0),
            1, 1, 1, 1, proteins[i], ligands[i],
            signal_E = 0.2, signal_S = 0.4, signal_ES_int = 1.7, signal_ES = 2.3)
        }
      })
      fitting <- reticulate::import("pykingenie.fitter_solution")$KineticsFitterSolution(
        "sample", signals, as.list(ligands), as.list(proteins),
        rep(list(reticulate::np_array(times)), length(proteins)))
      analyzer$add_fitting(fitting, "sample")
      state <- shiny::reactiveValues(fit_dataset_loaded = TRUE, traces_loaded = FALSE,
        surface_based_binding = FALSE, solution_model = "double", k_obs_plot_shown = TRUE)
      shiny::testServer(solutionFittingTestServer, args = list(state = state, pyKinetics = analyzer), {
        session$setInputs(`fit-fittingRegion` = "association_dissociation",
          `fit-fittingModel` = "one_to_one", `fit-fit_signal_E` = FALSE, `fit-fit_signal_S` = FALSE,
          `fit-fit_signal_ES` = FALSE, `fit-fit_t0` = FALSE, `fit-ESint_equals_ES` = FALSE)
        if (model == "one_binding_site") session$setInputs(`fit-submitOBS` = 1)
        else session$setInputs(`fit-submitIF` = 1)
        testthat::expect_null(fitting$fit_params_kinetics)
        testthat::expect_equal(state$solution_model, "double")
        session$setInputs(`fit-fit_signal_E` = TRUE, `fit-fit_signal_S` = TRUE, `fit-fit_signal_ES` = TRUE)
        for (invalid in c(0, NaN)) {
          fitting$lig_conc <- as.list(replace(ligands, 1, invalid))
          if (model == "one_binding_site") {
            session$setInputs(`fit-submitOBS` = session$input[["fit-submitOBS"]] + 1)
          } else {
            session$setInputs(`fit-submitIF` = session$input[["fit-submitIF"]] + 1)
          }
          testthat::expect_null(fitting$fit_params_kinetics)
          testthat::expect_equal(state$solution_model, "double")
        }
        fitting$lig_conc <- as.list(ligands)
        if (model == "one_binding_site") session$setInputs(`fit-submitOBS` = 4)
        else session$setInputs(`fit-submitIF` = 4)
        testthat::expect_true(state$kinetics_fit_done)
        testthat::expect_equal(state$solution_model, model)
        testthat::expect_false(state$k_obs_plot_shown)
        fitted <- fitting$signal_assoc_fit
        testthat::expect_length(fitted, length(signals))
        for (i in seq_along(signals)) {
          testthat::expect_length(fitted[[i]], length(times))
          testthat::expect_true(all(is.finite(fitted[[i]])))
        }
        names <- colnames(fitting$fit_params_kinetics)
        testthat::expect_true(any(grepl("free protein", names)))
        testthat::expect_true(any(grepl("free ligand", names)))
        if (model == "one_binding_site_if") {
          testthat::expect_true(any(grepl("intermediate", names, ignore.case = TRUE)))
        }
      })
    })
  })
}
