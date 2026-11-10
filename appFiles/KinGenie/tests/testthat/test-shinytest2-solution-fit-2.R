# AI-generated browser workflow following test-shinytest2-simulation-1.R.
# PR #4: https://github.com/SPC-Facility-EMBL-Hamburg/bindingKineticsApp/pull/4
library(shinytest2)
library(testthat)

TIMEOUT <- 60000

test_that("KinGenie solution induced fit fitting workflow", {
  # AI-generated: exercise the visible app and compare saved Plotly and download references.
  # PR #4: https://github.com/SPC-Facility-EMBL-Hamburg/bindingKineticsApp/pull/4
  app <- AppDriver$new(
    test_path("../.."),
    variant = "ci", name = "KinGenie-solution-fit-2",
    height = 756, width = 1400,
    load_timeout = 20000, timeout = TIMEOUT
  )
  on.exit(app$stop(), add = TRUE)

  accept_modal <- function() {
    # AI-generated adaptation of the existing browser-test helper.
    # PR #4: https://github.com/SPC-Facility-EMBL-Hamburg/bindingKineticsApp/pull/4
    app$wait_for_js("document.querySelector('.confirm') !== null")
    app$run_js("document.querySelector('.confirm').click()")
    app$wait_for_idle()
    invisible(NULL)
  }

  click <- function(button = NULL, selector = NULL) {
    # AI-generated adaptation of the existing browser-test helper.
    # PR #4: https://github.com/SPC-Facility-EMBL-Hamburg/bindingKineticsApp/pull/4
    if (!is.null(selector)) {
      app$click(selector = selector)
    } else {
      app$click(button)
    }
    invisible(NULL)
  }

  set_input <- function(inputs, ...) {
    # AI-generated adaptation of the existing browser-test helper.
    # PR #4: https://github.com/SPC-Facility-EMBL-Hamburg/bindingKineticsApp/pull/4
    do.call(app$set_inputs, c(inputs, list(...)))
    invisible(NULL)
  }

  wait_js <- function(element) {
    # AI-generated adaptation of the existing browser-test helper.
    # PR #4: https://github.com/SPC-Facility-EMBL-Hamburg/bindingKineticsApp/pull/4
    app$wait_for_js(
      paste0("document.getElementById('", element, "') !== null"),
      timeout = TIMEOUT
    )
    invisible(NULL)
  }

  wait_for_idle <- function() {
    # AI-generated adaptation of the existing browser-test helper.
    # PR #4: https://github.com/SPC-Facility-EMBL-Hamburg/bindingKineticsApp/pull/4
    app$wait_for_idle()
    invisible(NULL)
  }

  expect_values <- function(module_name,outputs) {
    # AI-generated adaptation of the existing browser-test helper.
    # PR #4: https://github.com/SPC-Facility-EMBL-Hamburg/bindingKineticsApp/pull/4
    app$expect_values(
        output = paste0(module_name, "-", outputs)
    )
  }

  clean_plotly_value <- function(json, digits = 10) {
    # AI-generated adaptation of the existing browser-test helper.
    # PR #4: https://github.com/SPC-Facility-EMBL-Hamburg/bindingKineticsApp/pull/4

    obj <- jsonlite::fromJSON(
      json,
      simplifyVector = FALSE
    )

    # Keep only the Plotly specification
    obj <- obj$x

    recurse <- function(x) {
      # AI-generated adaptation of the existing browser-test helper.
      # PR #4: https://github.com/SPC-Facility-EMBL-Hamburg/bindingKineticsApp/pull/4

      if (is.numeric(x))
        return(round(x, digits))

      if (!is.list(x))
        return(x)

      # Remove fields that change between versions
      x$dependencies <- NULL
      x$elementId <- NULL
      x$jsHooks <- NULL
      x$attrs <- NULL
      x$source <- NULL
      x$config <- NULL
      x$visdat <- NULL
      x$cur_data <- NULL
      x$cur_data_all <- NULL
      x$highlight <- NULL
      x$base_url <- NULL

      x$uid <- NULL
      x$key <- NULL
      x$frame <- NULL

      lapply(x, recurse)
    }

    recurse(obj)
  }


  expect_plotly_equal <- function(output, file,
                                  tolerance = 1e-7) {
    # AI-generated adaptation of the existing browser-test helper.
    # PR #4: https://github.com/SPC-Facility-EMBL-Hamburg/bindingKineticsApp/pull/4

    actual <- clean_plotly_value(
      app$get_value(output = output)
    )

    dir.create(dirname(file),
              recursive = TRUE,
              showWarnings = FALSE)

    if (!file.exists(file)) {

      saveRDS(actual, file)

      testthat::fail(
        paste(
          "Created reference:",
          file,
          "\nInspect it, commit it, then rerun the tests."
        )
      )
    }

    expected <- readRDS(file)

    # Helpful diagnostics if something changes
    diff <- waldo::compare(
      actual,
      expected,
      tolerance = tolerance
    )

    if (length(diff)) {

      cat(
        "\n========== Plotly differences ==========\n",
        paste(diff, collapse = "\n"),
        "\n=========================================\n"
      )

    }

    expect_equal(
      actual,
      expected,
      tolerance = tolerance
    )
  }

  accept_modal()

  click(selector = "a[data-value='menu_simulation']")
  app$wait_for_js(
    "document.querySelector('#shiny-tab-menu_simulation.active') !== null",
    timeout = TIMEOUT
  )
  wait_for_idle()

  set_input(list("modelSelectionAndRun-model_type_sim" = "solution"))
  wait_for_idle()
  set_input(list("modelSelectionAndRun-model_selected_sim" = "one_site_induced_fit"))
  wait_for_idle()
  set_input(list(
    "modelSelectionAndRun-time_step" = 0.05,
    "experimentalParams-total_time" = 8,
    "experimentalParams-init_lig_sim_solution" = 2.7,
    "experimentalParams-numb_dil_sim" = 2,
    "experimentalParams-lig_dil_factor_sim" = 3,
    "experimentalParams-protein_conc_sim" = 0.8,
    "experimentalParams-numb_dil_sim_prot" = 1,
    "experimentalParams-prot_dil_factor_sim" = 2,
    "bindingParams-signal_E" = 0,
    "bindingParams-signal_S" = 0,
    "bindingParams-kon_sim_1to1_adv" = 1,
    "bindingParams-koff_sim_1to1_adv" = 1,
    "bindingParams-kc_sim_1to1_adv" = 1,
    "bindingParams-krev_sim_1to1_adv" = 1,
    "bindingParams-signal_ES_int" = 1.7,
    "bindingParams-signal_ES" = 2.3
  ))
  click("modelSelectionAndRun-btn_cal_simulation")
  accept_modal()
  wait_for_idle()

  csv_file <- app$get_download("simResultsExport-btn_export_simulation")
  # Use the filename produced by the app outside Shiny test mode.
  simulation_file <- file.path(tempdir(), "simulation_KinGenie_solution.csv")
  file.copy(csv_file, simulation_file, overwrite = TRUE)
  click(selector = "a[data-value='menu_import']")
  wait_for_idle()
  app$upload_file("load-kineticFiles" = simulation_file)
  wait_js("load-submitLoadExperiment")
  set_input(list("load-newExperimentName" = "Solution"))
  click("load-submitLoadExperiment")
  wait_for_idle()

  click(selector = "a[data-value='menu_analyse']")
  wait_js("fitControls-triggerCreateDatasetSolution")
  click("fitControls-triggerCreateDatasetSolution")
  wait_for_idle()
  expect_plotly_equal(
    "fitResults-tracesFitSolution",
    test_path("reference", "traces_solution_induced_fit_input.rds")
  )

  set_input(list("fitControls-modelSelectionSolution" = "double"))
  click("fitControls-triggerFitSolution")
  accept_modal()
  app$wait_for_js(
    "document.querySelector(`a[data-value='Observed constants']`) !== null",
    timeout = TIMEOUT
  )

  set_input(list("fitControls-modelSelectionSolution" = "one_binding_site_if"))
  click("fitControls-triggerFitSolution")
  wait_js("fitControls-submitIF")
  set_input(list(
    "fitControls-fit_signal_E" = FALSE,
    "fitControls-fit_signal_S" = FALSE,
    "fitControls-fit_signal_ES" = FALSE
  ))
  click("fitControls-submitIF")
  accept_modal()
  wait_js("fitControls-submitIF")
  set_input(list(
    "fitControls-fit_signal_E" = FALSE,
    "fitControls-fit_signal_S" = FALSE,
    "fitControls-fit_signal_ES" = TRUE,
    "fitControls-ESint_equals_ES" = FALSE
  ))
  click("fitControls-submitIF")
  accept_modal()
  wait_for_idle()
  expect_plotly_equal(
    "fitResults-tracesFitSolution",
    test_path("reference", "traces_solution_induced_fit_fit.rds"),
    tolerance = 1e-4
  )
  app$wait_for_js(
    "document.querySelector(`a[data-value='Observed constants']`) === null",
    timeout = TIMEOUT
  )

  # Dismiss the completion notice if the progress notice was shown first.
  app$run_js("{ const button = document.querySelector('.sweet-alert.visible .confirm'); if (button) button.click(); }")
  wait_for_idle()
  click(selector = "a[data-value='Fitted params (Kinetics)']")
  wait_for_idle()
  expect_values("fitResults", c("fittingInfoKinetics"))

  app$stop()
})
