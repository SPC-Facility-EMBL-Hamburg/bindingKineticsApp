# AI-generated browser workflow following test-shinytest2-simulation-1.R.
# PR #4: https://github.com/SPC-Facility-EMBL-Hamburg/bindingKineticsApp/pull/4
library(shinytest2)
library(testthat)

TIMEOUT <- 60000

test_that("KinGenie mass-transport simulation and download workflow", {
  # AI-generated: exercise the visible app and compare saved Plotly and download references.
  # PR #4: https://github.com/SPC-Facility-EMBL-Hamburg/bindingKineticsApp/pull/4
  app <- AppDriver$new(
    test_path("../.."),
    variant = "ci", name = "KinGenie-mass-transport",
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

  set_input(list(
    "modelSelectionAndRun-model_selected_sim" = "one_site_mtl"
  ))
  wait_js("bindingParams-ktr_sim_1to1")
  set_input(list(
    "modelSelectionAndRun-time_step" = 0.5,
    "experimentalParams-init_lig_sim_surface" = 4,
    "experimentalParams-numb_dil_sim" = 2,
    "experimentalParams-lig_dil_factor_sim" = 2,
    "experimentalParams-protein_smax_sim" = 5,
    "experimentalParams-numb_dil_sim_prot" = 1,
    "experimentalParams-prot_dil_factor_sim" = 2,
    "experimentalParams-association_time" = 20,
    "experimentalParams-dissociation_time" = 30,
    "bindingParams-kd_sim_1to1" = 0.5,
    "bindingParams-koff_sim_1to1" = 0.01,
    "bindingParams-ktr_sim_1to1" = 0.25
  ))
  wait_for_idle()

  click("modelSelectionAndRun-btn_cal_simulation")
  accept_modal()
  wait_for_idle()
  expect_plotly_equal(
    "simPlotTabBox-signal_sim_plot",
    test_path("reference", "signal_sim_plot_41.rds")
  )
  app$expect_download("simResultsExport-btn_export_simulation")

  set_input(list("experimentalParams-is_single_cycle_sim" = TRUE))
  wait_for_idle()
  click("modelSelectionAndRun-btn_cal_simulation")
  accept_modal()
  wait_for_idle()
  expect_plotly_equal(
    "simPlotTabBox-signal_sim_plot",
    test_path("reference", "signal_sim_plot_42.rds")
  )
  app$expect_download("simResultsExport-btn_export_simulation")

  app$stop()
})
