library(shinytest2)
library(testthat)

TIMEOUT <- 10000

test_that("KinGenie simulation workflow 1", {
  app <- AppDriver$new(
    test_path("../.."),
    variant = "ci", name = "KinGenie-simulation",
    height = 756, width = 1400,
    load_timeout = 20000
  )
  
  accept_modal <- function() {
    app$wait_for_js("document.querySelector('.confirm') !== null")
    app$run_js("document.querySelector('.confirm').click()")
    app$wait_for_idle()
    invisible(NULL)
  }
  
  click <- function(button = NULL, selector = NULL) {
    if (!is.null(selector)) {
      app$click(selector = selector)
    } else {
      app$click(button)
    }
    invisible(NULL)
  }
  
  set_input <- function(inputs, ...) {
    do.call(app$set_inputs, c(inputs, list(...)))
    invisible(NULL)
  }
  
  wait_js <- function(element) {
    app$wait_for_js(
      paste0("document.getElementById('", element, "') !== null"),
      timeout = TIMEOUT
    )
    invisible(NULL)
  }
  
  wait_for_idle <- function() {
    app$wait_for_idle()
    invisible(NULL)
  }

  screenshot <- function() {
    app$expect_screenshot()
    invisible(NULL)
  }

  expect_values <- function(module_name,outputs) {
    app$expect_values(
        output = paste0(module_name, "-", outputs)
    )
  }

  clean_plotly_value <- function(json, digits = 10) {

    obj <- jsonlite::fromJSON(
      json,
      simplifyVector = FALSE
    )

    # Keep only the Plotly specification
    obj <- obj$x

    recurse <- function(x) {

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

  click("modelSelectionAndRun-btn_cal_simulation")
  accept_modal()
  wait_for_idle()

  expect_plotly_equal(
    "simPlotTabBox-signal_sim_plot",
    test_path("reference", "signal_sim_plot_1.rds")
  )

  set_input(list(
    "modelSelectionAndRun-model_type_sim" = "solution"
  ))
  wait_for_idle()
  set_input(
    list(
      "modelSelectionAndRun-model_selected_sim" = "one_site_conformational_selection"
    )
  )
  wait_for_idle()
  click("modelSelectionAndRun-btn_cal_simulation")
  accept_modal()
  wait_for_idle()
  expect_plotly_equal(
    "simPlotTabBox-signal_sim_plot",
    test_path("reference", "signal_sim_plot_2.rds")
  )

  app$expect_download("simResultsExport-btn_export_simulation")

  app$stop()
})


