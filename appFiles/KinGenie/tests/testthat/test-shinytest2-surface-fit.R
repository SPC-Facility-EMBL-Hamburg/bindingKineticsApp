library(shinytest2)
library(testthat)

TIMEOUT <- 60000

test_that("KinGenie load example and fit workflow", {
  app <- AppDriver$new(
    test_path("../.."),
    variant = "ci", name = "KinGenie-surface-fit",
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

  wait_plot_input_update <- function() {
    app$wait_for_js("
      const el = document.getElementById('plotInput-traces');
      el && !el.classList.contains('recalculating');
    ", timeout = TIMEOUT)
    wait_for_idle()
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
  click("load-loadExampleData")
  set_input(list("processing-operation" = "align_association"))
  click("processing-triggerProcessing")
  wait_js("processing-submitAlign")
  
  set_input(list(
    "processing-inPlaceAlignment" = TRUE,
    "processing-createNewSensorNames" = FALSE,
    "processing-keepRegeneration" = FALSE,
    "processing-keepLoading" = FALSE,
    "processing-keepBaseline" = FALSE,
    "processing-keepActivation" = FALSE,
    "processing-keepQuenching" = FALSE,
    "processing-keepCustom" = FALSE
  ))
  click("processing-submitAlign")
  wait_plot_input_update()
  set_input(list("processing-operation" = "correct_dissociation"))
  click("processing-triggerProcessing")
  wait_js("processing-submitCorrectDis")
  set_input(list("processing-nPointsCorrectDis" = 1, "processing-inPlaceCorrection" = TRUE))
  click(selector = "#processing-submitCorrectDis")
  wait_plot_input_update()
  
  set_input(list("processing-operation" = "subtract"))
  click("processing-triggerProcessing")
  wait_js("processing-submitSubtraction1")
  set_input(list("processing-inputBaseline" = "H1"))
  click("processing-submitSubtraction1")
  wait_for_idle()
  
  wait_js("processing-submitSubtraction2")
  set_input(list("processing-inPlaceSubtraction" = TRUE))
  click("processing-submitSubtraction2")
  wait_plot_input_update()
  expect_plotly_equal(
    "plotInput-traces",
    test_path("reference", "traces.rds")
  )
  
  click(selector = "a[data-value='menu_analyse']")
  app$wait_for_js(
    "document.querySelector('#shiny-tab-menu_analyse.active') !== null",
    timeout = TIMEOUT
  )
  wait_for_idle()
  
  click("fitData-quickSelect")
  wait_js("fitData-quickSelectOK")
  set_input(list("fitData-sample_ids_fit" = c("wt - imd")))
  click("fitData-quickSelectOK")
  wait_for_idle()
  
  click("fitControls-triggerCreateDataset")
  wait_for_idle()
  accept_modal()
  expect_plotly_equal(
    "fitResults-tracesAssDissFit",
    test_path("reference", "tracesAssDissFit1.rds")
  )
  
  set_input(list("fitControls-fittingRegion" = "association_dissociation"))
  wait_for_idle()
  set_input(list("fitControls-fittingModel" = "one_to_one_if"))
  click("fitControls-triggerFitting")
  wait_js("fitControls-submitKineticsFitting")
  set_input(list("fitControls-linkedRmax" = TRUE))
  click("fitControls-submitKineticsFitting")
  wait_for_idle()
  accept_modal()
  expect_plotly_equal(
    "fitResults-tracesAssDissFit",
    test_path("reference", "tracesAssDissFit2.rds")
  )
  app$click(selector = "a[data-value='kineticsParamsTable']")
  wait_for_idle()
  expect_values("fitResults",c("fittingInfoKinetics"))
  app$stop()
})


