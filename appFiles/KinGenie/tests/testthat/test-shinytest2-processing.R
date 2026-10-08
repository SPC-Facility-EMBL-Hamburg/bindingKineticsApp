library(shinytest2)
library(testthat)

TIMEOUT <- 60000

test_that("KinGenie proccessing workflow", {
  app <- AppDriver$new(
    test_path("../.."),
    variant = "ci", name = "KinGenie-surface-processing",
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

  load_file <- function(module_name,input,file_path) {

    labeled_input <- paste0(module_name, "-", input) 
    normalized_input <- path.expand(file_path)

    if (grepl("^/", normalized_input)) {
      resolved_path <- normalizePath(normalized_input, mustWork = TRUE)
    } else {
      app_root <- test_path("../..")
      relative_path <- sub("^[.]/", "", file_path)
      resolved_path <- normalizePath(
        file.path(app_root, relative_path),
        mustWork = TRUE
      )
    }

    named_path <- list()
    named_path[[labeled_input]] <- resolved_path
    do.call(app$upload_file, named_path)
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
  load_file(
    module_name = "load", 
    input = "kineticFiles", 
    file_path = "./www/test_bli_folder/230309_001.frd"
    )

  wait_js("load-submitLoadExperiment")

  set_input(list(
    "load-newExperimentName" = 'Test-1'
  ))

  click("load-submitLoadExperiment")

  wait_plot_input_update()

  expect_plotly_equal(
    "plotInput-traces",
    test_path("reference", "traces_processing.rds")
  )

  load_file(
    module_name = "load", 
    input = "kineticFiles", 
    file_path = "./www/test_bli_folder/230309_001.frd"
    )

  wait_js("load-submitLoadExperiment")

  set_input(list(
    "load-newExperimentName" = 'Test-2'
  ))

  click("load-submitLoadExperiment")

  wait_plot_input_update()

  expect_plotly_equal(
    "plotInput-traces",
    test_path("reference", "traces_processing_2.rds")
  )

  set_input(list("processing-selectedExperiment" = "All"))
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

  expect_plotly_equal(
    "plotInput-traces",
    test_path("reference", "traces_processing_3.rds")
  )

  app$stop()
})


