library(shinytest2)
library(testthat)

TIMEOUT <- 60000

test_that("KinGenie load example, create fitting dataset and try screening mode", {
  app <- AppDriver$new(
    test_path("../.."),
    variant = "ci", name = "KinGenie-surface-screening",
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

  edit_ligand_info_for_screening <- function() {
    wait_js("fitData-ligandInfo")
    hot_input <- app$get_value(input = "fitData-ligandInfo")

    if (is.null(hot_input)) {
      stop("fitData-ligandInfo input is NULL")
    }

    if (is.data.frame(hot_input)) {
      df <- hot_input
      n <- min(8, nrow(df))

      if (n > 0) {
        df$SampleID[seq_len(n)] <- as.character(seq_len(n))
      }
      df$Select <- FALSE
      if (n > 0) {
        df$Select[seq_len(n)] <- TRUE
      }

      set_input(list("fitData-ligandInfo" = df))
    } else if (is.list(hot_input) && !is.null(hot_input$data)) {
      headers <- unlist(hot_input$params$rColHeaders)
      sample_col <- match("SampleID", headers)
      select_col <- match("Select", headers)

      if (is.na(sample_col) || is.na(select_col)) {
        stop("SampleID or Select column not found in fitData-ligandInfo payload")
      }

      n_rows <- length(hot_input$data)
      n <- min(8, n_rows)

      for (i in seq_len(n_rows)) {
        row_i <- hot_input$data[[i]]
        if (i <= n) {
          row_i[[sample_col]] <- as.character(i)
          row_i[[select_col]] <- TRUE
        } else {
          row_i[[select_col]] <- FALSE
        }
        hot_input$data[[i]] <- row_i
      }

      if (!is.null(hot_input$changes)) {
        hot_input$changes$event <- "afterChange"
      } else {
        hot_input$changes <- list(event = "afterChange")
      }

      set_input(list("fitData-ligandInfo" = hot_input))
    } else {
      stop("Unexpected fitData-ligandInfo input structure")
    }

    app$run_js("(function(){
      const el = document.getElementById('fitData-ligandInfo');
      if (!el) throw new Error('fitData-ligandInfo not found');

      const widget = HTMLWidgets.getInstance(el);
      if (!widget || !widget.hot) throw new Error('Handsontable instance not found');

      const hot = widget.hot;
      const sampleCol = 2; // 0-based: third column is SampleID
      const selectCol = 3; // 0-based: fourth column is Select
      const nRows = hot.countRows();
      const limit = Math.min(8, nRows);

      for (let r = 0; r < limit; r++) {
        hot.setDataAtCell(r, sampleCol, String(r + 1), 'edit');
        hot.setDataAtCell(r, selectCol, true, 'edit');
      }

      for (let r = limit; r < nRows; r++) {
        hot.setDataAtCell(r, selectCol, false, 'edit');
      }
    })();")

    app$wait_for_idle()

    hot_after <- app$get_value(input = "fitData-ligandInfo")
    if (is.data.frame(hot_after)) {
      df_after <- hot_after
    } else if (is.list(hot_after) && !is.null(hot_after$data)) {
      headers_after <- unlist(hot_after$params$rColHeaders)
      mat_after <- do.call(rbind, lapply(hot_after$data, unlist))
      df_after <- as.data.frame(mat_after, stringsAsFactors = FALSE)
      if (length(headers_after) == ncol(df_after)) {
        colnames(df_after) <- headers_after
      }
      if ("Select" %in% colnames(df_after)) {
        df_after$Select <- as.logical(df_after$Select)
      }
    } else {
      stop("Unexpected fitData-ligandInfo input structure after edit")
    }

    n_after <- min(8, nrow(df_after))
    if (n_after > 0) {
      expect_true(all(as.character(df_after$SampleID[seq_len(n_after)]) == as.character(seq_len(n_after))))
      expect_true(all(df_after$Select[seq_len(n_after)]))
    }
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
  
  click(selector = "a[data-value='menu_analyse']")
  app$wait_for_js(
    "document.querySelector('#shiny-tab-menu_analyse.active') !== null",
    timeout = TIMEOUT
  )
  wait_for_idle()

  edit_ligand_info_for_screening()
  
  click("fitControls-triggerCreateDataset")
  wait_for_idle()
  accept_modal()

  click("visualizationConfigFit-configure_visualization")
  wait_js("visualizationConfigFit-visualization_screen_mode")
  set_input(list("visualizationConfigFit-visualization_screen_mode" = TRUE))
  click(selector = ".modal-footer .btn-default")
  wait_for_idle()

  # screenshot() # run this command only to verify the intended behaviour

  expect_plotly_equal(
    "fitResults-tracesAssDissFit",
    test_path("reference", "tracesAssDissScreeningMode.rds")
  )

  app$stop()
})


