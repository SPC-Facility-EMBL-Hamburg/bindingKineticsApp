box::use(
  reticulate[
    py_to_r
  ],
  rhandsontable[
    hot_col,
    hot_table,
    rhandsontable
  ],
  shiny[
    req
  ],
  shinyalert[
    shinyalert
  ]
)

#' @export
pandas_to_r <- function(py_df) {
  if (is.null(py_df)) {
    NULL
  }

  py_to_r(py_df)
}

#' @export
sample_type_to_letter <- function(string) {
  dict <- list(
    "KSAMPLE"        = "S",
    "SAMPLE"         = "S",
    "LOADING"        = "L",
    "KLOADING"       = "L",
    "NEUTRALIZATION" = "N",
    "REGENERATION"   = "R",
    "BASELINE"       = "B",
    "BUFFER"         = "B"
  )

  if (!string %in% names(dict)) {
    # Return the first character of the string
    substr(string, 1, 1)
  } else {
    dict[[string]]
  }
}

#' @export
subset_data <- function(data, max_points = 200) {
  total_points <- length(data)
  if (total_points > max_points) {
    step <- ceiling(total_points / max_points)
    data <- data[seq(1, total_points, by = step)]
  }
  data
}

#' @export
df_to_lines <- function(df) {
  numeric_columns <- sapply(df, is.numeric)

  original_colnames <- colnames(df)

  if (sum(numeric_columns) > 0) {
    df[, numeric_columns] <- signif(df[, numeric_columns], 5)

    # Sort from high to low using the first column
    df <- as.data.frame(df)
    df <- df[rev(order(df[, which(numeric_columns)[1]])), ]
  }

  df <- data.frame(lapply(df, as.character), stringsAsFactors = FALSE)

  # Assign original column names
  colnames(df) <- original_colnames

  # Create a character vector with custom text
  lines <- c(paste0(colnames(df), collapse = ","))

  for (i in 1:nrow(df)) {
    lines <- c(lines, paste0(df[i, ], collapse = ","))
  }

  # Specify the file and use the cat() function to create the text file
  lines
}

## Find the index of the element inside a list that has the lowest maximal value
#' @export
find_probable_baseline <- function(ys_lst) {
  first_element <- ys_lst[[1]]

  min_val <- max(unlist(first_element), na.rm = TRUE)
  min_idx <- 1

  if (length(ys_lst) == 1) {
    1
  }

  for (i in 2:length(ys_lst)) {
    current_element <- ys_lst[[i]]
    current_val <- max(unlist(current_element), na.rm = TRUE)

    if (current_val < min_val) {
      min_val <- current_val
      min_idx <- i
    }
  }

  min_idx
}

#' @export
add_all_option <- function(options) {
  if (length(options) == 1) {
    options
  } else {
    c(options, "All")
  }
}
