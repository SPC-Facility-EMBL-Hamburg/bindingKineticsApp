box::use(
  . / colors[
    get_palette
  ],
  rhandsontable[
    hot_col,
    hot_table,
    rhandsontable
  ],
  shiny[
    req
  ]
)

renderer <- "function(instance, td, row, col, prop, value, cellProperties) {
                Handsontable.renderers.TextRenderer.apply(this, arguments);
                if (instance.params) {
                    hcols = instance.params.col_highlight
                    hcols = hcols instanceof Array ? hcols : [hcols]
                    hrows = instance.params.row_highlight
                    hrows = hrows instanceof Array ? hrows : [hrows]

                    for (i = 0; i < hcols.length; i++) {
                        if (hcols[i] == col && hrows[i] == row) {
                            td.style.background = instance.getDataAtCell(row, col);
                        }
                    }
                }
  }"

rendererboolean <- "function(instance, td, row, col, prop, value, cellProperties) {
            Handsontable.renderers.CheckboxRenderer.apply(this, arguments);}"

#' @export
get_plotting_df <- function(ID, labels = NULL) {
  if (is.null(labels)) {
    labels <- ID
  }

  n_total <- length(labels)

  df <- data.frame(
    "Internal_ID" = ID,
    "Color" = get_palette(n_total),
    "Legend" = labels,
    "Show" = as.logical(rep(TRUE, n_total))
  )
  df
}

# Auxiliary function to render the legend dataframe
# Requires
# - The dataframe
# Output
# - The rendered table to show in the server side
## Auxiliary function to render the legend dataframe
#' @export
get_rtable_legend <- function(legend_df,
                              fix_column = 1,
                              hex_color_column = 2,
                              extra_columns = c(3),
                              boolean_column = 4) {
  req(legend_df)

  color_cells <- data.frame(col = hex_color_column, row = seq_len(nrow(legend_df)))

  rtable <- rhandsontable(legend_df,
    rowHeaders = NULL,
    col_highlight = color_cells$col - 1,
    row_highlight = color_cells$row - 1
  )

  rtable <- rtable |>
    hot_col(
      col = c(fix_column), readOnly = TRUE,
      renderer = renderer
    ) |>
    hot_col(col = c(hex_color_column, extra_columns), renderer = renderer) |>
    hot_col(col = c(boolean_column), renderer = rendererboolean) |>
    hot_table(stretchH = "all")

  rtable
}

#' @export
get_sensor_df <- function(names, exp = NULL) {
  df <- data.frame("ID" = names, "Select" = rep(TRUE, length(names)))
  if (!is.null(exp)) {
    df$Experiment <- rep(exp, length(names))
  }
  df
}

#' @export
get_rtable_processing <- function(names) {
  df <- get_sensor_df(names)
  rdf <- rhandsontable(df) |>
    hot_col("Select") |>
    hot_table(stretchH = "all") |>
    hot_col("ID", readOnly = TRUE)

  rdf
}

## Auxiliary function to render the ligand concentration dataframe
## for surface-based binding
#' @export
render_ligand_conc_df_surface <- function(df) {
  # Convert specified columns to integer
  columns_to_convert <- c(
    "Analyte_location", "Loading_location", "Replicate", "Smax_ID"
  )

  for (column in columns_to_convert) {
    if (column %in% colnames(df)) {
      df[[column]] <- as.integer(as.character(df[[column]]))
    }
  }

  # Create the vector of columns that can not be edited
  # These vector includes the columns 'Sensor', 'Analyte_location'
  # 'Loading_location', 'Experiment' and 'Replicate'
  possible_fixed_columns <- c(
    "Sensor", "Analyte_location", "Loading_location",
    "Experiment", "Replicate"
  )

  fixed_columns_ids <- c()
  for (col in possible_fixed_columns) {
    if (col %in% colnames(df)) {
      idx <- which(colnames(df) == col)
      fixed_columns_ids <- c(fixed_columns_ids, idx)
    }
  }

  # find the index of the column 'Concentration_micromolar'
  idx <- which(colnames(df) == "[Analyte] (μM)")

  # convert to rhandsontable and set read only to columns 1,2,5,6
  # do not show row index
  df <- rhandsontable(df, rowHeaders = FALSE) |>
    hot_col(col = fixed_columns_ids, readOnly = TRUE) |>
    hot_col(col = idx, type = "numeric", format = "0.00000")

  df
}

## Auxiliary function to render the ligand concentration dataframe
##  for in-solution binding
#' @export
render_ligand_conc_df_solution <- function(df) {
  # Create the vector of columns that can not be edited
  possible_fixed_columns <- c("Trace", "Experiment")

  fixed_columns_ids <- c()
  for (col in possible_fixed_columns) {
    if (col %in% colnames(df)) {
      idx <- which(colnames(df) == col)
      fixed_columns_ids <- c(fixed_columns_ids, idx)
    }
  }

  # find the index of the column 'Protein_concentration_micromolar'
  idx1 <- which(colnames(df) == "[Protein] (μM)")
  # find the index of the column 'Ligand_concentration_micromolar'
  idx2 <- which(colnames(df) == "[Ligand] (μM)")

  # convert to rhandsontable and set read only to the columns 1,5
  # dont show row index
  df <- rhandsontable(df, rowHeaders = FALSE) |>
    hot_col(col = fixed_columns_ids, readOnly = TRUE) |>
    hot_col(col = c(idx1, idx2), type = "numeric", format = "0.0000")

  df
}
