## Welcome message
welcomeMessage <- function() {
      shinyalert(
          paste("Welcome to", appName," <br><small>
          KinGenie is being developed by the Sample Preparation and Characterisation (SPC)
          Facility at EMBL Hamburg.
          Please contact us if you have any questions, suggestions or feature requests. All
          comments are welcome. You can reach out us via spc@embl-hamburg.de or
          oburastero@embl-hamburg.de<br>"
        ),
     imageWidth = 180,imageHeight = 180,closeOnClickOutside=FALSE,closeOnEsc=FALSE,
     confirmButtonText="I accept",size = "m",
     showCancelButton=TRUE,cancelButtonText="I decline",html=TRUE,
     confirmButtonCol="#8bb8e8",
     callbackR = function(x) {
       if (!x) welcomeMessage()
     })
}

# Return a string with the pattern %H:%M
get_hour_minute_sec <- function() {
    time_str <- as.character(Sys.time())
    time_str <- strsplit(time_str,' ')[[1]][2]
    return(time_str)
}

popUpWarning <- function(string) shinyalert(text = string,type = "warning",closeOnEsc = T,closeOnClickOutside = T,html=T)
popUpInfo    <- function(string,size='m') shinyalert(text = string,type = "info",   closeOnEsc = T,closeOnClickOutside = T,html=T,size=size)
popUpSuccess <- function(string) shinyalert(text = string,type = "success",closeOnEsc = T,closeOnClickOutside = T,html=T)

sample_type_to_letter <- function(string) {

    dict <- list(
        "KSAMPLE"        = "S",
        "SAMPLE"         = "S",
        "LOADING"        = "L",
        "KLOADING"       = "L",
        "NEUTRALIZATION" = "N",
        "REGENERATION"   = "R",
        "BASELINE"       = "B",
        "BUFFER"         = "B")

    if (!string %in% names(dict)) {
        # Return the first character of the string
        return(substr(string,1,1))
    }

    return(dict[[string]])

}

subset_data <- function(data, max_points = 200) {
    total_points <- length(data)
    if (total_points > max_points) {
        step <- ceiling(total_points / max_points)
        data <- data[seq(1, total_points, by = step)]
    }
    return(data)
}

# Obtain a set of n colours (nColors) based on the Set 3 from the Brewer palette
get_palette <- function(n_colors) {

    if (n_colors <= 9) return(RColorBrewer::brewer.pal(9, "Set1")[1:n_colors])
    if (n_colors <= 12) return(RColorBrewer::brewer.pal(12, "Set3")[1:n_colors])

    return( colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_colors) )

}

get_plotting_df <- function(ID,labels=NULL) {

    if (is.null(labels)) {
        labels <- ID
    }

    n_total        <- length(labels)

    df <- data.frame("Internal_ID"= ID,
                   "Color"        = get_palette(n_total),
                   "Legend"       = labels,
                   "Show"         = as.logical(rep(TRUE,n_total))
                   )

    return(df)

}

# Auxiliary function to render the legend dataframe
# Requires
# - The dataframe
# Output
# - The rendered table to show in the server side
render_RHandsontable <- function(legend_df,fix_column=1,hex_color_column=2,extra_columns=c(3),boolean_column=4) {

  color_cells <- data.frame(col=hex_color_column,row=1:nrow(legend_df))

  rtable <- rhandsontable(legend_df,rowHeaders=NULL,
                          col_highlight = color_cells$col - 1,
                          row_highlight = color_cells$row - 1
  )

  rtable <- rtable %>% hot_col(col = c(fix_column), readOnly=TRUE,
                               renderer = myrenderer) %>%
    hot_col(col = c(hex_color_column,extra_columns),renderer = myrenderer) %>%
    hot_col(col = c(boolean_column),renderer = myrendererBoolean) %>%
    hot_table(stretchH='all')

  return(renderRHandsontable({rtable}))

}

get_sensor_df <- function(names,exp=NULL) {
  df  <- data.frame('ID' = names,'Select' = rep(TRUE,length(names)))
  if (!is.null(exp)) {
    df$Experiment <- rep(exp,length(names))
  }
  return(df)
}

get_rtable_processing <- function(names) {

    df  <- get_sensor_df(names)
    rdf <- rhandsontable(df) %>%
        hot_col('Select') %>%
        hot_table(stretchH='all')  %>%
        hot_col('ID',readOnly=TRUE)

    return(rdf)
}

## Convert dataframe to lines
df_to_lines <- function(df) {

  numeric_columns <- sapply(df, is.numeric)

  original_colnames <- colnames(df)

  if (sum(numeric_columns) > 0 ){

    df[, numeric_columns] <- signif(df[, numeric_columns], 5)

    # Sort from high to low using the first column
    df   <- as.data.frame(df)
    df <- df[rev(order(df[,which(numeric_columns)[1]])),]

  }

  df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE)

  # Assign original column names
  colnames(df) <- original_colnames

  # Create a character vector with custom text
  lines <- c(paste0(colnames(df),collapse = ','))

  for (i in 1:nrow(df)) {
    lines <- c(lines,paste0(df[i,],collapse = ','))
  }

  # Specify the file and use the cat() function to create the text file
  return(lines)

}

# Find the index of the element inside a list that has the lowest maximal value
# The list elements are sublists
find_probable_baseline <- function(ys_lst) {

    first_element <- ys_lst[[1]]

    min_val <- max(unlist(first_element),na.rm=TRUE)
    min_idx <- 1

    if (length(ys_lst) == 1) return(1)

    for (i in 2:length(ys_lst)) {

        current_element <- ys_lst[[i]]
        current_val     <- max(unlist(current_element),na.rm=TRUE)

        if (current_val < min_val) {

            min_val <- current_val
            min_idx <- i

        }

    }

    return(min_idx)

}

# Print message to user regarding the simulation model
print_model_message <- function(model_type_sim,model_selected_sim) {

    c1 <- model_type_sim == 'surface'

    model_typ_text <- ifelse(c1,'a surface-based binding model (constant ligand concentration)',
    'an in solution-based binding model (ligand depletion)')

    c2 <- model_selected_sim %in% c('one_site_induced_fit','one_site_conformational_selection')

    extra_txt <- ifelse(!c1 && c2,
    'Relaxation rates are computed using Equations taken from Fabian Paul, Thomas R. Weikl, 2016.<br>','')

    extra_txt2 <- ifelse(c1 && model_selected_sim == 'one_site_induced_fit',
    'In this model we assume that the induced and intermediate complex produce the same signal.<br>','')

    if (model_selected_sim == 'one_site') {

        popUpInfo(
            paste0("The simulation will be performed using ",model_typ_text,"<br>
            The chemical equation is P + L ⇄ PL.
            Please wait until the simulation is completed."),
            size = 'l'
        )

    }

    if (model_selected_sim == 'heterogeneous_analyte') {

        popUpInfo(
            paste0("The simulation will be performed using ",model_typ_text,"<br>
            The chemical equations are
            L + P<sub>1</sub> ⇄ P<sub>1</sub>L (binding to the analyte P<sub>1</sub>) and
            L + P<sub>2</sub> ⇄ P<sub>2</sub>L (binding to the analyte P<sub>2</sub>). <br>
            The chemical reactions are: <br>
            <span style='display: inline-block; text-align: center;'>
            k<sub>on,1</sub><br>
            L + P<sub>1</sub> → P<sub>1</sub>L <br>
            k<sub>off,1</sub><br>
            P<sub>1</sub>L → L + P<sub>1</sub> <br>
            k<sub>on,2</sub><br>
            L + P<sub>2</sub> → P<sub>2</sub>L <br>
            k<sub>off,2</sub><br>
            P<sub>2</sub>L → L + P<sub>2</sub> <br>
            </span> <br>
            Please wait until the simulation is completed."),
            size = 'l'
        )

    }

    if (model_selected_sim == 'one_site_mtl') {

        popUpInfo(
            paste0('<p style="line-height: 2;">
            The simulation will be performed using ',model_typ_text,'<br>
            The chemical equations are
            L<sub>bulk</sub> ⇄ L<sub>surface</sub> (mass transport step) and
            L<sub>surface</sub> + P ⇄ PL (binding step). <br>
            The chemical reactions are: <br>
            <span style="display: inline-block; text-align: center;">
            k<sub>tr</sub><br>
            L<sub>bulk</sub> → L<sub>surface</sub> <br>
            k<sub>tr</sub><br>
            L<sub>surface</sub> → L<sub>bulk</sub> <br>
            k<sub>on</sub><br>
            L<sub>surface</sub> + P → PL <br>
            k<sub>off</sub><br>
            PL → L<sub>surface</sub> + P <br>
            </span> <br>
            ',extra_txt,'
            Please wait until the simulation is completed.
            </p>'),
            size = 'l'
        )
    }

    if (model_selected_sim == 'one_site_induced_fit') {

        popUpInfo(
            paste0('<p style="line-height: 2;">
            The simulation will be performed using ',model_typ_text,'<br>
            The chemical equations are
            P + L ⇄ PL<sub>int</sub> (binding step) and PL<sub>int</sub> ⇄ PL (induced fit step). <br>
            The chemical reactions are: <br>
            <span style="display: inline-block; text-align: center;">
            k<sub>on</sub><br>
            P + L → PL<sub>int</sub> <br>
            k<sub>off</sub><br>
            PL<sub>int</sub> → P + L <br>
            k<sub>c</sub><br>
            PL<sub>int</sub> → PL <br>
            k<sub>rev</sub><br>
            PL → PL<sub>int</sub> <br>
            </span> <br>
            ',extra_txt,',
            ',extra_txt2,'
            Please wait until the simulation is completed.
            </p>'),
            size = 'l'
        )
    }

    if (model_selected_sim == 'one_site_conformational_selection') {

        popUpInfo(
            paste0('<p style="line-height: 2;">
            The simulation will be performed using ',model_typ_text,'<br>
            The chemical equations are
            P<sub>inactive</sub> ⇄ P<sub>active</sub> (conformational change step) and P<sub>active</sub> + L ⇄ PL (binding step). <br>
            The chemical reactions are: <br>
            <span style="display: inline-block; text-align: center;">
            k<sub>on</sub><br>
            P<sub>active</sub> + L → PL <br>
            k<sub>off</sub><br>
            PL → P<sub>active</sub> + L <br>
            k<sub>c</sub><br>
            P<sub>inactive</sub> → P<sub>active</sub> <br>
            k<sub>rev</sub><br>
            P<sub>active</sub> → P<sub>inactive</sub> <br>
            </span> <br>
            ',extra_txt,'
            Please wait until the simulation is completed.
            </p>'),
            size = 'l'
        )
    }


    return(NULL)

}

render_combined_ligand_conc_df <- function(df) {

    # Convert specified columns to integer
    columns_to_convert <- c("Analyte_location", "Loading_location", "Replicate","Smax_ID")
    for (column in columns_to_convert) {
        if (column %in% colnames(df)) {
            df[[column]] <- as.integer(as.character(df[[column]]))
        }
    }

    # Create the vector of columns that can not be edited
    # These vector includes the columns 'Sensor', 'Analyte_location'
    # 'Loading_location', 'Experiment' and 'Replicate'
    possible_fixed_columns <- c('Sensor','Analyte_location','Loading_location','Experiment','Replicate')

    fixed_columns_ids <- c()
    for (col in possible_fixed_columns) {

        if (col %in% colnames(df)) {

            idx               <- which(colnames(df) == col)
            fixed_columns_ids <- c(fixed_columns_ids, idx)
        }

    }

    # find the index of the column 'Concentration_micromolar'
    idx <- which(colnames(df) == '[Analyte] (μM)')

    # convert to rhandsontable and set read only to columns 1,2,5,6
    # do not show row index
    df <- rhandsontable(df,rowHeaders = FALSE) %>%
        hot_col(col = fixed_columns_ids, readOnly = TRUE) %>%
        hot_col(col = idx, type = 'numeric', format = '0.00000')

    return(df)

}

render_combined_ligand_conc_df_solution <- function(df) {

    # Create the vector of columns that can not be edited
    possible_fixed_columns <- c('Trace','Experiment')

    fixed_columns_ids <- c()
    for (col in possible_fixed_columns) {

        if (col %in% colnames(df)) {

            idx               <- which(colnames(df) == col)
            fixed_columns_ids <- c(fixed_columns_ids, idx)
        }

    }

    # find the index of the column 'Protein_concentration_micromolar'
    idx1 <- which(colnames(df) == '[Protein] (μM)')
    # find the index of the column 'Ligand_concentration_micromolar'
    idx2 <- which(colnames(df) == '[Ligand] (μM)')


    # convert to rhandsontable and set read only to the columns 1,5
    # dont show row index
    df <- rhandsontable(df,rowHeaders = FALSE) %>%
        hot_col(col = fixed_columns_ids, readOnly = TRUE) %>%
        hot_col(col = c(idx1,idx2), type = 'numeric', format = '0.0000')

    return(df)

}

add_all_option <- function(options) {

  if (length(options) == 1) {
    return(options)
  } else {
    return(c(options,'All'))
  }

}


