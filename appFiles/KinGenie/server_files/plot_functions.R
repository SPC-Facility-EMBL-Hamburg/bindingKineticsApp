plot_ticks_length_cst <- 8

get_colors_from_numeric_values <- function(values,minVal,maxVal,useLogScale=TRUE) {

    viridis <- c("#440154","#471063","#481d6f","#472a7a",
               "#414487","#3c4f8a","#375a8c",
               "#32648e","#2a788e","#26828e",
               "#228b8d","#1f958b","#22a884",
               "#2cb17e","#3bbb75","#4ec36b",
               "#7ad151","#95d840","#b0dd2f","#cae11f",
               "#fde725")

    if (useLogScale) {

        minVal <- log10(minVal)
        maxVal <- log10(maxVal)

    }

    seq <- seq(minVal,maxVal,length.out = 21)

    if (useLogScale) {
        idx <- sapply(values,function(v) which.min(abs(log10(v) - seq)))
    } else {
        idx <- sapply(values,function(v) which.min(abs(v - seq)))
    }

  return(viridis[idx])
}

config_fig <- function(fig,nameForDownload,plot_type,plot_width,plot_height) {
    fig %>%  config(
        toImageButtonOptions = list(
        format = plot_type,
        filename = nameForDownload,
        width = plot_width * 50,
        height = plot_height * 50
        ), displaylogo = FALSE,
        modeBarButtonsToRemove = list(
            'hoverClosestCartesian',
            'hoverCompareCartesian',
            'lasso2d','select2d')
        )
}

plot_list_to_fig <- function(plot_name,plot_list,title_text,axis_size,
                             plot_type,plot_width,plot_height,nrows=2,share_axis = TRUE) {

    i <- length(plot_list)

    if (i > 1) {

        if (share_axis) {

            fig <- subplot(plot_list,nrows = nrows, margin = c(0.03,0.03,0.1,0.2),shareY = TRUE, shareX = TRUE)

        } else {

            fig <- subplot(plot_list,nrows = nrows,margin = c(0.03,0.03,0.05,0.05),
            shareY = FALSE, shareX = FALSE,titleX = TRUE,titleY = TRUE)

        }

    } else {

        fig <- plot_list[[1]] %>%  layout(title=list(text=title_text,font = list(size = axis_size) ) )

    }

    fig <- config_fig(fig,plot_name,plot_type,plot_width,plot_height)

    return(fig)

}

# This function will:
#  add the x axis to the plot           argument 'xaxis'
#  add the y axis to the plot           argument 'yaxis'
#  add the legend                       argument 'leg'
#  at a certain position                argument 'x_pos_annot'
#  if there is more than one condition  argument 'tot_cond'
# set the size font                     argument 'axis_size'
add_layout_to_subplot <- function(fig,xaxis,yaxis,leg,tot_cond,axis_size,
                                  showlegend=TRUE) {

    fig <- fig %>%
        layout(xaxis = xaxis,yaxis=yaxis,
                legend = list(font = list(size = axis_size-1)),
                annotations = list(x = 0.5 , y = 1.12,
                text = ifelse(tot_cond > 1,leg,""), showarrow = F,
                xref='paper', yref='paper',font = list(size = axis_size)),
                showlegend=showlegend)

  return(fig)

}

plot_plate_info <- function(sample_row,sample_column,sample_type,
                            sample_id,sample_conc_labeled,
                            experiment_name,
                            font_size = 18) {

    sample_type_simple <- sapply(sample_type,sample_type_to_letter)

    sample_matrix    <- data.frame(sample_row,sample_column,sample_type_simple)

    # Set marker colors based on sample_type_simple
    marker_colors <- ifelse(sample_matrix$sample_type_simple == 'S', 'pink', 'rgba(0,0,0,0)')

    fig <- plot_ly(data = sample_matrix,
                x = ~sample_column,
                y = ~sample_row,
                type = 'scatter',
                mode = 'markers+text',
                marker = list(size = 26, color = marker_colors, line = list(color = 'black', width = 1)),
                text = ~sample_type_simple,
                textposition = 'middle center',
                textfont = list(size = font_size-2),
                hoverinfo = 'text',
                hovertext = ~paste(
                    'ID:', sample_id,
                    '<br>Type:', sample_type,
                    '<br>',      sample_conc_labeled)
                )

    fig <- fig %>% layout(xaxis = list(title = '',tickfont = list(size = font_size)),
                          yaxis = list(title = '',tickfont = list(size = font_size),autorange = 'reversed'),
                          title = list(text = paste0("Sample plate layout (",experiment_name,")"),
                                       font = list(size = font_size+1)))

    return(fig)

}

plot_traces <- function(xs,ys,legends,colors,show,
                        font_size = 18,
                        show_grid_x = FALSE,
                        show_grid_y = FALSE,
                        marker_size = 1,
                        line_width = 2,
                        vertical_lines = NULL
                        ) {

    fig <- plot_ly()

    total_traces         <- sum(show)

    max_points_per_trace <- floor(4000 / total_traces)

    for (i in 1:length(xs)) {

        if (!show[i]) next

        x <- unlist(xs[[i]])
        y <- unlist(ys[[i]])

        x <- subset_data(x,max_points = max_points_per_trace)
        y <- subset_data(y,max_points = max_points_per_trace)

        fig <- fig %>%
            add_trace(x = x, y = y,
                      type = 'scatter', mode = 'lines+markers',
                      name = legends[i], color = I(colors[i]),
                      marker = list(size = marker_size),
                      line = list(width = line_width)
            )

    }

    fig <- fig %>% layout(
        xaxis = list(
            title = 'Time (s)',
            tickfont = list(size = font_size),
            titlefont = list(size = font_size),
            tickformat = "digit",
            showgrid = show_grid_x,
            showline = TRUE,
            zeroline = FALSE,
            ticks = "outside",
            tickwidth = 2,
            ticklen = plot_ticks_length_cst
        ),
        yaxis = list(
            title = 'Signal (a.u.)',
            tickfont = list(size = font_size),
            titlefont = list(size = font_size),
            showgrid = show_grid_y,
            showline = TRUE,
            zeroline = FALSE,
            ticks = "outside",
            tickwidth = 2,
            ticklen = plot_ticks_length_cst
        ),
        title = list(
            text = "Sensor traces",
            font = list(size = font_size + 1)
        ),
        legend = list(font = list(size = font_size)),
        margin = list(t = 40)
    )

    return(fig)

}

plot_traces_all <- function(all_xs,all_ys,legends,colors,show,
                            plot_width=12,
                            plot_height=10,
                            plot_type='png',
                            font_size = 18,
                            show_grid_x = FALSE,
                            show_grid_y = FALSE,
                            marker_size = 1,
                            line_width = 2,
                            vertical_lines = NULL) {

    cnt <- 1

    plot_list <- list()

    for (i in 1:length(all_xs)) {

        xs <- all_xs[[i]]
        ys <- all_ys[[i]]

        n <- length(xs)

        legends_tmp <- legends[cnt:(cnt+n-1)]
        colors_tmp  <- colors[cnt:(cnt+n-1)]
        show_tmp    <- show[cnt:(cnt+n-1)]

        cnt <- cnt + n

        fig <- plot_traces(xs,ys,legends_tmp,colors_tmp,show_tmp,
                            font_size = font_size,
                            show_grid_x = show_grid_x,
                            show_grid_y = show_grid_y,
                            marker_size = marker_size,
                            line_width = line_width,
                            vertical_lines = vertical_lines)

        plot_list[[i]] <- fig

    }

    fig <- plot_list_to_fig('Sensor_traces',plot_list,'Sensor traces',font_size,
                            plot_type,plot_width,plot_height,nrows=2,share_axis = FALSE)

    return(fig)

}

plot_steady_state <- function(pyKinetics_fittings,
                              plot_width=12,
                              plot_height=10,
                              plot_type='png',
                              font_size = 18,
                              show_grid_x = FALSE,
                              show_grid_y = FALSE,
                              marker_size = 1,
                              line_width = 2,
                              plot_fit   = FALSE) {


    fig_lst <- list()

    for (fit in pyKinetics_fittings) {

        fig <- plot_ly()

        yFit   <- fit$signal_ss_fit

        names <- fit$names

        for (i in 1:length(fit$lig_conc_lst_per_id)) {

            x <- as.numeric(fit$lig_conc_lst_per_id[[i]])
            y <- fit$signal_ss[[i]]
            fig <- add_trace(fig,x = x,y = y,type = 'scatter',name=names[i],
            marker = list(size = marker_size))

            # plot the fitted signal as lines
            if (!is.null(yFit) && plot_fit) {

                df_fit <- data.frame('x'=x,'y'=yFit[[i]])
                # Sort by x values
                df_fit <- df_fit[order(df_fit$x),]

                fig <- add_trace(
                    fig,data = df_fit,x = ~x,y = ~y,
                    mode = 'lines',name=names[i],
                    line = list(color = 'black',width = line_width),
                    showlegend=FALSE
                    )

            }

            fig <- fig %>% layout(
                xaxis = list(
                    title = 'Ligand concentration (μM)',
                    tickfont = list(size = font_size),
                    titlefont = list(size = font_size),
                    exponentformat = 'power',
                    showgrid = show_grid_x,
                    showline = TRUE,
                    zeroline = FALSE,
                    ticks = "outside",
                    tickwidth = 2,
                    ticklen = plot_ticks_length_cst,
                    type = 'log'
                ),
                yaxis = list(
                    title = 'Signal (a.u.)',
                    tickfont = list(size = font_size),
                    titlefont = list(size = font_size),
                    showgrid = show_grid_y,
                    showline = TRUE,
                    zeroline = FALSE,
                    ticks = "outside",
                    tickwidth = 2,
                    ticklen = plot_ticks_length_cst
                ),
                legend = list(font = list(size = font_size)),
                margin = list(t = 40)
            )

        }

        fig_lst[[length(fig_lst)+1]] <- fig

    }

    fig <- plot_list_to_fig('Steady-state_plot',fig_lst,'Steady-state',font_size,plot_type,plot_width,plot_height,2,FALSE)

    return(fig)

}

plot_steady_state_residuals <- function(
        pyKinetics_fittings,
        plot_width=12,
        plot_height=10,
        plot_type='png',
        font_size = 18,
        show_grid_x = FALSE,
        show_grid_y = FALSE,
        marker_size = 1,
        line_width = 2,
        plot_fit   = FALSE) {

    fig_lst <- list()

    for (fit in pyKinetics_fittings) {

        fig <- plot_ly()

        yFit <- fit$signal_ss_fit

        names <- fit$names

        for (i in 1:length(fit$lig_conc_lst_per_id)) {

            # plot the fitted signal as lines
            if (!is.null(yFit) && plot_fit) {

                x <- as.numeric(fit$lig_conc_lst_per_id[[i]])
                y <- unlist(fit$signal_ss[[i]])

                y_fit <- unlist(yFit[[i]])

                y_res <- y_fit - y

                fig <- add_trace(fig,x = x,y = y_res,
                type = 'scatter',marker = list(size = marker_size),
                name = names[i],showlegend=TRUE)

            }

            fig <- fig %>% layout(
                xaxis = list(
                    title = 'Ligand concentration (μM)',
                    tickfont = list(size = font_size),
                    titlefont = list(size = font_size),
                    exponentformat = 'power',
                    showgrid = show_grid_x,
                    showline = TRUE,
                    zeroline = FALSE,
                    ticks = "outside",
                    tickwidth = 2,
                    ticklen = plot_ticks_length_cst,
                    type = 'log'
                ),
                yaxis = list(
                    title = 'Yfit - Yexp',
                    tickfont = list(size = font_size),
                    titlefont = list(size = font_size),
                    showgrid = show_grid_y,
                    showline = TRUE,
                    zeroline = FALSE,
                    ticks = "outside",
                    tickwidth = 2,
                    ticklen = plot_ticks_length_cst
                ),
                legend = list(font = list(size = font_size)),
                margin = list(t = 40)
            )

        }

        fig_lst[[length(fig_lst)+1]] <- fig

    }

    fig <- plot_list_to_fig('Steady-state_residuals_plot',fig_lst,'Steady-state residuals',
    font_size,plot_type,plot_width,plot_height,2,FALSE)

    return(fig)

}

plot_association_dissociation <- function(pyKinetics_fittings,
                                          plot_width=12,
                                          plot_height=10,
                                          plot_type='png',
                                          font_size = 18,
                                          show_grid_x = FALSE,
                                          show_grid_y = FALSE,
                                          marker_size = 1,
                                          line_width = 2,
                                          plot_assoc = TRUE,
                                          plot_disso = TRUE,
                                          plot_fit   = FALSE,
                                          split_by_smax_id = TRUE,
                                          max_points_per_plot= 2000,
                                          smooth_curves_fit= FALSE,
                                          rolling_window=0.1) {

    fig_lst <- list()

    all_lig_conc <- unlist(lapply(pyKinetics_fittings,function(fit) fit$lig_conc_lst))

    min_lig <- min(all_lig_conc)
    max_lig <- max(all_lig_conc)

    need_color_bar <- min_lig != max_lig

    count <- 0

    for (fit in pyKinetics_fittings) {

        if (!split_by_smax_id) fig <- plot_ly() # Create a new figure for all replicates

        lig_conc <- fit$lig_conc_lst_per_id

        fitted_curves_assoc         <- fit$signal_assoc_fit
        fitted_curves_disso         <- fit$signal_disso_fit

        we_have_fitted_curves1 <- !is.null(fitted_curves_assoc)
        we_have_fitted_curves2 <- !is.null(fitted_curves_disso)

        fc_counter  <- 1

        time_assoc       <- fit$time_assoc_lst
        time_disso       <- fit$time_disso_lst
        raw_curves_assoc <- fit$assoc_lst
        raw_curves_disso <- fit$disso_lst

        for (i in 1:length(lig_conc)) {

            if (split_by_smax_id) fig <- plot_ly() # Create a new figure per replicate

            count <- count + 1

            l   <- as.numeric(lig_conc[[i]])

            n_traces <- length(l)

            max_points_per_trace <- floor(max_points_per_plot / n_traces)

            for (j in 1:n_traces) {

                l_conc    <- l[j]
                hex_color <- get_colors_from_numeric_values(l_conc,min_lig,max_lig)

                if (plot_assoc) {

                    x <- time_assoc[[fc_counter]]
                    y <- raw_curves_assoc[[fc_counter]]

                    # Apply median smoothing
                    if (smooth_curves_fit) {
                        y   <- median_filter(y,x,rolling_window)
                    }

                    x_temp <- subset_data(x,max_points = max_points_per_trace)
                    y_temp <- subset_data(y,max_points = max_points_per_trace)

                    # Use markers for non-smoothed data
                    if (!smooth_curves_fit) {

                        fig <- add_trace(fig,x = x_temp,y = y_temp,type = 'scatter',mode = 'markers',
                        name=paste0(i,' ',l[j]),color=I(hex_color),showlegend=FALSE,
                        marker = list(size = marker_size))

                    } else {
                        # use lines for smoothed data
                        fig <- add_trace(fig,x = x_temp,y = y_temp,type = 'scatter',mode = 'lines',
                        name=paste0(i,' ',l[j]),color=I(hex_color),showlegend=FALSE,
                        line=list(width = line_width))

                    }

                    if (we_have_fitted_curves1 && plot_fit) {

                        y_temp <- subset_data(fitted_curves_assoc[[fc_counter]],max_points = max_points_per_trace)

                        fig <- add_trace(fig,x = x_temp,y = y_temp,type = 'scatter',mode = 'lines',
                                         name=paste0(i,' ',l[j]),
                                         line=list(color='black',width = line_width),showlegend=FALSE,
                                         inherit = FALSE)

                    }

                }

                if (plot_disso) {

                    x <- time_disso[[fc_counter]]
                    y <- raw_curves_disso[[fc_counter]]

                    # Apply median smoothing
                    if (smooth_curves_fit) {
                        y   <- median_filter(y,x,rolling_window)
                    }

                    x_temp <- subset_data(x,max_points = max_points_per_trace)
                    y_temp <- subset_data(y,max_points = max_points_per_trace)

                    # Use markers for non-smoothed data
                    if (!smooth_curves_fit) {

                        fig <- add_trace(fig,x = x_temp,y = y_temp,type = 'scatter',mode = 'markers',
                        name=paste0(i,' ',l[j]),color=I(hex_color),showlegend=FALSE,
                        marker = list(size = marker_size))

                    } else {
                        # use lines for smoothed data
                        fig <- add_trace(fig,x = x_temp,y = y_temp,type = 'scatter',mode = 'lines',
                        name=paste0(i,' ',l[j]),color=I(hex_color),showlegend=FALSE,
                        line=list(width = line_width))

                    }

                    if (we_have_fitted_curves2 && plot_fit) {

                        y_temp <- subset_data(fitted_curves_disso[[fc_counter]],max_points = max_points_per_trace)

                        fig <- add_trace(fig,x = x_temp,y = y_temp,type = 'scatter',mode = 'lines',
                                         name=paste0(i,' ',l[j]),
                                         line=list(color='black',width = line_width),showlegend=FALSE,
                                         inherit = FALSE)

                    }

                }

                fc_counter <- fc_counter + 1

            }

            if (count == 1 && need_color_bar) {

                df <- data.frame(x=0,y=0,values = log10(c(min_lig,max_lig)))

                fig <- add_trace(fig,data=df,x = ~x,y = ~y,type = 'scatter',
                                 mode = 'markers',color=~values,
                                 marker = list(size = 0.001),
                                 showlegend=FALSE)

                log10_min_lig <- log10(min_lig)
                log10_max_lig <- log10(max_lig)

                tickvals <- c(log10_min_lig,log10_min_lig + (log10_max_lig - log10_min_lig)/2,log10_max_lig)

                tickvals_pow <- 10^tickvals
                ticktext     <- signif(tickvals_pow,3)

            }

            fig <- fig %>% layout(
                xaxis = list(
                    title = 'Time (seconds)',
                    tickfont = list(size = font_size),
                    titlefont = list(size = font_size),
                    showgrid = show_grid_x,
                    showline = TRUE,
                    zeroline = FALSE,
                    ticks = "outside",
                    tickwidth = 2,
                    ticklen = plot_ticks_length_cst
                ),
                yaxis = list(
                    title = 'Signal (a.u.)',
                    tickfont = list(size = font_size),
                    titlefont = list(size = font_size),
                    showgrid = show_grid_y,
                    showline = TRUE,
                    zeroline = FALSE,
                    ticks = "outside",
                    tickwidth = 2,
                    ticklen = plot_ticks_length_cst
                ),
                legend = list(font = list(size = font_size)),
                margin = list(t = 40)
            )

            if (count == 1 && need_color_bar) {

                fig <- fig %>% colorbar(
                    title = list(text='Ligand \n concentration (μM)',font=list(size=font_size-1)),
                    tickvals = tickvals,  # Ticks from max to min, rounded to two decimal places
                    ticktext = ticktext,  # Use the same tick values as labels
                    tickfont = list(size = font_size-2),  # Font size of the ticks
                    len = 0.6,  # Length of the color bar
                    outlinewidth = 0)

            }

            if (split_by_smax_id) {

                fig_lst[[length(fig_lst)+1]] <- fig

            }

        }

        if (!split_by_smax_id) fig_lst[[length(fig_lst)+1]] <- fig

    }

    fig <- plot_list_to_fig('Association-traces_plot',fig_lst,'Kinetic traces',font_size,plot_type,plot_width,plot_height,2,FALSE)

    return(fig)

}


# To plot the interaction traces of solution-based experiments
plot_interactions <- function(

            pyKinetics_fittings,
            plot_width=12,
            plot_height=10,
            plot_type='png',
            font_size = 18,
            show_grid_x = FALSE,
            show_grid_y = FALSE,
            marker_size = 1,
            line_width = 2,
            plot_fit   = FALSE,
            max_points_per_plot= 2000,
            smooth_curves_fit= FALSE,
            rolling_window=0.1

            ) {

    fig_lst <- list()

    all_lig_conc  <- unlist(lapply(pyKinetics_fittings,function(fit) fit$lig_conc))

    min_lig <- min(all_lig_conc)
    max_lig <- max(all_lig_conc)

    need_color_bar <- min_lig != max_lig

    # Create a vector of 10 different shapes that can be used for the markers inside the plot.
    shapes <- c('circle','square','diamond','cross','x','triangle-up',
                'triangle-down','triangle-left','triangle-right','pentagon')

    count <- 0

    for (fit in pyKinetics_fittings) {

        lig_conc  <- fit$lig_conc
        prot_conc <- fit$prot_conc
        unq_prot  <- unique(prot_conc)

        raw_curves_assoc <- fit$assoc_lst
        time_assoc       <- fit$time_assoc_lst

        n_traces <- length(lig_conc)

        fig <- plot_ly()

        max_points_per_trace <- floor(max_points_per_plot / n_traces)

        for (i in 1:n_traces) {

            l <- as.numeric(lig_conc[i])
            p <- as.numeric(prot_conc[i])
            hex_color <- get_colors_from_numeric_values(l,min_lig,max_lig)

            shape_idx <- which(unq_prot == p)

            shape <- shapes[shape_idx]

            x <- time_assoc[[i]]
            y <- raw_curves_assoc[[i]]

            # Apply median smoothing
            if (smooth_curves_fit) {
                y   <- median_filter(y,x,rolling_window)
            }

            x <- subset_data(x,max_points = max_points_per_trace)
            y <- subset_data(y,max_points = max_points_per_trace)

            fig <- fig %>%
                add_trace(x = x, y = y,
                          type = 'scatter', mode = 'markers',
                          name = paste0('Ligand: ', l, ' μM, Protein: ', p, ' μM'),
                          marker = list(size = marker_size, color = hex_color, symbol = shape),
                          showlegend = FALSE)


            # If we have fitted curves, plot them as lines
            if (plot_fit) {

                y_fit <- fit$signal_assoc_fit[[i]]

                if (!is.null(y_fit)) {

                    y_fit <- subset_data(y_fit,max_points = max_points_per_trace)

                    fig <- fig %>%
                        add_trace(x = x, y = y_fit,
                                  type = 'scatter', mode = 'lines',
                                  name = paste0('Ligand: ', l, ' μM, Protein: ', p, ' μM'),
                                  line = list(color = 'black', width = line_width),
                                  showlegend = FALSE)
                }
            }
        }
        ## Add a legend for the shapes only
        for (i in 1:length(unq_prot)) {

            fig <- fig %>%
                add_trace(x = c(0,0), y = c(0,0),
                          type = 'scatter', mode = 'markers',
                          name = paste0('Protein: ', unq_prot[i], ' μM'),
                          marker = list(size = 0.001, symbol = shapes[i], color = 'black'),
                          showlegend = TRUE)

        }

        count <- count + 1

        if (count == 1 && need_color_bar) {

            df <- data.frame(x=0,y=0,values = log10(c(min_lig,max_lig)))

            fig <- add_trace(fig,data=df,x = ~x,y = ~y,type = 'scatter',
                             mode = 'markers',color=~values,
                             marker = list(size = 0.001),
                             showlegend=FALSE)

            log10_min_lig <- log10(min_lig)
            log10_max_lig <- log10(max_lig)

            tickvals <- c(log10_min_lig,log10_min_lig + (log10_max_lig - log10_min_lig)/2,log10_max_lig)

            tickvals_pow <- 10^tickvals
            ticktext     <- signif(tickvals_pow,3)

        }

        fig <- fig %>% layout(
            xaxis = list(
            title = 'Time (seconds)',
            tickfont = list(size = font_size),
            titlefont = list(size = font_size),
            showgrid = show_grid_x,
            showline = TRUE,
            zeroline = FALSE,
            ticks = "outside",
            tickwidth = 2,
            ticklen = plot_ticks_length_cst
        ),
        yaxis = list(
            title = 'Signal (a.u.)',
            tickfont = list(size = font_size),
            titlefont = list(size = font_size),
            showgrid = show_grid_y,
            showline = TRUE,
            zeroline = FALSE,
            ticks = "outside",
            tickwidth = 2,
            ticklen = plot_ticks_length_cst
        ),
        legend = list(font = list(size = font_size)),
        margin = list(t = 40)
        )

        if (count == 1 && need_color_bar) {

            fig <- fig %>% colorbar(
            title = list(text='Ligand \n concentration (μM)',font=list(size=font_size-1)),
            tickvals = tickvals,  # Ticks from max to min, rounded to two decimal places
            ticktext = ticktext,  # Use the same tick values as labels
            tickfont = list(size = font_size-2),  # Font size of the ticks
            len = 0.6,  # Length of the color bar
            outlinewidth = 0)

        }

        fig_lst[[length(fig_lst)+1]] <- fig

    }

    fig <- plot_list_to_fig('Interaction_traces_plot',fig_lst,'Interaction traces',
                            font_size,plot_type,plot_width,plot_height,2,FALSE)

    return(fig)

}

plot_interactions_residuals <- function(
            pyKinetics_fittings,
            plot_width=12,
            plot_height=10,
            plot_type='png',
            font_size = 18,
            show_grid_x = FALSE,
            show_grid_y = FALSE,
            marker_size = 1,
            max_points_per_plot= 2000
            ) {

    fig_lst <- list()

    all_lig_conc  <- unlist(lapply(pyKinetics_fittings,function(fit) fit$lig_conc))

    min_lig <- min(all_lig_conc)
    max_lig <- max(all_lig_conc)

    # Create a vector of 10 different shapes that can be used for the markers inside the plot.
    shapes <- c('circle','square','diamond','cross','x','triangle-up',
                'triangle-down','triangle-left','triangle-right','pentagon')

    count <- 0

    for (fit in pyKinetics_fittings) {

        lig_conc  <- fit$lig_conc
        prot_conc <- fit$prot_conc
        unq_prot  <- unique(prot_conc)

        raw_curves_assoc <- fit$assoc_lst
        time_assoc       <- fit$time_assoc_lst

        n_traces <- length(lig_conc)

        fig <- plot_ly()

        max_points_per_trace <- floor(max_points_per_plot / n_traces)

        for (i in 1:n_traces) {

            l <- as.numeric(lig_conc[i])
            p <- as.numeric(prot_conc[i])
            hex_color <- get_colors_from_numeric_values(l,min_lig,max_lig)

            shape_idx <- which(unq_prot == p)

            shape <- shapes[shape_idx]

            x <- time_assoc[[i]]
            y <- raw_curves_assoc[[i]]

            x <- subset_data(x,max_points = max_points_per_trace)
            y <- subset_data(y,max_points = max_points_per_trace)

            y_fit <- fit$signal_assoc_fit[[i]]

            if (!is.null(y_fit)) {

                y_fit <- subset_data(y_fit,max_points = max_points_per_trace)
                y_res <- y_fit - y

                fig <- fig %>%
                    add_trace(x = x, y = y_res,
                              type = 'scatter', mode = 'markers',
                              name = paste0('Ligand: ', l, ' μM, Protein: ', p, ' μM'),
                              marker = list(size = marker_size, color = hex_color, symbol = shape),
                              showlegend = FALSE)

            }

        }
        ## Add a legend for the shapes only
        for (i in 1:length(unq_prot)) {

            fig <- fig %>%
                add_trace(x = c(0,0), y = c(0,0),
                          type = 'scatter', mode = 'markers',
                          name = paste0('Protein: ', unq_prot[i], ' μM'),
                          marker = list(size = 0.001, symbol = shapes[i], color = 'black'),
                          showlegend = TRUE)

        }

        count <- count + 1

        if (count == 1) {

            df <- data.frame(x=0,y=0,values = log10(c(min_lig,max_lig)))

            fig <- add_trace(fig,data=df,x = ~x,y = ~y,type = 'scatter',
                             mode = 'markers',color=~values,
                             marker = list(size = 0.001),
                             showlegend=FALSE)

            log10_min_lig <- log10(min_lig)
            log10_max_lig <- log10(max_lig)

            tickvals <- c(log10_min_lig,log10_min_lig + (log10_max_lig - log10_min_lig)/2,log10_max_lig)

            tickvals_pow <- 10^tickvals
            ticktext     <- signif(tickvals_pow,3)

        }

        fig <- fig %>% layout(
            xaxis = list(
            title = 'Time (seconds)',
            tickfont = list(size = font_size),
            titlefont = list(size = font_size),
            showgrid = show_grid_x,
            showline = TRUE,
            zeroline = FALSE,
            ticks = "outside",
            tickwidth = 2,
            ticklen = plot_ticks_length_cst
        ),
        yaxis = list(
            title = 'Signal (a.u.)',
            tickfont = list(size = font_size),
            titlefont = list(size = font_size),
            showgrid = show_grid_y,
            showline = TRUE,
            zeroline = FALSE,
            ticks = "outside",
            tickwidth = 2,
            ticklen = plot_ticks_length_cst
        ),
        legend = list(font = list(size = font_size)),
        margin = list(t = 40)
        )

        if (count == 1) {

            fig <- fig %>% colorbar(
            title = list(text='Ligand \n concentration (μM)',font=list(size=font_size-1)),
            tickvals = tickvals,  # Ticks from max to min, rounded to two decimal places
            ticktext = ticktext,  # Use the same tick values as labels
            tickfont = list(size = font_size-2),  # Font size of the ticks
            len = 0.6,  # Length of the color bar
            outlinewidth = 0)

        }

        fig_lst[[length(fig_lst)+1]] <- fig

    }

    fig <- plot_list_to_fig('Interaction_traces_residuals_plot',fig_lst,'Residuals',
                            font_size,plot_type,plot_width,plot_height,2,FALSE)

    return(fig)

}

plot_association_dissociation_residuals <- function(

            pyKinetics_fittings,
            plot_width=12,
            plot_height=10,
            plot_type='png',
            font_size = 18,
            show_grid_x = FALSE,
            show_grid_y = FALSE,
            marker_size = 1,
            line_width = 2,
            plot_assoc = TRUE,
            plot_disso = TRUE,
            plot_fit   = FALSE,
            split_by_smax_id = TRUE,
            max_points_per_plot= 2000

            ) {

    fig_lst <- list()

    all_lig_conc <- unlist(lapply(pyKinetics_fittings,function(fit) fit$lig_conc_lst_per_id))

    min_lig <- min(all_lig_conc)
    max_lig <- max(all_lig_conc)

    count <- 0

    for (fit in pyKinetics_fittings) {

        if (!split_by_smax_id) fig <- plot_ly() # Create a new figure for all replicates

        lig_conc <- fit$lig_conc_lst_per_id

        fitted_curves_assoc         <- fit$signal_assoc_fit
        fitted_curves_disso         <- fit$signal_disso_fit

        we_have_fitted_curves1 <- !is.null(fitted_curves_assoc)
        we_have_fitted_curves2 <- !is.null(fitted_curves_disso)

        raw_curves_assoc <- fit$assoc_lst
        raw_curves_disso <- fit$disso_lst
        time_assoc       <- fit$time_assoc_lst
        time_disso       <- fit$time_disso_lst

        fc_counter  <- 1

        for (i in 1:length(lig_conc)) {

            if (split_by_smax_id) fig <- plot_ly() # Create a new figure per replicate

            count <- count + 1

            l   <- as.numeric(lig_conc[[i]])

            n_traces <- length(l)

            max_points_per_trace <- floor(max_points_per_plot / n_traces)

            for (j in 1:n_traces) {

                l_conc    <- l[j]
                hex_color <- get_colors_from_numeric_values(l_conc,min_lig,max_lig)

                if (plot_assoc) {

                    if (we_have_fitted_curves1 && plot_fit) {

                        x_temp <- subset_data(time_assoc[[fc_counter]],max_points = max_points_per_trace)
                        y_temp <- subset_data(raw_curves_assoc[[fc_counter]],max_points = max_points_per_trace)
                        y_fit <- subset_data(fitted_curves_assoc[[fc_counter]],max_points = max_points_per_trace)

                        y_res <- y_fit - y_temp

                        fig <- add_trace(fig,x = x_temp,y = y_res,type = 'scatter',mode = 'markers',
                                        name=paste0(i,' ',l[j]),color=I(hex_color),
                                        showlegend=FALSE,marker = list(size = marker_size))

                    }

                }

                if (plot_disso) {

                    if (we_have_fitted_curves2 && plot_fit) {

                        x_temp <- subset_data(time_disso[[fc_counter]],max_points = max_points_per_trace)
                        y_temp <- subset_data(raw_curves_disso[[fc_counter]],max_points = max_points_per_trace)
                        y_fit <- subset_data(fitted_curves_disso[[fc_counter]],max_points = max_points_per_trace)

                        y_res <- y_fit - y_temp

                        fig <- add_trace(fig,x = x_temp,y = y_res,type = 'scatter',mode = 'markers',
                                         name=paste0(i,' ',l[j]),color=I(hex_color),showlegend=FALSE,marker = list(size = marker_size))

                    }

                }

                fc_counter <- fc_counter + 1

            }

            if (count == 1) {

                df <- data.frame(x=0,y=0,values = log10(c(min_lig,max_lig)))

                fig <- add_trace(fig,data=df,x = ~x,y = ~y,type = 'scatter',
                                 mode = 'markers',color=~values,
                                 marker = list(size = 0.001),
                                 showlegend=FALSE)

                log10_min_lig <- log10(min_lig)
                log10_max_lig <- log10(max_lig)

                tickvals <- c(log10_min_lig,log10_min_lig + (log10_max_lig - log10_min_lig)/2,log10_max_lig)

                tickvals_pow <- 10^tickvals
                ticktext     <- signif(tickvals_pow,3)

            }

            fig <- fig %>% layout(
                xaxis = list(
                    title = 'Time (seconds)',
                    tickfont = list(size = font_size),
                    titlefont = list(size = font_size),
                    showgrid = show_grid_x,
                    showline = TRUE,
                    zeroline = FALSE,
                    ticks = "outside",
                    tickwidth = 2,
                    ticklen = plot_ticks_length_cst
                ),
                yaxis = list(
                    title = 'Yfit - Yexp',
                    tickfont = list(size = font_size),
                    titlefont = list(size = font_size),
                    showgrid = show_grid_y,
                    showline = TRUE,
                    zeroline = FALSE,
                    ticks = "outside",
                    tickwidth = 2,
                    ticklen = plot_ticks_length_cst
                ),
                legend = list(font = list(size = font_size)),
                margin = list(t = 40)
            )

            if (count == 1) {

                fig <- fig %>% colorbar(
                    title = list(text='Ligand \n concentration (μM)',font=list(size=font_size-1)),
                    tickvals = tickvals,  # Ticks from max to min, rounded to two decimal places
                    ticktext = ticktext,  # Use the same tick values as labels
                    tickfont = list(size = font_size-2),  # Font size of the ticks
                    len = 0.6,  # Length of the color bar
                    outlinewidth = 0)

            }

            if (split_by_smax_id) {

                fig_lst[[length(fig_lst)+1]] <- fig

            }

        }

        if (!split_by_smax_id) fig_lst[[length(fig_lst)+1]] <- fig

    }

    fig <- plot_list_to_fig('residuals_plot',fig_lst,'Residuals',font_size,plot_type,plot_width,plot_height,2,FALSE)

    return(fig)

}

plot_simulation <- function(time_assoc_per_cycle,
                            time_disso_per_cycle,
                            signal_per_cycle_assoc,
                            signal_per_cycle_disso,
                            protein_concs,
                            ligand_concs,
                            plot_width=12,
                            plot_height=10,
                            plot_type='png',
                            font_size = 18,
                            show_grid_x = FALSE,
                            show_grid_y = FALSE,
                            marker_size = 1,
                            is_single_cycle = FALSE) {

    fig_lst <- list()

    min_lig <- min(ligand_concs)
    max_lig <- max(ligand_concs)

    n_pcs <- length(protein_concs)

    n_traces <- n_pcs * length(ligand_concs)

    max_points_per_trace <- floor(8000 / n_traces)

    max_signal <- max(unlist(signal_per_cycle_assoc))

    for (i in 1:n_pcs) {

        fig <- plot_ly()

        for (cycle in 1:length(signal_per_cycle_assoc)) {

            time_assoc <- time_assoc_per_cycle[[cycle]]
            time_disso <- time_disso_per_cycle[[cycle]]

            curves_assoc  <- signal_per_cycle_assoc[[cycle]][[i]]

            if (!is.null(time_disso)) {
                curves_disso  <- signal_per_cycle_disso[[cycle]][[i]]
            }

            for (ii in 1:length(curves_assoc)) {

                if (is_single_cycle) {

                    l_conc    <- ligand_concs[cycle]

                } else {

                    l_conc    <- ligand_concs[ii]

                }

                hex_color <- get_colors_from_numeric_values(l_conc,min_lig,max_lig)

                y   <- curves_assoc[[ii]]

                fig <- add_trace(fig,x = as.numeric(time_assoc),y = y,type = 'scatter',mode = 'markers',
                                name=paste0(l_conc),color=I(hex_color),
                                showlegend=FALSE,marker = list(size = marker_size))

                if (!is.null(time_disso)) {

                    y   <- curves_disso[[ii]]

                    fig <- add_trace(fig,x = as.numeric(time_disso),y = y,type = 'scatter',mode = 'markers',
                                    name=paste0(l_conc),color=I(hex_color),
                                    showlegend=FALSE,marker = list(size = marker_size))

                }

            }

        }

        if (i == 1) {

            df <- data.frame(x=0,y=0,values = log10(c(min_lig,max_lig)))

            fig <- add_trace(fig,data=df,x = ~x,y = ~y,type = 'scatter',
                             mode = 'markers',color=~values,
                             marker = list(size = 0.001),
                             showlegend=FALSE)

            log10_min_lig <- log10(min_lig)
            log10_max_lig <- log10(max_lig)

            tickvals <- c(log10_min_lig,log10_min_lig + (log10_max_lig - log10_min_lig)/2,log10_max_lig)

            tickvals_pow <- 10^tickvals
            ticktext     <- signif(tickvals_pow,3)

        }

        leg <- paste0(signif(protein_concs[i],4))

        xaxis <- list(title = 'Time (seconds)',
                    tickfont = list(size = font_size),
                    titlefont = list(size = font_size),
                    showgrid = show_grid_x,
                    showline = TRUE,
                    zeroline = FALSE,
                    ticks = "outside",
                    tickwidth = 2,
                    ticklen = plot_ticks_length_cst
                    )

        yaxis = list(
                    title = 'Signal (a.u.)',
                    tickfont = list(size = font_size),
                    titlefont = list(size = font_size),
                    showgrid = show_grid_y,
                    showline = TRUE,
                    zeroline = FALSE,
                    ticks = "outside",
                    tickwidth = 2,
                    ticklen = plot_ticks_length_cst,
                    range = c(0,max_signal*1.05)
                    )

        fig <- add_layout_to_subplot(fig,xaxis,yaxis,leg,n_pcs,font_size)

        if (i == 1) {

            fig <- fig %>% colorbar(
                title = list(text='Ligand \n concentration (μM)',font=list(size=font_size-1)),
                tickvals = tickvals,  # Ticks from max to min, rounded to two decimal places
                ticktext = ticktext,  # Use the same tick values as labels
                tickfont = list(size = font_size-2),  # Font size of the ticks
                len = 0.6,  # Length of the color bar
                outlinewidth = 0)

        }

        fig_lst[[length(fig_lst)+1]] <- fig

    }

    fig <- plot_list_to_fig('kingenie_simulation_plot',fig_lst,'Kinetic traces',
    font_size,plot_type,plot_width,plot_height,2,TRUE)

    return(fig)

}

# Plot the dominant relaxation rate versus total ligand concentration
# k_obs_dominant_per_pc, k_obs_non_dominant_per_pc, and ligand_concs_per_pc are lists of vectors
# Each vector corresponds to a different protein concentration
plot_relaxation_rates <- function(k_obs_dominant_per_pc,
                                  k_obs_non_dominant_per_pc,
                                  protein_concs,
                                  ligand_concs_per_pc,
                                  plot_width=12,
                                  plot_height=10,
                                  plot_type='png',
                                  font_size = 18,
                                  show_grid_x = FALSE,
                                  show_grid_y = FALSE,
                                  marker_size = 1) {

    fig <- plot_ly()

    n_pcs <- length(protein_concs)

    # Obtain a list of n-colors for the protein concentrations
    colors <- get_colors_from_numeric_values(protein_concs,min(protein_concs),max(protein_concs))

    # Find if ligand_concs_per_pc is a list or a vector
    # If it is a vector, convert it to a list by replicating it n_pcs times
    # We do this to ensure that each protein concentration has its own ligand concentration vector
    if (!is.list(ligand_concs_per_pc)) {
        ligand_concs_per_pc <- rep(list(ligand_concs_per_pc), n_pcs)
    }

    for (i in 1:n_pcs) {

        k_obs_dominant     <- k_obs_dominant_per_pc[[i]]
        k_obs_non_dominant <- k_obs_non_dominant_per_pc[[i]]
        ligand_concs       <- ligand_concs_per_pc[[i]]

        hex_color <- colors[i]

        fig <- add_trace(fig,x = ligand_concs,y = k_obs_dominant,type = 'scatter',mode = 'markers',
                        name=paste0('[P]<sub>0</sub> = ',signif(protein_concs[i],4), ' (μM) ; k<sub>1</sub> '),color = I(hex_color),
                        showlegend=TRUE,marker = list(size = marker_size)
                        )

        fig <- add_trace(fig,x = ligand_concs,y = k_obs_non_dominant,type = 'scatter',mode = 'markers',
                        name=paste0('[P]<sub>0</sub> = ',signif(protein_concs[i],4), ' (μM) ; k<sub>2</sub> '),color = I(hex_color),
                        showlegend=TRUE,marker = list(size = marker_size,symbol="diamond")
                        )

    }

    xaxis <- list(title = '[L]<sub>0</sub> (μM)',
            tickfont = list(size = font_size),
            titlefont = list(size = font_size),
            showgrid = show_grid_x,
            showline = TRUE,
            zeroline = FALSE,
            ticks = "outside",
            tickwidth = 2,
            ticklen = plot_ticks_length_cst
    )

    yaxis = list(
                title = 'Relaxation rate (s<sup>-1</sup>)',
                tickfont = list(size = font_size),
                titlefont = list(size = font_size),
                showgrid = show_grid_y,
                showline = TRUE,
                zeroline = FALSE,
                ticks = "outside",
                tickwidth = 2,
                ticklen = plot_ticks_length_cst
    )

    fig <- add_layout_to_subplot(fig,xaxis,yaxis,"Dominant relaxation rates",1,font_size)

    fig <- config_fig(fig,"Dominant_relaxation_rates",plot_type,plot_width,plot_height)

    return(fig)

}

plot_many_relaxation_rates <- function(
    k_obs_dominant_per_pc_lst,
    k_obs_non_dominant_per_pc_lst,
    protein_concs_lst,
    ligand_concs_per_pc_lst,
    plot_width=12,
    plot_height=10,
    plot_type='png',
    font_size = 18,
    show_grid_x = FALSE,
    show_grid_y = FALSE,
    marker_size = 1) {

    fig_lst <- list()

    n_plots <- length(k_obs_dominant_per_pc_lst)

    for (i in 1:n_plots) {

        k_obs_dominant_per_pc     <- k_obs_dominant_per_pc_lst[[i]]
        k_obs_non_dominant_per_pc <- k_obs_non_dominant_per_pc_lst[[i]]
        protein_concs             <- protein_concs_lst[[i]]
        ligand_concs_per_pc       <- ligand_concs_per_pc_lst[[i]]

        fig <- plot_relaxation_rates(k_obs_dominant_per_pc,
                                     k_obs_non_dominant_per_pc,
                                     protein_concs,
                                     ligand_concs_per_pc,
                                     plot_width=plot_width,
                                     plot_height=plot_height,
                                     plot_type=plot_type,
                                     font_size=font_size,
                                     show_grid_x=show_grid_x,
                                     show_grid_y=show_grid_y,
                                     marker_size=marker_size)

        fig_lst[[length(fig_lst)+1]] <- fig

    }

    fig <- plot_list_to_fig('Dominant_relaxation_rates_plot',fig_lst,'Dominant relaxation rates',
                            font_size,plot_type,plot_width,plot_height,2,FALSE)

    return(fig)

}


# Plot the observed constant versus ligand concentration, color depends on the experiment name
# Requires a three column data frame
# The first column is the ligand concentration
# The second column is the observed constant
# The third column is the experiment name
diagnostic_plot <- function(df,
                            plot_width=12,
                            plot_height=10,
                            plot_type='png',
                            font_size = 18,
                            show_grid_x = FALSE,
                            show_grid_y = FALSE,
                            marker_size = 1,
                            line_width = 2,
                            add_fit = TRUE) {

    fig_lst <- list()
    unq_exp <- unique(df[,3])

    n_exps <- length(unq_exp)

    i <- 0
    for (exp in unq_exp) {

        i <- i + 1
        fig <- plot_ly()

        df_exp <- df[df[,3] == exp,]

        x   <- as.numeric(df_exp[,1])
        y   <- as.numeric(df_exp[,2])

        fig <- add_trace(fig,x = x,y = y,type = 'scatter',mode = 'markers',
                        name=exp,color = I('blue'),
                        showlegend=FALSE,marker = list(size = marker_size+1)
                        )

        # Fit a line and add the line to the plot (black line)
        if (add_fit) {

            lm_fit <- lm(y ~ x)
            x_fit <- seq(min(x), max(x), length.out = 2)
            y_fit <- predict(lm_fit, newdata = data.frame(x = x_fit))

            show_fit_legend <- i == 1

            fig <- add_trace(fig,x = x_fit,y = y_fit,type = 'scatter',mode = 'lines',
                            name='Linear fit',color = I('black'),
                            showlegend=show_fit_legend,line = list(width = line_width)
                            )

        }

        xaxis <- list(title = '[L]<sub>0</sub> (μM)',
                tickfont = list(size = font_size),
                titlefont = list(size = font_size),
                showgrid = show_grid_x,
                showline = TRUE,
                zeroline = FALSE,
                ticks = "outside",
                tickwidth = 2,
                ticklen = plot_ticks_length_cst
        )

        y_axis_label <- 'k<sub>obs</sub> (s<sup>-1</sup>)'

        # include the experiment name, if we have more than one
        if (n_exps > 1) {

            y_axis_label <- paste0(y_axis_label, ' (',exp,')')

        }

        yaxis = list(
                    title = y_axis_label,
                    tickfont  = list(size  = font_size),
                    titlefont = list(size = font_size),
                    showgrid = show_grid_y,
                    showline = TRUE,
                    zeroline = FALSE,
                    ticks = "outside",
                    tickwidth = 2,
                    ticklen = plot_ticks_length_cst
        )

        fig <- add_layout_to_subplot(fig,xaxis,yaxis,exp,1,font_size)

        fig_lst[[length(fig_lst)+1]] <- fig

    }

    fig <- plot_list_to_fig('Observed_constants_plot',fig_lst,'Observed constant',font_size,
                            plot_type,plot_width,plot_height,2,FALSE)

    return(fig)

}