box::use(
  shiny[
    moduleServer,
    reactiveVal
  ]
)

#' @export
plotConfigServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    reactiveVal(list(
      x_label = "Time (s)", # x-axis label
      y_label = "Signal (a.u.)", # y-axis label
      width = 1000, # plot width in pixels
      height = 600, # plot height in pixels
      type = "png", # plot type
      axis_size = 16, # axis size
      x_legend_pos = 1, # x legend position
      y_legend_pos = 1, # y legend position
      color_bar_length = 0.5, # color bar length
      show_colorbar = TRUE, # show colorbar
      show_grid_x = FALSE, # show x grid
      show_grid_y = FALSE, # show y grid
      marker_size = 5, # marker size
      line_width = 2, # line width
      max_points = 4000, # max points to display per subplot
      n_xticks = 4, # number of ticks on x-axis
      n_yticks = 3, # number of ticks on y-axis
      tick_length = 8, # tick length
      tick_width = 2, # tick width
      split_by_smax = TRUE,
      smooth_curves = FALSE,
      rolling_window_size = 0.5
    ))
  })
}
