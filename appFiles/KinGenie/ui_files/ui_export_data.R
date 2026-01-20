box(title = "Data", width = 2, solidHeader = T, status = "primary",
    fluidRow(
        column(12, p(style = "font-size: 120%",HTML(""),downloadLink('download_curves','Raw data (kinetics)'))),
        column(12, p(style = "font-size: 120%",HTML(""),downloadLink('download_fitted_curves','Fitted data (kinetics)'))),
    ))



