source("ui_files/theme.R")
source("ui_files/logo.R")
source("ui_files/busy_indicator.R")

shinyUI(dashboardPage(

    title = paste0(appName),

    dashboardHeader(  title = logo_grey_light, titleWidth = 200), #logo_grey_light is described in logo.R
    dashboardSidebar( collapsed = F,width = 200,

    sidebarMenu(

        style = "white-space: normal;",
        menuItem("1. Import data",  icon = icon("file-circle-plus"),      tabName = "menu_import"),

        menuItem(div(tags$img(src = "test.svg", width="20px"),"2. Analyse"), tabName = "menu_analyse"),

        menuItem("3. Export data",  icon = icon("file-export"),      tabName = "menu_export"),

        menuItem("Simulation",      icon = icon("magnifying-glass-chart"),   tabName = "menu_simulation"),
        menuItem("User guide",      icon = icon("user-astronaut"),   tabName = "menu_user_guide"),
        #menuItem("Tutorial",        icon = icon("book-open"),        tabName = "menu_tutorial"),
        menuItem("About",           icon = icon("circle-info"),      tabName = "menu_about")

    )),

    dashboardBody(theme_grey_light,

    tags$head(
    tags$style(HTML("
        /* place modal left edge at 1/3 of viewport width */
        .modal-dialog {
        position: fixed !important;
        left: 13.3333vw !important;
        top: 15vh !important;
        transform: none !important;
        margin: 0 !important;
        width: 400px !important;
        max-width: calc(100vw - 40px) !important;
        }
        .modal-content {
        border-radius: 12px !important;
        }
        /* enable vertical scrolling inside the modal */
        .modal-body {
        max-height: calc(100vh - 200px) !important;
        overflow-y: auto !important;
        }
    "))
    ),

        tabItems(

            tabItem(

                tabName = "menu_import",

                    fluidRow(

                        column(4,

                            source("ui_files/1ui_load_input_box.R",local = TRUE)$value,
                            source("ui_files/1ui_processing.R",local = TRUE)$value
                        ),

                        column(8,

                            tags$head(
                                tags$style(HTML("
                                    /* Custom left modal */
                                    .modal-dialog {
                                    position: fixed !important;
                                    left: 80px;        /* margin from the left */
                                    top: 80px;         /* margin from the top */
                                    width: 300px;      /* wider modal */
                                    margin: 0;         /* no default Bootstrap margin */
                                    }
                                    .modal-content {
                                    border-radius: 12px;
                                    }
                                "))
                            ),

                            #Custom CSS to increase plot height
                            tags$head(tags$style("
                            #traces{height:600px !important;}
                            #tracesInSolution{height:600px !important;}
                            #tracesAssDiss{height:600px !important;}
                            ")),

                            conditionalPanel('output.surface_based_binding',
                                  # TabBox to plot the CD spectra and the associated voltage
                                tabBox(title = "", width = 12,id = "tabBoxImport",
                                    tabPanel("Traces",plotlyOutput("traces")),
                                    tabPanel("Assoc. and Diss. traces",plotlyOutput("tracesAssDiss")),
                                    tabPanel("Steps",DTOutput("stepsInfo"))
                                )
                            ),

                            conditionalPanel('!output.surface_based_binding',
                                  # TabBox to plot the CD spectra and the associated voltage
                                tabBox(title = "", width = 12,id = "tabBoxImportInSolution",
                                    tabPanel("Traces",plotlyOutput("tracesInSolution"))
                                )
                            )

                        )
                    ),

                    fluidRow(

                        column(12,

                            source("ui_files/1ui_plotting_box.R",local = TRUE)$value

                        )
                    )
            ),

            tabItem(

                tabName = "menu_analyse",

                    fluidRow(

                        column(5,

                            source("ui_files/2i_fit_data_box.R",local = TRUE)$value

                        ),

                        column(7,

                            conditionalPanel('output.surface_based_binding',
                                source("ui_files/2ui_fitting_box.R",local = TRUE)$value
                            ),

                            conditionalPanel('!output.surface_based_binding',
                                source("ui_files/2ui_fitting_box_in_solution.R",local = TRUE)$value
                            )


                        ),

                        column(12,

                            #Custom CSS to increase plot height
                            tags$head(tags$style("
                            #steady_state{height:700px !important;}
                            #tracesAssDissFit{height:700px !important;}
                            #residuals{height:700px !important;}
                            #diagnostic_plot{height:700px !important;}
                            #tracesFitSolution{height:700px !important;}
                            #residualsSolution{height:700px !important;}
                            #kobs_plot{height:700px !important;}
                            ")),

                            conditionalPanel('output.surface_based_binding',

                                  # TabBox to plot the CD spectra and the associated voltage
                                tabBox(title = "", width = 12,id = "tabBoxFit",
                                    tabPanel("Assoc. and diss. traces",plotlyOutput("tracesAssDissFit")),
                                    tabPanel("Fitted params (boundaries)",DTOutput("fittedParamsBoundaries")),
                                    tabPanel("Residuals",plotlyOutput("residuals"))

                                )
                            ),

                            conditionalPanel('!output.surface_based_binding',

                                  # TabBox to plot the CD spectra and the associated voltage
                                tabBox(title = "", width = 12,id = "tabBoxFitInSolution",
                                    tabPanel("Traces",plotlyOutput("tracesFitSolution")),
                                    tabPanel("Residuals",plotlyOutput("residualsSolution"))

                                )
                            ),

                            source("ui_files/2ui_plotting_box.R",local = TRUE)$value
                        )
                    )

            ),

            tabItem(tabName = "menu_export",

                fluidRow(

                    source("ui_files/ui_export_data.R"         , local=T)$value,
                    source("ui_files/ui_export_logbook.R"      , local=T)$value

                )

            ),

            tabItem(

                tabName = "menu_simulation",

                    fluidRow(

                        source("ui_files/ui_model_selection_box_simulate.R",    local=T)$value,
                        source("ui_files/ui_experimental_parameters_simulate.R",local=T)$value,
                        source("ui_files/ui_binding_parameters_simulate.R",     local=T)$value,
                        #source("ui_files/ui_explore_param_simulate.R"      ,local=T)$value,

                            #Custom CSS to increase plot height
                        tags$head(tags$style("
                        #signal_sim_plot{height:600px !important;}
                        ")),

                        tabBox(title = "", width = 12,id = "tabset_sim",
                            #tabPanel("Fraction of occupied binding sites",plotlyOutput("fractionOccupied_sim_plot")),
                            tabPanel("Signal",plotlyOutput("signal_sim_plot"))
                        ),


                        source("ui_files/ui_plot_options_box_sim.R",        local = TRUE)$value
                    )
            ),
      tabItem(tabName = "menu_user_guide", includeHTML("./shiny_docs/user_guide.html")),
      #tabItem(tabName = "menu_tutorial"  , includeHTML("www/docs/tutorial.html"  )),
      tabItem(tabName = "menu_about"     , includeHTML("./shiny_docs/about.html"     ))
            
      ))

))
