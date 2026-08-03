box::use(
  . / box / logo[
    get_logo
  ],
  . / box / modules / menu / menu_analyse[
    menuAnalyseUI
  ],
  . / box / modules / menu / menu_export[
    menuExportUI
  ],
  . / box / modules / menu / menu_import[
    menuImportUI
  ],
  . / box / modules / menu / menu_simulation[
    menuSimulationUI
  ],
  . / box / theme[
    theme_grey_light
  ],
  shiny[
    column,
    fluidRow,
    icon,
    includeHTML,
    shinyUI,
    tags
  ],
  shinydashboard[
    dashboardBody,
    dashboardHeader,
    dashboardPage,
    dashboardSidebar,
    menuItem,
    sidebarMenu,
    tabItem,
    tabItems
  ]
)

logo_grey_light <- get_logo(appName)

ui <- shinyUI(dashboardPage(
  title = paste0(appName),
  dashboardHeader(title = logo_grey_light, titleWidth = 200), # logo_grey_light is described in logo.R
  dashboardSidebar(
    collapsed = F, width = 200,
    sidebarMenu(
      style = "white-space: normal;",
      menuItem("1. Import data", icon = icon("file-circle-plus"), tabName = "menu_import"),
      menuItem(div(tags$img(src = "test.svg", width = "20px"), "2. Analyse"), tabName = "menu_analyse"),
      menuItem("3. Export data", icon = icon("file-export"), tabName = "menu_export"),
      menuItem("Simulation", icon = icon("magnifying-glass-chart"), tabName = "menu_simulation"),
      menuItem("User guide", icon = icon("user-astronaut"), tabName = "menu_user_guide"),
      menuItem("About", icon = icon("circle-info"), tabName = "menu_about")
    )
  ),
  dashboardBody(
    theme_grey_light,
    tags$head(
      tags$link(
        rel = "stylesheet",
        type = "text/css",
        href = "css/style.css"
      )
    ),
    tabItems(
      tabItem(
        tabName = "menu_import", menuImportUI("menu_import")
      ),
      tabItem(
        tabName = "menu_analyse",
        menuAnalyseUI("menu_analyse")
      ),
      tabItem(
        tabName = "menu_export",
        menuExportUI("menu_export")
      ),
      tabItem(
        tabName = "menu_simulation",
        menuSimulationUI("menu_simulation")
      ),
      tabItem(tabName = "menu_user_guide", includeHTML("./docs/user_guide.html")),
      tabItem(tabName = "menu_about", includeHTML("./docs/about.html"))
    )
  )
))
