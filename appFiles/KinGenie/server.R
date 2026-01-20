options(shiny.maxRequestSize=100*1024^2)
options(stringsAsFactors = F)

# Import the Python package
pykingenie <- import("pykingenie")

# List of R script files to source
r_scripts <- c('helpers.R','plot_functions.R')

# Source the R helper functions 
for (script in r_scripts) source(paste0('server_files/', script))

### End of variables to change

function(input, output, session) {

    user  <- Sys.info()['user']
    if (user != 'os') welcomeMessage() # helpers.R

    # To handle the general processing, the unfolding models,
    # the secondary structure calculation, and the custom models
    pyKinetics                <- pykingenie$KineticsAnalyzer()

    source(paste0(base_dir,"reactives/reactives_values.R"                 ), local = T)
    source(paste0(base_dir,"reactives/reactives.R"                        ), local = T)
    source(paste0(base_dir,"reactives/plot_reactives.R"                   ), local = T)
    source(paste0(base_dir,"reactives/table_reactives.R"                  ), local = T)
    source(paste0(base_dir,"reactives/analysis_reactives.R"               ), local = T)
    source(paste0(base_dir,"reactives/analysis_reactives_solution.R"      ), local = T)
    source(paste0(base_dir,"reactives/download_reactives.R"               ), local = T)
    source(paste0(base_dir,"reactives/simulation_reactives.R"             ), local = T)

}

