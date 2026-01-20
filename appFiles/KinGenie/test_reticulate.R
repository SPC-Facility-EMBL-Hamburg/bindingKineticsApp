library(reticulate)

# Import the Python package
pykingenie <- import("pykingenie")
pyKinetics <- pykingenie$KineticsAnalyzer()

# List of R script files to source
r_scripts <- c('helpers.R','plot_functions.R')

# Source the R helper functions
for (script in r_scripts) source(paste0('server_files/', script))

gt <- pykingenie$KinGenieCsvSolution("Experiment")
gt$read_csv('/home/os/Downloads/simulation_KinGenie_2025-07-15.csv')

pyKinetics$add_experiment(gt, 'Experiment')

pyKinetics$merge_conc_df_solution()
print(pyKinetics$combined_conc_df)

pyKinetics$init_fittings()
pyKinetics$generate_fittings_solution(pyKinetics$combined_conc_df)

model <- "single"  # or "double" or "one_binding_site"
pyKinetics$submit_fitting_solution(model)

# Print the name of the fitting objects
print(pyKinetics$fittings_names)

name = pyKinetics$fittings_names[1]

# Extract the fitted parameters for the first fitting
fitting <- pyKinetics$fittings[[name]]
fitting$fit_one_binding_site()

print(fitting$fit_params_kinetics)