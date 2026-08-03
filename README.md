# The KinGenie app

Last time updated: July 2026

## Introduction

KinGenie is an R shiny application to analyse binding kinetics data. The workflow of KinGenie is divide into four steps:

1) Data Importing: Upload CSV files or raw data from Octet or Gator instruments. 

2) Preprocessing: For surface-based binding, the association phase and dissociation phase can be aligned. For solution-based binding, a time cut-off can be used. For both types of experiments, the reference traces can be subtracted.

3) Analysis: Apply a binding model to extract kinetic parameters such as $k_{\text{on}}$ and $k_{\text{off}}$.

4) Export: Export the fitted curves.

The data import, processing and fitting functions are taken from the Python package [pyKinGenie](https://github.com/osvalB/pykingenie)

## Fast installation

To try KinGenie, you can install it via [docker](https://www.docker.com)

``` bash
docker build -t kingenie:latest -f Dockerfile_deployment .
docker run -p 3838:3838 kingenie:latest
```

## Getting started

To run the app locally you will need R and Python.

R package dependencies are managed with [`renv`](https://rstudio.github.io/renv/), with the package versions declared in `renv.lock`. Python packages are managed with `uv`, with the package versions declared in requirements.txt

1) Clone the repository and open R in the project root

``` bash
git clone https://github.com/SPC-Facility-EMBL-Hamburg/bindingKineticsApp
cd bindingKineticsApp
```

2) Restore the R package environment

``` bash
Rscript -e 'renv::restore()'
```

This reads `renv.lock` and installs all required R packages into a project-local library.

3) Create a Python virtual environment with [`uv`](https://docs.astral.sh/uv/)

If you don't have `uv` installed yet, see the [installation instructions](https://docs.astral.sh/uv/getting-started/installation/).

Then create the virtual environment:

``` bash
uv venv /home/$(whoami)/myenv --python 3.12
```

4) Install the required Python packages

```bash
uv pip install --python /home/$(whoami)/myenv/bin/python -r ./appFiles/requirements.txt
```

5) Run KinGenie

``` bash 
cd appFiles/KinGenie
Rscript -e 'shiny::runApp()'
```

### Updating R dependencies

If you add a new R package to the app, install it as usual and then update the lockfile so the change is reproducible for everyone else:

```R
install.packages("newPackageName")
renv::snapshot()
```

Commit the updated `renv.lock` along with your changes.

### Running tests

Tests use [`shinytest2`](https://rstudio.github.io/shinytest2/):

``` bash
cd ./appFiles/KinGenie/
NOT_CRAN=true Rscript -e "shinytest2::test_app()"
```

## Acknowledgments

The KinGenie app is possible thanks to:

R language: R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

R package shiny:   Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2020). shiny: Web Application Framework for R. R package version 1.4.0.2. https://CRAN.R-project.org/package=shiny

R package tidyverse: Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686

R package shinydashboard:   Winston Chang and Barbara Borges Ribeiro (2018). shinydashboard: Create Dashboards with 'Shiny'. R package version 0.7.1. https://CRAN.R-project.org/package=shinydashboard

R package ggplot2:   H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

R package reshape2:   Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.

R package tippy:   John Coene (2018). tippy: Add Tooltips to 'R markdown' Documents or 'Shiny' Apps. R package version 0.0.1. https://CRAN.R-project.org/package=tippy

R package shinyalert:   Pretty Popup Messages (Modals) in 'Shiny'. R package version 1.1. https://CRAN.R-project.org/package=shinyalert

R package plotly:   C. Sievert. Interactive Web-Based Data Visualization with R, plotly, and shiny. Chapman and Hall/CRC Florida, 2020.

R package rhandsontable:   Jonathan Owen (2018). rhandsontable: Interface to the 'Handsontable.js' Library. R package version 0.3.7. https://CRAN.R-project.org/package=rhandsontable

R package shinyjs:   Dean Attali (2020). shinyjs: Easily Improve the User Experience of Your Shiny Apps in Seconds. R package version 1.1. https://CRAN.R-project.org/package=shinyjs

R package reticulate:   Kevin Ushey, JJ Allaire and Yuan Tang (2020). reticulate: Interface to 'Python'. R package version 1.16. https://CRAN.R-project.org/package=reticulate

R package shinycssloaders:   Andras Sali and Dean Attali (2020). shinycssloaders: Add CSS Loading Animations to 'shiny' Outputs. R package version 0.3. https://CRAN.R-project.org/package=shinycssloaders

R package DT: Xie Y, Cheng J, Tan X (2022). DT: A Wrapper of the JavaScript Library 'DataTables'. R package version 0.25, https://CRAN.R-project.org/package=DT.

R package colourpicker: Attali D (2021). _colourpicker: A Colour Picker Tool for Shiny and for Selecting Colours in Plots_. R package version 1.1.1, <https://CRAN.R-project.org/package=colourpicker>.

Python Software Foundation. Python Language Reference, version 3.12. Available at https://www.python.org

Python package numpy: Travis E, Oliphant. A guide to NumPy, USA: Trelgol Publishing, (2006). Stéfan van der Walt, S. Chris Colbert, and Gaël Varoquaux. The NumPy Array: A Structure for Efficient Numerical Computation, Computing in Science & Engineering, 13, 22-30 (2011), DOI:10.1109/MCSE.2011.37

Python package pandas: Wes McKinney. Data Structures for Statistical Computing in Python, Proceedings of the 9th Python in Science Conference, 51-56 (2010)

Python package scipy: Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.
