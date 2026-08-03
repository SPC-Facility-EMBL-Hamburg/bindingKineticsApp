# This dockerfile is used to test KinGenie in the CI/CD pipeline, 
# For running KinGenie in your computer use Dockerfile_deployment!
# ## From Linux
# $ docker build -t binding-kinetics:latest -f Dockerfile .
## From MacOS
# $ colima start
# $ docker build -t binding-kinetics:latest --file Dockerfile . --no-cache --platform linux/amd64
FROM rocker/r2u:24.04 AS base

LABEL maintainer="oburastero@gmail.com"

COPY appFiles/Rprofile.site /usr/lib/R/etc/Rprofile.site

ARG SYSTEM_PACKAGES="wget python3 python3-pip python3-venv python3-dev libssl-dev gnupg ca-certificates"
ARG R_PACKAGES="reticulate shiny shinyjs shinydashboard shinycssloaders  \
shinyalert tippy rhandsontable plotly reshape2 tidyverse DT colourpicker box testthat shinytest2 chromote"

RUN apt-get update && \
    apt-get install --no-install-recommends --yes ${SYSTEM_PACKAGES} && \
    wget -q https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb && \
    apt-get install -y ./google-chrome-stable_current_amd64.deb && \
    rm google-chrome-stable_current_amd64.deb && \
    install.r ${R_PACKAGES} && \
    apt-get autoclean && \
    rm -rf /var/lib/apt/lists/* && \
    useradd -ms /bin/bash shiny

COPY appFiles/requirements.txt /tmp/requirements.txt

USER shiny

RUN python3 -m venv /home/shiny/myenv && \
    /home/shiny/myenv/bin/pip install --prefer-binary --no-cache-dir -r /tmp/requirements.txt

FROM base AS ci

COPY --chown=shiny:shiny ./appFiles/  /home/shiny/

USER root

EXPOSE 3838

USER shiny
ENV CHROMOTE_CHROME=/usr/bin/google-chrome
ENV CHROMOTE_CHROME_ARGS="--no-sandbox --disable-dev-shm-usage"

CMD ["R", "-e", "shiny::runApp('/home/shiny/KinGenie')"]
