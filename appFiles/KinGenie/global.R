packages <- c("shinydashboard","shinycssloaders","rhandsontable","plotly",
              "shinyalert","reticulate","DT","reshape2","tidyverse",
              "colourpicker","shinyjs")

invisible(lapply(packages, library, character.only = TRUE))

appName     <- "KinGenie"
user        <- Sys.info()['user']


# Find if we are in a Mac computer
isMac <- Sys.info()['sysname'] == 'Darwin'

if (isMac) {
  reticulate::use_python("/Users/oburastero/myenv/bin/python", required = TRUE)
  base_dir <- paste0("/Users/",user,"/Desktop/arise/bindingKineticsApp/appFiles/",appName,"/")
  
} else {
  reticulate::use_python(paste0("/home/",user,"/myenv/bin/python"), required = TRUE)
  base_dir <- paste0("/home/",user,"/bindingKineticsApp/appFiles/",appName,"/")
}

# developer path

# path for the docker user
if (user == 'shiny') base_dir <- paste0("/home/shiny/",appName,'/')

colorPalette8 <- c("#E41A1C","#377EB8","#4DAF4A",
                   "#984EA3","#FF7F00","#FFFF33",
                   "#A65628","#F781BF")

myrenderer <- "function(instance, td, row, col, prop, value, cellProperties) {
                Handsontable.renderers.TextRenderer.apply(this, arguments);
                if (instance.params) {
                    hcols = instance.params.col_highlight
                    hcols = hcols instanceof Array ? hcols : [hcols]
                    hrows = instance.params.row_highlight
                    hrows = hrows instanceof Array ? hrows : [hrows]
                    
                    for (i = 0; i < hcols.length; i++) { 
                        if (hcols[i] == col && hrows[i] == row) {
                            td.style.background = instance.getDataAtCell(row, col);
                        }
                    }
                }
  }"     

myrendererBoolean <- "function(instance, td, row, col, prop, value, cellProperties) {
            Handsontable.renderers.CheckboxRenderer.apply(this, arguments);}"
