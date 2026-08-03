box::use(
  reticulate[
    use_python
  ]
)

appName <- "KinGenie"
user <- Sys.info()["user"]

# Find if we are in a Mac computer
isMac <- Sys.info()["sysname"] == "Darwin"

py <- Sys.getenv("RETICULATE_PYTHON", "")

if (nzchar(py)) {
  use_python(py, required = TRUE)
} else if (isMac) {
  use_python("/Users/oburastero/myenv/bin/python", required = TRUE)
} else {
  use_python(paste0("/home/", user, "/myenv/bin/python"), required = TRUE)
}