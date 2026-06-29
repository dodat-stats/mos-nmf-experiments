## This makes sure that R loads the workflowr package
## automatically, everytime the project is loaded
if (Sys.getenv("RSTUDIO_PANDOC") == "") {
  rstudio_pandoc = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64"
  if (file.exists(file.path(rstudio_pandoc, "pandoc"))) {
    Sys.setenv(RSTUDIO_PANDOC = rstudio_pandoc)
  }
}

if (requireNamespace("workflowr", quietly = TRUE)) {
  message("Loading .Rprofile for the current workflowr project")
  library("workflowr")
} else {
  message("workflowr package not installed, please run install.packages(\"workflowr\") to use the workflowr functions")
}
