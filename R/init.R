.onAttach <- function(...) {
  welcome <- paste("",
                   "----------------------------",
                   " Welcome to missDeaths v2.8",
                   "----------------------------",
                   sep = "\n")
  packageStartupMessage(welcome)
}

.onUnload <- function (libpath) {
  library.dynam.unload("missDeaths", libpath)
}