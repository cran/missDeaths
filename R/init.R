.onAttach <- function(...) {
  welcome <- paste("",
                   "----------------------------",
                   " Welcome to missDeaths v2.5",
                   "----------------------------",
                   sep = "\n")
  packageStartupMessage(welcome)
}

.onUnload <- function (libpath) {
  library.dynam.unload("missDeaths", libpath)
}