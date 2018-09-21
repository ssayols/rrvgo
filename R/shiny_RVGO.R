#' shiny_RVGO
#'
#' Launch an interactive web interface
#' @param ... other params sent to shiny::runApp()
shiny_RVGO <- function(...) {
  shiny::runApp(appDir=system.file("shiny_RVGO", package="RVGO"), ...)
}
