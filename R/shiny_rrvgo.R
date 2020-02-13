#' shiny_rrvgo
#' Launch an interactive web interface.
#'
#' @param ... other params sent to shiny::runApp().
shiny_rrvgo <- function(...) {
  shiny::runApp(appDir=system.file("shiny_rrvgo", package="rrvgo"), ...)
}
