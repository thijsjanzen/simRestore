#' runs shiny app locally
#' @description
#' This function allows for local execution of the shiny app. Alternatively,
#' an online version of this app
#' \href{https://thijsjanzen.shinyapps.io/simRestoreApp/}{can be found here}.
#'
#' @return No return value
#' @export
run_shiny_app <- function() {
  appDir <- system.file("shiny-examples",
                        "simRestoreApp",
                        package = "simRestore")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `simRestore`.",
         call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
