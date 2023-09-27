#' runs shiny version
#' @examples
#' run_shiny_app()
#' @export
run_shiny_app <- function() {
  appDir <- system.file("shiny-examples",
                        "simRestoreApp",
                        package = "simRestore")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.",
         call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
