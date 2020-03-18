#' launches the cellcuratoR shiny app
#'
#' @export launchApp
#'
#' @return shiny application object
#'
#' @example \dontrun {launchApp()}
#'
#' @import shiny


# wrapper for shinyApp()
# adapted from MangoTheCat/shinyAppDemo (https://github.com/MangoTheCat/shinyAppDemo)
launchApp <- function() {
  shinyApp(ui = shinyAppUI,
           server = shinyAppServer)
}
