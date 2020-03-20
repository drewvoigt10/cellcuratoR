#' launches the cellcuratoR shiny app
#'
#' @export launchApp
#'
#' @return shiny application object
#'
#' @examples
#' \dontrun{
#' launchApp()
#' }
#'
#' @importFrom shiny shinyApp


# wrapper for shinyApp()
# adapted from MangoTheCat/shinyAppDemo (https://github.com/MangoTheCat/shinyAppDemo)
launchApp <- function() {
  shinyApp(ui = shinyAppUI,
           server = shinyAppServer)
}
