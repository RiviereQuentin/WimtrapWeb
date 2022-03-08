#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom bslib bs_themer
#' @noRd
app_server <- function( input, output, session ) {
  # Your application server logic 
  mod_AdaptativeSelection_server("AdaptativeSelection_ui_1")
}
