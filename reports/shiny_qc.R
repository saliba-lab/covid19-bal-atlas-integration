#
# Shiny app for quality control
#

# Source functions
suppressMessages({
  source("~/local/covid19-bal-atlas-integration/bin/_functions.R")
  modules <- "~/local/covid19-bal-atlas-integration/reports/shiny_modules/"
  for (i in paste0(modules, list.files(modules))) {
    source(i)
  }
})

# Variables
file_path <- "~/local/covid19-bal-atlas-integration/data/"

# Define UI
ui <- shiny::fluidPage(
  
  shiny::navbarPage(
    title = "Quality Control"
    ,
    shiny::tabPanel(
      title = "Overview"
      ,
      datasetInput("file", file_path)
      ,
      qc_global_Input("file_qc_global")
      )
    ,
    shiny::tabPanel("Samples")
    )
  )

# Define server
server <- function(input, output, session) {
  
  # Close session when app closed
  shiny::onSessionEnded(
    session = session,
    fun     = stopApp
  )
  
  # Create data storage
  rv <- shiny::reactiveValues()
  
  datasetServer("file", file_path, rv)
  
  qc_global_Server("file_qc_global", rv, "file")
  
}

# Run the application 
shiny::shinyApp(ui = ui, server = server)
