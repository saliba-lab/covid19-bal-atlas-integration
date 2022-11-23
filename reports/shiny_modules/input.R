#' Load data in a top column - UI
#' 
#' @id Shiny module key
#' 
datasetInput <- function(id, path) {
  # Row layout
  shiny::fluidRow(
    # Column
    shiny::column(
      width = 1, offset = 0
      ,
      shiny::actionButton(
        inputId = shiny::NS(id, "file_load"),
        label   = "Load data", 
        icon = shiny::icon("database")
      )
    )
    ,
    # Column
    shiny::column(
      width = 2, offset = 0
      ,
      shiny::selectInput(shiny::NS(id, "file"),
                         "Select file", 
                         choices = list.files(path, pattern = "h5ad")
                         )
    )
    ,
    # Column
    shiny::column(
      width = 4, offset = 0.5
      ,
      shinycssloaders::withSpinner(
        shiny::verbatimTextOutput(shiny::NS(id, "shape")), 
        type = 5, size = 0.5, proxy.height = "50px"
      )
    )
  )
}

#' Load data in a top column - Server
#' 
#' @id Shiny module key
#' @datapath Path to data
#' @rv Reactive values to store data
#' 
datasetServer <- function(id, path, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    # Load data on click
    shiny::observeEvent(input$file_load, {
      if (is.null(rv[[id]])) {
        rv[[id]] <- read_h5ad(paste0(path, input$file))
      }
    })
    # Show object dimensions
    output$shape <- shiny::renderPrint({
      shiny::req(input$file_load)
      glue::glue("{dim(rv[[id]])[1]} genes across {dim(rv[[id]])[2]} cells")
    })
  })
}
