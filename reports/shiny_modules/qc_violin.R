#' Violin plots for QC metrics
#' 
#' @id Shiny module key
#' 
qcviolinInput <- function(id, width=12) {
  # Row layout
  shiny::fluidRow(
    # Column
    shiny::column(
      width = width, offset = 0
      ,
      shiny::plotOutput(shiny::NS(id, "plot"))
      ,
      shiny::selectInput(shiny::NS(id, "ytransform"),
                         "Transform y-axis with log10?",
                         choices = c(FALSE, TRUE)
                         )
    )
  )
}

#' Violin plots for QC metrics
#' 
#' @id Shiny module key
#' @rv Reactive values to store data
#' 
qcviolinServer <- function(id, rv, key) {
  shiny::moduleServer(id, function(input, output, session) {
    
    output$plot <- shiny::renderPlot({
      shiny::req(rv[[key]])
      obs <- rv[[key]]@colData
      
      df <- obs[, c("libsize", "nfeatures", "percent.mt")]
      df <- as.data.frame(df)
      df <- tidyr::gather(df, "metric", "value")
      
      if (input$ytransform) {
        ytransform <- ggplot2::scale_y_continuous(trans = "log10")
      } else {
        ytransform <- NULL
      }
      
      ggplot2::ggplot(df, ggplot2::aes("", value)) +
        ggplot2::geom_point(position = "jitter", shape = 21) +
        ggplot2::geom_violin() +
        ggplot2::facet_wrap(~metric, scales = "free") +
        ggplot2::theme_classic(20) +
        ggplot2::labs(x = NULL) +
        ytransform
    })
    
  })
}