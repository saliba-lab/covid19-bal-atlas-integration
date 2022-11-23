#' Scatter plots for QC metrics
#' 
#' @id Shiny module key
#' 
qc_global_Input <- function(id, width=12) {
  # Row layout
  shiny::fluidRow(
    # Column
    shiny::column(
      width = 6, offset = 0
      ,
      shiny::plotOutput(shiny::NS(id, "violin"))
      )
    ,
    # Column
    shiny::column(
      width = 6, offset = 0
      ,
      shiny::plotOutput(shiny::NS(id, "scatter"))
      )
    ,
    shiny::column(
      width = 3,
      shiny::uiOutput(shiny::NS(id, "transform"))
      )
    ,
    shiny::column(
      width = 3,
      shiny::uiOutput(shiny::NS(id, "slider_features")),
      shiny::verbatimTextOutput(shiny::NS(id, "quality_table"))
      )
    ,
    shiny::column(
      width = 3,
      shiny::uiOutput(shiny::NS(id, "slider_counts")))
    ,
    shiny::column(
      width = 3,
      shiny::uiOutput(shiny::NS(id, "slider_percent.mt"))
      )
    )
}

#' Scatter plots for QC metrics
#' 
#' @id Shiny module key
#' @rv Reactive values to store data
#' 
qc_global_Server <- function(id, rv, key) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # Add QC metrics
    shiny::observeEvent(rv[[key]], {
      ds <- rv[[key]]
      
      ds$libsize <- Matrix::colSums(ds@assays@data$counts)
      ds$nfeatures <- Matrix::colSums(ds@assays@data$counts > 0)
      index <- which(stringr::str_detect(rownames(ds), "^MT-"))
      ds$percent.mt <- round(
        Matrix::colSums(ds@assays@data$counts[index, ]) / ds$libsize, 3
      ) * 100
      
      rv[[key]]@colData <- ds@colData
    })
    
    # Axis transform
    output$transform <- shiny::renderUI({
      shiny::req(rv[[key]])
      shiny::selectInput(
        shiny::NS(id, "transform"), 
        "Transform axes with log10?", 
        choices = c(FALSE, TRUE)
        )
    })
    
    # Sliders
    rv[["sl"]] <- list()
    
    output$slider_features <- shiny::renderUI({
      shiny::req(rv[[key]])
      ds <- rv[[key]]
      max <- max(ds$nfeatures)
      shiny::sliderInput(
        shiny::NS(id, "slider_features"),
        "Select range of features", 
        min=0, max=max, step=1, value=c(0, max)
        )
    })
    shiny::observe(rv[["sl"]]$f <- input$slider_features)
    
    output$slider_counts <- shiny::renderUI({
      shiny::req(rv[[key]])
      ds <- rv[[key]]
      max <- max(ds$libsize)
      shiny::sliderInput(
        shiny::NS(id, "slider_counts"),
        "Select range of counts", 
        min=0, max=max, step=1, value=c(0, max)
      )
    })
    shiny::observe(rv[["sl"]]$c <- input$slider_counts)
    
    output$slider_percent.mt <- shiny::renderUI({
      shiny::req(rv[[key]])
      ds <- rv[[key]]
      max <- max(ds$percent.mt)
      shiny::sliderInput(
        shiny::NS(id, "slider_percent.mt"),
        "Select range of mitochondrial percentage", 
        min=0, max=max, step=.1, value=c(0, max)
      )
    })
    shiny::observe(rv[["sl"]]$p <- input$slider_percent.mt)
    
    # Violin
    output$violin <- shiny::renderPlot({
      shiny::req(rv[[key]])
      ds <- rv[[key]]@colData
      rv$quality <- c("TRUE"="good", "FALSE"="bad")[as.character(
        dplyr::between(ds$nfeatures, rv$sl$f[1], rv$sl$f[2]) &
          dplyr::between(ds$libsize, rv$sl$c[1], rv$sl$c[2]) &
          dplyr::between(ds$percent.mt, rv$sl$p[1], rv$sl$p[2])
      )]
      
      obs <- rv[[key]]@colData
      
      df <- obs[, c("libsize", "nfeatures", "percent.mt")]
      df <- as.data.frame(df)
       df$quality <- rv$quality
      df <- tidyr::gather(df, "metric", "value", -quality)
      
      if (input$transform) {
        ytransform <- ggplot2::scale_y_continuous(trans = "log10")
      } else {
        ytransform <- NULL
      }
      
      plot <- ggplot2::ggplot(df, ggplot2::aes("", value)) +
        ggplot2::geom_point(
          ggplot2::aes(col = quality), position = "jitter", shape = 21) +
        ggplot2::geom_violin(fill = NA, size = 2) +
        ggplot2::facet_wrap(~metric, scales = "free") +
        ggplot2::theme_classic(20) +
        ggplot2::labs(x = NULL) +
        ytransform
      
      return(plot)
    })
    
    # Scatter
    output$scatter <- shiny::renderPlot({
      shiny::req(rv[[key]])
      obs <- rv[[key]]@colData
      
      df <- obs[, c("libsize", "nfeatures", "percent.mt")]
      df <- as.data.frame(df)
      df$quality <- rv$quality
      df <- tidyr::gather(df, "metric", "value", -libsize, -quality)
      
      if (input$transform) {
        ytransform <- ggplot2::scale_y_continuous(trans = "log10")
        xtransform <- ggplot2::scale_x_continuous(trans = "log10")
      } else {
        ytransform <- NULL
        xtransform <- NULL
      }
      
      plot <- ggplot2::ggplot(df, ggplot2::aes(libsize, value, col=quality)) +
        ggplot2::geom_point(position = "jitter", shape = 21) +
        ggplot2::facet_wrap(~metric, scales = "free") +
        ggplot2::theme_classic(20) +
        ggplot2::labs(x = NULL) +
        ytransform +
        xtransform
      
      return(plot)
    })
    
    output$quality_table <- shiny::renderText({
      shiny::req(rv$quality)
      
      freq <- table(rv$quality)
      paste0(names(freq), ": ", freq)
    })
    
  })
}