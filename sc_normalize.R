
#' UI function for statistics module
sc_normalizeUI <- function(id) {
  ns <- NS(id)
  
  panels.ui <- tagList(
    plotOutput(ns("plot_gene_expression")),
    helpText("Histogram of mean gene expression across all cells.")
    )
  
  options.ui <- 
    tagList(
      h4("Options"),
      selectInput(ns("sel_method"), label = "Method", 
                  choices = c("log-normalization"), selected = "log-normalization")
    )
  
  fluidRow(
    box(options.ui, width = 4),
    box(panels.ui, width = 8)
  )
}

#' Server function for statistics module
#' 
#' @return A dataframe as a reactive value.
sc_normalizeServer <- function(input, output, session, sessionData) {
  
  sobj_norm <- reactive({
    sobj <- sessionData$sobj_filtered()
    
    print("Normalizing...")
    
    withProgress(message = 'Normalizing data...', {
      NormalizeData(sobj, 
                    normalization.method = "LogNormalize", 
                    scale.factor = median(sobj@meta.data$nUMI))
    })      
  })

  output$plot_gene_expression <- renderPlot({
    sobj <- sobj_norm()

    withProgress(message = 'Making gene expression plot...', {
      exprs <- rowMeans(as.matrix(sobj@data))

      hist(exprs, breaks=100, xlab="Mean normalized expression (log)")
    })
  })
  
  sessionData$sobj_norm <- sobj_norm
  
  return(sessionData)
}


