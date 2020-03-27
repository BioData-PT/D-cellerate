
sc.normalize <- function(sobj, normalize.params) {
  sobj <- NormalizeData(sobj, 
                        normalization.method = normalize.params$method, 
                        scale.factor = normalize.params$scale.factor)
}

sc.normalize.histogram <- function(sobj) {
  exprs <- rowMeans(as.matrix(sobj[["RNA"]]@data))
  
  hist(exprs, breaks=100, xlab="Mean normalized expression (log)")
}

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

  status <- sessionData$status
  
  normalize.params <- reactive({
    # when changing normalization params, invalidate further analyses
    status$vargenes_ready <- FALSE
    status$pca_ready <- FALSE
    status$clustering_ready <- FALSE
    
    sobj <- sessionData$sobj_filtered()
    
    list(
      method = "LogNormalize",
      scale.factor = 1e4
      )
  })
  
  sobj_norm <- reactive({
    sessionData$status$normalize_ready <- FALSE
    
    req(sessionData$sobj_filtered())
    
    print("Normalizing...")

    sobj <- NULL
    
    withProgress(message = 'Normalizing data...', {
      sobj <- sc.normalize(sessionData$sobj_filtered(), normalize.params())
    })      
    
    sessionData$status$normalize_ready <- TRUE
    
    return (sobj)
  })

  output$plot_gene_expression <- renderPlot({
    sc.normalize.histogram(sobj_norm())
  })
  
  sessionData$sobj_norm <- sobj_norm
  sessionData$normalize.params <- normalize.params
  sessionData$normalize.fun <- reactive(sc.normalize)
  sessionData$normalize.histogram.fun <- reactive(sc.normalize.histogram)

    
  return(sessionData)
}


