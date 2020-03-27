
# filter
sc.filter <- function(mat, filter.params) {
  umi.per.barcode <- colSums(mat)
  genes.per.barcode <- colSums(mat > 0)
  
  sobj <- with(filter.params, {
    mat <- mat[ , which(
      umi.per.barcode >= min.umi 
      & umi.per.barcode <= max.umi 
      & genes.per.barcode >= min.genes
      & genes.per.barcode <= max.genes), drop=FALSE ]
    
    cells.per.gene <- rowSums(mat)
    
    mat <- mat[ which(cells.per.gene >= min.cells),, drop=FALSE ]
    
    sobj <- CreateSeuratObject(counts = mat)
    
    if (mito.pattern != "") {
      sobj[["percent.mito"]] <- PercentageFeatureSet(sobj, pattern = mito.pattern)
      
      sobj <- sobj[ , which(sobj@meta.data$percent.mito <= max.mito) ]
    }
    
    if (!is.na(sample.size)) {
      sobj <- sobj[ , sample(1:dim(sobj)[2], sample.size) ]
    }
    
    sobj
  })
  
  return(sobj)
}

# barcodes plot
sc.plot.barcodes <- function(mat, filter.params) {
  with(filter.params, {
    umi.per.barcode <- colSums(mat)
    x <- sort(umi.per.barcode, decreasing = TRUE)
    plot(x, log="xy", type="l", xlab="Barcodes", ylab="UMI counts")
    
    abline(h=c(min.umi, max.umi), lty="dashed", col="darkred")
    
    w <- which(x >= min.umi & x <= max.umi)
    lines(w, x[w], col="green2", lwd=2)
  })
}

# violin plots
sc.plot.violins <- function(sobj) {
  VlnPlot(sobj, features = intersect(c("nCount_RNA", "nFeature_RNA", "percent.mito"), colnames(sobj@meta.data)), pt.size = 0.2)
}




#' UI function for statistics module
sc_filterUI <- function(id) {
  ns <- NS(id)

  summary.ui <- tagList(
    h4("Filtering summary:"),
    tableOutput(ns("table_summary"))
  )

  barcodes.ui <- tagList(
    plotOutput(ns("plot_barcodes")),
    helpText("Plot showing total UMI per barcode. Barcodes selected for further",
             "analysis are displayed in green."))
  
  distributions.ui <- tagList(
    plotOutput(ns("plot_violins")),
    helpText("Violin plots showing distributions of total UMI per cell, and number of detected genes per cell."))
    
    
  panels.ui <- tabsetPanel(type="pills",
                           tabPanel("Barcode Plot", barcodes.ui),
                           tabPanel("Distributions", distributions.ui),
                           tabPanel("Summary", summary.ui))
  
  mito.ui <- tagList(
    helpText("Enter a regular expression to define the subset of mitochondrial. (e.g. ^mt- for mouse or ^MT- for human. Leave blank to skip this analysis."),
    textInput(ns("text_mito_pattern"), label = "Pattern for mitochondrial genes", value = ""),
    textOutput(ns("text_num_mito"))
  )
  
  more.mito <- tagList(
    helpText("Filter cells based on percentage of mitochondrial RNA."),
    numericInput(ns("num_max_mito"), "Max Mitochodrial RNA Ratio", value = 100, min = 0)
  )
  
  options.ui <- 
    tagList(
      h4("Filtering options"),
      helpText("Filter cells based on total UMI count."),
      numericInput(ns("num_min_umi"), "Min UMI", value = 500, min = 1),
      numericInput(ns("num_max_umi"), "Max UMI", value = 1e6, min = 1),
      helpText("Filter cells based on number of genes detected."),
      numericInput(ns("num_min_genes"), "Min Genes", value = 1, min = 1),
      numericInput(ns("num_max_genes"), "Max Genes", value = 1e6, min = 1),
      helpText("Filter genes based on number of cells where it is expressed."),
      numericInput(ns("num_min_cells"), "Min Cells", value = 5, min = 0),
      h4("Mitochodrial RNA Options"),
      #helpText("Check below to filter cells based on percentage of mitochodrial RNA."),
      #checkboxInput(ns("check_mito"), label = "Analyse mitochidrial RNA", value = FALSE),
      #conditionalPanel(condition = paste0("input['", ns("check_mito"), "']", " == true"), mito.ui)
      mito.ui,
      conditionalPanel(condition = paste0("input['", ns("text_mito_pattern"), "']", ' != ""'), more.mito),
      h4("Sampling options"),
      helpText("Randomly sample a subset of cells passing filters. This allows for quicker",
               "exploration, before using the full dataset."),
      checkboxInput(ns("check_sample"), label = "Sample cells", value=FALSE),
      conditionalPanel(condition = paste0("input['", ns("check_sample"), "']", ' == true'), 
                       numericInput(ns("num_sample"), label = "Number of cells", value = 1000))
    )
  
  fluidRow(
    box(options.ui, width = 4),
    box(panels.ui, width = 8)
  )
}

#' Server function for statistics module
#' 
#' @return A dataframe as a reactive value.
sc_filterServer <- function(input, output, session, sessionData) {

  status <- sessionData$status
  
  # observe({
  #   req(sessionData$dataframe())
  #   
  #   df <- sessionData$dataframe()
  #   
  #   max.umi <- max(colSums(df))
  #   max.genes <- max(colSums(df > 0))
  #   
  #   min.umi <- min(colSums(df))
  #   min.genes <- min(colSums(df > 0))
  #   
  #   updateSliderInput(session, "slider_umi", min = min.umi, max = max.umi, value = c(min.umi, max.umi))
  #   updateSliderInput(session, "slider_genes", min = min.genes, max = max.genes, value = c(min.genes, max.genes))
  #   
  #   return(df)
  # })
  
  filter.params <- reactive({
    # when changing filter params, invalidate further analyses
    status$vargenes_ready <- FALSE
    status$pca_ready <- FALSE
    status$clustering_ready <- FALSE
    
    list(
      min.umi = input$num_min_umi,
      max.umi = input$num_max_umi,
      min.genes = input$num_min_genes,
      max.genes = input$num_max_genes,
      min.cells = input$num_min_cells,
      mito.pattern = input$text_mito_pattern,
      max.mito = input$num_max_mito,
      sample.size = ifelse(input$check_sample == TRUE, input$num_sample, NA)
    )
  })
  
  sobj_raw <- reactive({
    df <- sessionData$dataframe()
    
    print("Creating raw Seurat object...")
    
    withProgress(message = 'Creating Seurat object...', {
      CreateSeuratObject(raw.data = df)
    })
  })
  
  sobj_filtered <- reactive({
    sessionData$status$filter_ready <- FALSE
    sessionData$status$normalize_ready <- FALSE
    
    req(sessionData$dataframe(), filter.params())
    
    print("Creating filtered Seurat object...")
    
    withProgress(message = 'Creating Seurat object...', {
      sobj <- sc.filter(sessionData$dataframe(), filter.params())
    })
    
    sessionData$status$filter_ready <- TRUE
    
    return(sobj)
  })
  
    
  output$table_summary <- renderTable({
    df <- sessionData$dataframe()
    sobj <- sobj_filtered()
    
    dfsum <- with(filter.params(), {
      umi.per.barcode <- colSums(df)
      genes.per.barcode <- colSums(df > 0)
      
      snames <- c(
        "Total barcodes",
        "Total genes",
        "Barcodes below min UMI threshold",
        "Barcodes above max UMI threshold",
        "Barcodes below min gene threshold",
        "Barcodes above max gene threshold",
        "Barcodes passing filters",
        "Genes passing filters")
      
      svalues <- c(
        ncol(df),
        nrow(df),
        sum(umi.per.barcode <= min.umi),
        sum(umi.per.barcode >= max.umi),
        sum(genes.per.barcode <= min.genes),
        sum(genes.per.barcode >= max.genes),
        dim(sobj)[2],
        dim(sobj)[1])
      
      data.frame(row.names = snames, Value=svalues)
    })
    
    return(dfsum)
  }, rownames = TRUE, colnames = FALSE)
    
  output$plot_barcodes <- renderPlot({
    sc.plot.barcodes(sessionData$dataframe(), filter.params())
  })
  
  output$plot_violins <- renderPlot({
    sc.plot.violins(sobj_filtered())
  })

  output$text_num_mito <- renderText({
    sobj <- sobj_filtered()
    
    if (input$text_mito_pattern != "") {
      mito.genes <- grep(input$text_mito_pattern, rownames(sobj), value = TRUE)
      return(paste0("Number of mitochondrial genes: ", length(mito.genes)))
    } else {
      return(paste0("Skip."))
    }
  })
  
  sessionData$sobj_filtered <- sobj_filtered
  sessionData$filter.params <- filter.params
  sessionData$filter.fun <- reactive(sc.filter)
  sessionData$barcode.plot.fun <- reactive(sc.plot.barcodes)
  sessionData$violin.plot.fun <- reactive(sc.plot.violins)
  
  return(sessionData)
}


