# DE module


#' UI function for statistics module
sc_markersUI <- function(id) {
  ns <- NS(id)
  
  # cluster.ui <- tagList(
  #   plotOutput(ns("plot_clusters")),
  #   verbatimTextOutput(ns("text_clusters"))
  #   )
  # 
  markers.ui <- tagList(
    plotOutput(ns("plot_top_markers")),
    helpText('Heatmap showing gene expression for the top markers in each cluster.'),
    hr(),
    dataTableOutput(ns('table_markers'))
  )
  
  # viz.ui <- tagList(
  #   plotOutput(ns("plot_tsne_gene")),
  #   inputPanel(
  #     textInput(ns("text_gene_id"), label="Gene ID")
  #   )
  # )

  vizmarkers.ui <- tagList(
    plotOutput(ns("plot_genes")),
    inputPanel(
      selectInput(ns("sel_type"), label = "Plot", choices = c("Feature plots", "Violin plots")),
      selectInput(ns("sel_cluster"), label = "Cluster", choices = NULL),
      selectInput(ns("sel_genes"), label = "Genes", choices = NULL, multiple = TRUE)
    )
  )
  
  panels.ui <- tabsetPanel(type="pills",
                           tabPanel("Markers", markers.ui),
                           tabPanel("Visualize Markers", vizmarkers.ui)
                           )
  
  sidepanel.ui <- tagList(
      h4("Differential Expression Options"),
      helpText("Select options for marker gene discovery."),
      selectInput(ns("sel_test"), label="Test to use", 
                  choices=c("Wilcoxon rank sum"="wilcox", 
                            "Student's t-test"="ttest",
                            "Standard AUC classifier"="roc")),
      numericInput(ns("num_logfc"), label = "LogFC threshold", value=0.25),
      helpText("Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25 Increasing logfc.threshold speeds up the function, but can miss weaker signals."),
      numericInput(ns("num_minpct"), label = "Minimum fraction of cells", value=0.1),
      helpText("Only test genes that are detected in a minimum fraction of cells in either of the two populations.", 
               "Meant to speed up the function by not testing genes that are very infrequently expressed.")
      #h4("More Options"),
      #helpText("Moooooore options! Mooooore!")
    )

  fluidRow(
    box(sidepanel.ui, width = 4),
    box(panels.ui, width = 8)
  )
  
}

#' Server function for statistics module
#' 
#' @return A dataframe as a reactive value.
sc_markersServer <- function(input, output, session, sessionData) {
  
  # clusters.observer <- observe({
  #   #cluster.names <- sessionData$cluster_names()
  #   
  #   #print("YA!")
  #   #print(cluster.names)
  #   
  #   
  #   suspend()
  # }, suspended = TRUE)
  
  all_markers <- reactive({
    sobj <- sessionData$sobj_cluster()
    
    logfc.threshold <- input$num_logfc
    min.pct <- input$num_minpct
    
    clusters <- as.character(sort(unique(sobj@ident)))
    n <- length(clusters)
    
    print("Finding markers for all clusters...")
    
    withProgress(message = 'Finding markers for all clusters...', value = 0, max = n+1, {
      #markers <- FindAllMarkers(object = sobj, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
      markers <- lapply(clusters, function(cl) {
        incProgress(1, detail = paste("Testing cluster", cl))
        
        df <- FindMarkers(sobj, ident.1 = cl, only.pos = TRUE, 
                          min.pct = min.pct, logfc.threshold = logfc.threshold, 
                          thresh.use = 0.25)
        df$cluster <- cl
        df$gene <- rownames(df)
        
        return(df)
      })
      names(markers) <- clusters
      
    })
    
    updateSelectInput(session, "sel_cluster", label = "Cluster", choices = clusters, selected = clusters[1])
    
    return(markers)
  })
  
  # available.genes <- reactive({
  #   markers <- all_markers()
  #   markers <- markers[[ input$sel_cluster ]]
  # 
  #   genes <- markers$gene
  #   
  #   updateSelectInput(session, "sel_genes", choices=genes, selected=NULL)
  #   
  #   return(available.genes)
  # })
  
  observeEvent(input$sel_cluster, {
    ready <- sessionData$clustering_stats$ready
    
    if (ready == TRUE) {
      print("updating genes...")
      
      markers <- all_markers()
      markers <- markers[[ input$sel_cluster ]]
      
      genes <- markers$gene
      
      updateSelectInput(session, "sel_genes", choices=genes, selected=genes[1])
    }
  })
  
  output$plot_top_markers <- renderPlot({
    sobj <- sessionData$sobj_tsne_cluster()
    markers <- all_markers()
    
    top.markers <- do.call(rbind, lapply(markers, head))
    DoHeatmap(sobj, genes.use = top.markers$gene, slim.col.label = TRUE, remove.key = TRUE)
  })
  
  output$table_markers <- renderDataTable({
    markers <- all_markers()
    
    df <- do.call(rbind, markers)
    
    #print(df)
    
    return(DT::datatable(df, options = list(scrollX = TRUE), filter = "top"))
  })
  
  # output$plot_tsne_gene <- renderPlot({
  #   req(input$sel_cluster)
  #   
  #   markers <- all_markers()
  #   markers <- markers[[ input$sel_cluster ]]
  #   
  #   sobj <- sessionData$sobj_tsne_cluster()
  #   
  #   sobj@meta.data$markers.mean <- colMeans(sobj@data[ rownames(markers), ])
  #   
  #   FeaturePlot(sobj, features.plot = "markers.mean", no.legend = FALSE)
  # })

  output$plot_genes <- renderPlot({
    req(input$sel_type)
    req(input$sel_cluster)
    req(input$sel_genes)
    
    all_markers()
    #markers <- markers[[ input$sel_cluster ]]
    
    genes <- input$sel_genes
    
    sobj <- sessionData$sobj_tsne_cluster()
    
    #sobj@meta.data$markers.mean <- colMeans(sobj@data[ markers$gene, ])
    
    if (input$sel_type == "Feature plots") {
      FeaturePlot(sobj, features.plot = genes, no.legend = FALSE)
    } else { 
      VlnPlot(sobj, features.plot = genes, point.size.use=0.5)
    }
  })
  
  sessionData$all_markers <- all_markers
  
  return(sessionData)
}


