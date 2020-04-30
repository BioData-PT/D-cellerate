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
  
  de.panel <- tagList(
    selectInput(ns("sel_de_cluster1"), label = "Cluster 1", choices = NULL),
    selectInput(ns("sel_de_cluster2"), label = "Cluster 2", choices = NULL),
    dataTableOutput(ns('table_de')),
    plotOutput(ns("plot_de_genes"))
  )
  
  panels.ui <- tabsetPanel(type="pills",
                           tabPanel("Differential expression", de.panel),
                           tabPanel("Markers", markers.ui),
                           tabPanel("Visualize Markers", vizmarkers.ui)
                           )
  
  sidepanel.ui <- tagList(
      h4("Differential Expression Options"),
      helpText("Select options for marker gene discovery."),
      selectInput(ns("sel_test"), label="Test to use", 
                  choices=c("Wilcoxon rank sum"="wilcox", 
                            "Student's t-test"="t", #aggs
                            "Standard AUC classifier"="roc")),
      numericInput(ns("num_logfc"), label = "LogFC threshold", value=0.25),
      helpText("Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25 Increasing logfc.threshold speeds up the function, but can miss weaker signals."),
      numericInput(ns("num_minpct"), label = "Minimum fraction of cells", value=0.25),
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
  
  marker_list <- reactiveValues()
  cluster_list <- reactiveValues()

  all_markers <- reactive({
    sobj <- sessionData$sobj_cluster()
    
    logfc.threshold <- input$num_logfc
    min.pct <- input$num_minpct
    
    allclusters <- split(Idents(sobj), Idents(sobj))

    new_markers <- reactiveValues()
    
    # determine which clusters have changed    
    for (cl in names(allclusters)) {
      saved.cells <- as.character(cluster_list[[ cl ]])
      new.cells <- as.character(allclusters[[ cl ]])
      
      if (length(saved.cells) == length(new.cells)) {
        new_markers[[ cl ]] <- marker_list[[ cl ]]
      } 
    }
    
    cluster_list <<- do.call(reactiveValues, allclusters)
    
    print("Finding markers...")
    
    n <- length(allclusters)
    clusters <- names(allclusters)
    
    withProgress(message = 'Finding markers...', value = 0, max = n+1, {
      for (cl in names(allclusters)) {
        incProgress(1, detail = paste("Testing cluster", cl))
        
        if (is.null(new_markers[[ cl ]])) {
          df <- FindMarkers(sobj, ident.1 = cl, only.pos = TRUE, 
                            min.pct = min.pct, logfc.threshold = logfc.threshold, 
                            test.use = input$sel_test)
          df <- df[ df$p_val_adj < 0.01, ]
          df$cluster <- cl
          df$gene <- rownames(df)
          
          new_markers[[ cl ]] <- df
        }
      }
    })

    marker_list <<- new_markers
    
    updateSelectInput(session, "sel_cluster", label = "Cluster", choices = clusters, selected = clusters[1])

    return(reactiveValuesToList(new_markers))
  })

  observeEvent(input$sel_cluster, {
    ready <- sessionData$status$clustering_ready
    
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
    DoHeatmap(sobj, features = top.markers$gene) 
  })
  
  output$table_markers <- renderDataTable({
    markers <- all_markers()
    
    df <- do.call(rbind, markers)
    
    #print(df)
    
    return(DT::datatable(df, options = list(scrollX = TRUE), filter = "top"))
  })
  
  output$plot_genes <- renderPlot({
    req(input$sel_type)
    req(input$sel_cluster)
    req(input$sel_genes)
    
    all_markers()

    genes <- input$sel_genes
    
    sobj <- sessionData$sobj_tsne_cluster()
    
    if (input$sel_type == "Feature plots") {
      FeaturePlot(sobj, features = genes)
    } else { 
      VlnPlot(sobj, features = genes, pt.size=0.5)
    }
  })
  
  #
  # Differential expression
  #
  
  observe({
    if (sessionData$status$clustering_ready) {
      sobj <- sessionData$sobj_cluster()
      
      clusters <- sort(unique(Idents(sobj)))
      
      updateSelectInput(session, "sel_de_cluster1", choices = clusters, selected = clusters[1])
      updateSelectInput(session, "sel_de_cluster2", choices = clusters, selected = clusters[2])
    }
  })
  
  de_result <- reactive({
    
    #sessionData$status$marker_ready <- FALSE  #aggs 
    
    sobj <- sessionData$sobj_cluster()

    req(input$sel_de_cluster1, input$sel_de_cluster2, input$sel_test)
    
    withProgress(message = 'Finding differentially expressed genes...', value = 0, {
      df <- FindMarkers(sobj, ident.1 = input$sel_de_cluster1, ident.2 = input$sel_de_cluster2, 
                        only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
      df <- df[ df$p_val_adj < 0.01, ]
    })
    sessionData$status$marker_ready <- TRUE  #aggs:  at this point of the analysis change icon color to green
    return(df)
  })
  
  output$table_de <- renderDataTable({
    df <- de_result()
    
    return(DT::datatable(df, options = list(scrollX = TRUE), filter = "top"))
  })
  
  output$plot_de_genes <- renderPlot({
    
    req(input$table_de_rows_selected)
    
    sobj <- sessionData$sobj_cluster()
    
    w <- input$table_de_rows_selected
    
    if (length(w) > 0) {
      genes <- rownames(de_result())[ w ]
  
      VlnPlot(sobj, features = genes, pt.size=0.5)
      
    }
    
  })
  
  sessionData$all_markers <- all_markers
  
  return(sessionData)
}



