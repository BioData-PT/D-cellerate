sc_pcavizUI <- function(id, show.scree=TRUE) {
  ns <- NS(id)
  
  if (show.scree == TRUE) {
    scree.ui <- tagList(
      plotOutput(ns("plot_scree")),
      helpText('PCA "scree plot" (or "elbow plot") showing the proportion',
               'of variance explained by each proncipal component.'),
      hr())
  } else {
    scree.ui <- NULL
  }
  
  tagList(
    scree.ui,
    plotOutput(ns("plot_pca"), height="auto"),
    inputPanel(
      numericInput(ns("num_pc1"), "Dim1", value = 1),
      numericInput(ns("num_pc2"), "Dim2", value = 2),
      numericInput(ns("num_pc3"), "Dim3", value = 3)
    ),
    helpText('Scatterplots of the output of the principal component',
             'analysis. Three components may be displayed. The second seclected component',
             'will be reused in the Y-axis of the second plot.')
  )
}

sc_pcavizServer <- function(input, output, session, sobj, cluster.name=NA) {
  
  output$plot_pca <- renderPlot({
    req(sobj())
 
    .sobj <- sobj()
    
    if (!is.na(cluster.name)) {
      .sobj <- SetAllIdent(.sobj, cluster.name)
      print("New ident")
    }
    
    table(.sobj@ident)
    
    dim1 <- input$num_pc1
    dim2 <- input$num_pc2
    dim3 <- input$num_pc3
    
    p1 <- PCAPlot(object = .sobj, dim.1 = dim1, dim.2 = dim2, do.return=TRUE) + theme(legend.pos="none")
    p2 <- PCAPlot(object = .sobj, dim.1 = dim3, dim.2 = dim2, do.return=TRUE) 
    grid.arrange(p1, p2, ncol=2)
  }, height = function() {
    session$clientData[[ paste0("output_", session$ns("plot_pca"), "_width") ]] / 2
  })
  
  output$plot_scree <- renderPlot({
    req(sobj())
    
    eigs <- sobj()@dr$pca@sdev**2
    props <- eigs / sum(eigs)
    plot(props, ylab="Proportion of variance", xlab="Principal Component", pch=20, cex=1.5)
  })
  
}


sc_tsnevizUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    plotOutput(ns("plot_tsne"), width = "100%"),
    inputPanel(
      radioButtons(ns("radio_col_by"), label = "Color by", choices=c("Cluster", "Feature")),
      selectInput(ns("sel_col_by"), label = "", choices=NULL),
      checkboxInput(ns("check_legend"), label = "Show legend", value = TRUE)
    )
  )
}

sc_tsnevizServer <- function(input, output, session, status, sobj, cluster.name=NA) {
  
  observe({
    if (!status$tsne_ready)
      return()
    
    if (input$radio_col_by == "Cluster") {
      choices <- "Unclustered"
      selected <- "Unclustered"
      
      if (status$clustering_ready == TRUE) {
        choices <- c("Original Clusters", "Renamed Clusters")
        
        selected <- "Renamed Clusters"
      }
    } else {
        choices <- c("nUMI", "nGene")
        selected <- "nUMI"
        
        if (status$pca_ready == TRUE) {
          req(sobj())
          
          choices <- c(choices, colnames(sobj()@dr$pca@cell.embeddings))
        }
    }
    
    updateSelectInput(session, "sel_col_by", choices=choices, selected = selected)
  })
  
  output$plot_tsne <- renderPlot({
    req(sobj())

    choice <- input$sel_col_by
    
    if ((choice %in% c("nGene", "nUMI")) | grepl("^PC", choice)) {
      FeaturePlot(sobj(), features.plot = choice, reduction.use = "tsne", no.legend = !input$check_legend)
    } else {
      .sobj <- sobj()
      
      if (choice == "Original Clusters") {
        .sobj <- SetAllIdent(.sobj, "original.clusters")
      } else if (choice == "Renamed Clusters") {
        .sobj <- SetAllIdent(.sobj, "new.clusters")
      }
      
      TSNEPlot(.sobj, do.label = TRUE)
    }
        
    # if (!is.na(cluster.name)) {
    #   .sobj <- SetAllIdent(.sobj, cluster.name)
    #   
    #   TSNEPlot(.sobj, do.label = TRUE)
    # } else {
    #   
    # }

  })
  #, height = function() {
  #  session$clientData[[ paste0("output_", session$ns("plot_tsne"), "_width") ]] 
  #})
  
}



