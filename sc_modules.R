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
    plotOutput(ns("plot_tsne"), width = "100%")
  )
}

sc_tsnevizServer <- function(input, output, session, sobj, cluster.name=NA) {
  
  output$plot_tsne <- renderPlot({
    .sobj <- sobj()
    
    if (!is.na(cluster.name)) {
      .sobj <- SetAllIdent(.sobj, cluster.name)
    }
    
    TSNEPlot(.sobj, do.label = TRUE)
  })
  #, height = function() {
  #  session$clientData[[ paste0("output_", session$ns("plot_tsne"), "_width") ]] 
  #})
  
}



