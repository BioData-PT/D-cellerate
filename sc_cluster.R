# Clustering module

source("sc_modules.R")

# TODO: other clustering algorithms
# TODO: allow using genes instead of PCs for clustering and tSNE
# TODO: allow to specify which dimensions to use (instead of how many)

#' UI function for statistics module
sc_clusterUI <- function(id) {
  ns <- NS(id)
  
  cluster.ui <- tagList(
    plotOutput(ns("plot_clusters")),
    helpText("Number of cells per original cluster."),
    verbatimTextOutput(ns("text_clusters")),
    conditionalPanel(condition = paste0("input['", ns("check_rename"), "']", " == true"),
                     hr(),
                     plotOutput(ns("plot_clusters_merged")),
                     helpText("Number of cells per renamed cluster."))
    )
  
  tsne.panel <- tagList(
    sc_tsnevizUI(ns("sc_tsneviz")),
    helpText("t-SNE projection displaying original clusters."),
    conditionalPanel(condition = paste0("input['", ns("check_rename"), "']", " == true"),
                     hr(),
                     sc_tsnevizUI(ns("sc_tsneviz_merged")),
                     helpText("t-SNE projection displaying renamed clusters."))
    )

  pca.panel <- tagList(
    sc_pcavizUI(ns("sc_pcaviz"), show.scree=FALSE),
    #helpText("t-SNE projection displaying original clusters."),
    conditionalPanel(condition = paste0("input['", ns("check_rename"), "']", " == true"),
                     hr(),
                     sc_pcavizUI(ns("sc_pcaviz_merged"), show.scree=FALSE))
                     #helpText("t-SNE projection displaying renamed clusters."))
  )
  
  panels.ui <- tabsetPanel(type="pills",
                           tabPanel("Clustering", cluster.ui),
                           tabPanel("PCA", pca.panel),
                           tabPanel("t-SNE", tsne.panel))
  
  sidepanel.ui <- tagList(
    h4("Clustering Options"),
    helpText(""),
    selectInput(ns("sel_cluster_method"), label = "Algorithm", choices = c("Louvain (original)")),
    numericInput(ns("num_resolution"), "Resolution", value = 0.8, min = 0),
    numericInput(ns("num_pca"), "Number of PCs to use", value = 20, min = 0),
    #h4("Aggregate Clusters"),
    checkboxInput(ns("check_rename"), label = "Rename/ Merge Clusters", value=FALSE),
    #uiOutput(ns("ui_rename"))
    conditionalPanel(condition = paste0("input['", ns("check_rename"), "']", " == true"), uiOutput(ns("ui_rename")))
    #dataTableOutput(ns("dt_names"))
  )
  
  fluidRow(
    box(sidepanel.ui, width = 4),
    box(panels.ui, width = 8)
  )
  
}

#' Server function for statistics module
#' 
#' @return A dataframe as a reactive value.
sc_clusterServer <- function(input, output, session, sessionData) {
  
  # outputOptions(output, "ui_rename", suspendWhenHidden=FALSE)
  clustering_stats <- reactiveValues(ready = FALSE)
  
  sobj_cluster <- reactive({
    req(sessionData$sobj_pca())
    
    sobj <- sessionData$sobj_pca()
    resolution <- input$num_resolution    
    ndims <- input$num_pca
    
    print("Finding clusters...")
    
    withProgress(message = 'Finding clusters...', {
      sobj <- FindClusters(sobj, reduction.type = "pca", dims.use = 1:ndims, 
                           resolution = resolution, print.output = 0, save.SNN = FALSE)
    })
    
    sobj@meta.data$original.clusters <- sobj@ident
  
    clustering_stats$ready <- TRUE

    return(sobj)
  })

  sobj_cluster_renamed <- reactive({
    req(sobj_cluster())
    
    sobj <- sobj_cluster()
    
    if (input$check_rename == TRUE) {
      clusters <- cluster_names()
      
      sobj@meta.data$new.clusters <- plyr::mapvalues(sobj@meta.data$original.clusters, from=clusters$from, to=clusters$to)
    } else {    
      # clusters <- data.frame(from=sort(unique(sobj@meta.data$original.clusters)), 
      #                        to=sort(unique(sobj@meta.data$original.clusters)), 
      #                        stringsAsFactors = FALSE)
      sobj@meta.data$new.clusters <- sobj@meta.data$original.clusters
    }
    
    sobj <- SetAllIdent(sobj, "new.clusters")
    
    return(sobj)
  })
  
  sobj_tsne_cluster <- reactive({
    req(sobj_cluster(), sessionData$sobj_tsne())
    
    sobj.cluster <- sobj_cluster_renamed()
    sobj <- sessionData$sobj_tsne()
    
    print("Adding clustering to t-SNE...")
    
    sobj@meta.data$original.clusters <- sobj.cluster@meta.data$original.clusters
    sobj@meta.data$new.clusters <- sobj.cluster@meta.data$new.clusters

    sobj <- SetAllIdent(sobj, "new.clusters")
    
    # clustering_stats$ready <- TRUE
    
    return(sobj)
  })

  # Useless
  output$dt_names <- renderDT({
    req(sobj_cluster())
    
    sobj <- sobj_cluster()
    
    from <- sort(unique(sobj@ident))
    to <- from
    
    df <- data.frame(from=from, to=to, stringsAsFactors = FALSE)
    
    return(df)
  }, selection = 'none', server = F, editable = T)
  
  # cluster_names <- reactive({
  #   req(sobj_cluster())
  # 
  #   sobj <- sobj_cluster()
  # 
  #   from <- sort(unique(sobj@ident))
  #   to <- sapply(paste0("text_cluster_to_", from), function(x) input[[ x ]])
  # 
  #   do.call(req, lapply(paste0("text_cluster_to_", from), function(x) input[[ x ]]))
  # 
  #   df <- data.frame(from=from, to=to, stringsAsFactors = FALSE)
  # 
  #   print("Cluster names")
  #   print(df)
  # 
  #   return(df)
  # })
  
  cluster_names <- reactive({
    req(sobj_cluster())

    sobj <- sobj_cluster()

    from <- sort(unique(sobj@ident))
    to <- sapply(paste0("text_cluster_to_", from), function(x) input[[ x ]])

    do.call(req, lapply(paste0("text_cluster_to_", from), function(x) input[[ x ]]))

    df <- data.frame(from=from, to=to, stringsAsFactors = FALSE)

    return(df)
  })
  # 
  output$ui_rename <- renderUI({
    req(sobj_cluster())

    print("Updating cluster rename inputs...")

    ns <- session$ns
    sobj <- sobj_cluster()

    clusters.original <- sort(unique(sobj@ident))

    boxes <- lapply(clusters.original, function(x) {
      textInput(ns(paste0("text_cluster_to_", x)), label = paste0("Cluster ", x), value = x)
    })
    boxes <- do.call(flowLayout, boxes)

    return(tagList(boxes))
  })
   
  output$plot_clusters <- renderPlot({
    sobj <- sobj_cluster()

    sobj <- SetAllIdent(sobj, "original.clusters")
    
    df <- data.frame(Cell=names(sobj@ident), Cluster=sobj@ident)
    
    ggplot(df, aes(x=Cluster, fill=Cluster)) + 
      geom_bar() + 
      xlab("Cluster") + 
      ylab("Number of cells")
  })
 
  output$plot_clusters_merged <- renderPlot({
    sobj <- sobj_cluster_renamed()
    
    sobj <- SetAllIdent(sobj, "new.clusters")
    
    df <- data.frame(Cell=names(sobj@ident), Cluster=sobj@ident)
    
    ggplot(df, aes(x=Cluster, fill=Cluster)) + 
      geom_bar() + 
      xlab("Cluster") + 
      ylab("Number of cells")
  })
  
  # output$text_clusters <- renderText({
  #   sobj <- sobj_cluster()
  # 
  #   table(sobj@ident)
  # })
  
  # output$plot_pca <- renderPlot({
  #   sobj <- sobj_cluster()
  #   
  #   dim1 <- input$num_pc1
  #   dim2 <- input$num_pc2
  #   dim3 <- input$num_pc3
  #   
  #   p1 <- PCAPlot(object = sobj, dim.1 = dim1, dim.2 = dim2, do.return=TRUE) + theme(legend.pos="none")
  #   p2 <- PCAPlot(object = sobj, dim.1 = dim3, dim.2 = dim2, do.return=TRUE) + theme(legend.pos="none")
  #   grid.arrange(p1, p2, ncol=2)
  # }, height = function() {
  #   session$clientData[[ paste0("output_", session$ns("plot_pca"), "_width") ]] / 2
  # })
  
  callModule(sc_pcavizServer, "sc_pcaviz", sobj_cluster_renamed, "original.clusters")
  callModule(sc_pcavizServer, "sc_pcaviz_merged", sobj_cluster_renamed, "new.clusters")
  callModule(sc_tsnevizServer, "sc_tsneviz", sobj_tsne_cluster, "original.clusters")
  callModule(sc_tsnevizServer, "sc_tsneviz_merged", sobj_tsne_cluster, "new.clusters")
  
  sessionData$sobj_cluster <- sobj_cluster_renamed
  sessionData$sobj_tsne_cluster <- sobj_tsne_cluster
  sessionData$clustering_stats <- clustering_stats
  #sessionData$cluster_names <- cluster_names
  
  return(sessionData)
}


