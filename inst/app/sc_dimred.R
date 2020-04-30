
source("sc_modules.R")

sc.vargenes <- function(sobj, vargenes.params) {
  # sobj <- FindVariableGenes(sobj, mean.function = ExpMean, dispersion.function = LogVMR,  
  #                           x.low.cutoff = vargenes.params$xmin, 
  #                           x.high.cutoff = vargenes.params$xmax, 
  #                           y.cutoff = vargenes.params$ycutoff,
  #                           do.plot = FALSE)
  
  sobj <- FindVariableFeatures(sobj, 
                               selection.method = vargenes.params$method, 
                               nfeatures = vargenes.params$num.genes)
}

sc.vargenes.plot <- function(sobj, vargenes.params) {
  # FindVariableGenes(sobj, mean.function = ExpMean, dispersion.function = LogVMR,  
  #                   x.low.cutoff = vargenes.params$xmin, 
  #                   x.high.cutoff = vargenes.params$xmax, 
  #                   y.cutoff = vargenes.params$ycutoff)
  # abline(v=c(vargenes.params$xmin, vargenes.params$xmax), h=vargenes.params$ycutoff, col="darkred", lty="dashed")
  top10 <- head(VariableFeatures(sobj), 10)
  
  plot1 <- VariableFeaturePlot(sobj)
  LabelPoints(plot = plot1, points = top10, repel = TRUE)
}

sc.scale <- function(sobj) {
  sobj <- ScaleData(sobj, features = rownames(sobj))
}

sc.pca <- function(sobj, pca.params) {
  if (pca.params$genes.use == "All genes") {
    pc.genes <- rownames(sobj)
  } else {
    pc.genes <- VariableFeatures(sobj)
  }

  # regress vars. TODO: make this ugly piece of code a little prettier
  # rv <- character()
  # if (input$check_reg_umi) {
  #   rv <- append(rv, "nUMI")
  # }
  # if (input$check_reg_mito) {
  #   rv <- append(rv, "percent.mito")
  # }
  # 
  # if (length(rv) == 0) {
  #   rv <- NULL
  #   detail.text <- NULL
  # } else {
  #   detail.text <- paste0("Regressing out ", paste(rv, collapse=", "), ".")
  # }

  sobj <- RunPCA(object = sobj, features = pc.genes, npcs = pca.params$num.pcs)

  sobj@meta.data$original.clusters <- rep("None", dim(sobj)[2])
  sobj@meta.data$new.clusters <- rep("None", dim(sobj)[2])
  
  return(sobj)

}

#-------------------------------------------------------------------------------------------
#
## AGGS code:

sc.cluster <- function( sobj, cluster.params ) {
  if ( cluster.params$sel_cluster_method == "0" ) {
    sobj <- sobj 
  } else {
    resolution <- cluster.params$resolution
    ndims <- cluster.params$ndims
    algorithm <- cluster.params$algorithm
    sobj <- FindNeighbors( sobj, reduction = "pca", dims = 1:ndims )
    sobj <- FindClusters( sobj, resolution = resolution, algorithm = algorithm ) 
    sobj@meta.data$original.clusters <- Idents( sobj )
  }
  return( sobj )
}

#sc.cluster.names <- function( sobj ) {
#  from <- sort( unique( sobj@meta.data$original.clusters ) )
#  to <- sapply( paste0( "text_cluster_to_", from ), function(x) input[[ x ]] )
#  do.call( req, lapply( paste0("text_cluster_to_", from ), function(x) input[[ x ]] ) )
#  df <- data.frame( from = from, to = to, stringsAsFactors = FALSE )
#  return( df )
#}

sc.cluster.renamed <- function( sobj, cluster.params, cluster.renamed.params ) {
  if ( cluster.params$renamed == TRUE) {
    df <- data.frame( from = cluster.renamed.params$from, 
                      to = cluster.renamed.params$to, stringsAsFactors = FALSE )
    clusters <- df
    sobj@meta.data$new.clusters <- plyr::mapvalues( sobj@meta.data$original.clusters, 
                                                    from = clusters$from, to = clusters$to )
  } else { 
    sobj@meta.data$new.clusters <- sobj@meta.data$original.clusters
  }
  Idents( sobj ) <- sobj@meta.data$new.clusters
  return( sobj )
}

sc.tsne <- function( sobj, tsne.params ) {
  sobj <- RunTSNE(sobj, 
                  dims = 1:tsne.params$pcnum, 
                  perplexity = tsne.params$perplexity, 
                  seed.use = tsne.params$seed)
  return(sobj)
}

sc.tsne.cluster <- function( sobj, sobj.cluster ) {
  sobj@meta.data$original.clusters <- sobj.cluster@meta.data$original.clusters
  sobj@meta.data$new.clusters <- sobj.cluster@meta.data$new.clusters
  Idents(sobj) <- sobj@meta.data$new.clusters
  return(sobj)
}

#-------------------------------------------------------------------------------------------

#' UI function for statistics module
sc_dimredUI <- function(id) {
  ns <- NS(id)
  
  vargenes.ui <- tagList(
    plotOutput(ns("plot_vargenes")),
    verbatimTextOutput(ns("text_hvginfo"))
    )
  
  vg.box <- box(
    collapsible = TRUE,
    width = NULL,
    title = "Variable Genes Options",
    helpText("Identify features that are outliers on a 'mean variability plot'."),
    #numericInput(ns("num_ycutoff"), "Dispersion cutoff", value = 0.5, min = 0),
    #numericInput(ns("num_xmin"), "Expression min", value = 0.05, min = 0),
    #numericInput(ns("num_xmax"), "Expression max", value = 6, min = 0)
    selectInput(ns("sel_vargene_method"), label = "Selection method", choices=c("vst", "dispersion"), selected = "vst"),
    numericInput(ns("num_vargenes"), label="Number of genes", value=2000, min=1)
  )
  
  pca.box <- box(
    collapsible=TRUE,
    width = NULL,
    title = "PCA Options",
    helpText("Select the set of genes to use and the number of principal components to compute."),
    selectInput(ns("sel_pca_genes"), label = "Genes to use", 
                choices = c("All genes", "Variable genes"), selected = "Variable genes"),
    checkboxInput(ns("check_reg_umi"), label = "Regress nUMI", value = FALSE),
    checkboxInput(ns("check_reg_mito"), label = "Regress percent.mito", value = FALSE),
    numericInput(ns("num_pcs"), "Number of PCs to compute", value = 40, min = 0)
  )
  
  tsne.box <- box(
    collapsible=TRUE,
    width = NULL,
    title = "t-SNE Options",
    numericInput(ns("num_tsne_pca"), "Number of PCs to use", value = 20, min = 0),
    numericInput(ns("num_perplexity"), "Perplexity", value = 30, min = 1),
    numericInput(ns("num_tse_seed"), "RNG Seed", value = 42, min = 1)
  )
  
  clust.box <- box(
    collapsible=TRUE,
    width = NULL,
    title = "Clustering Options",
    helpText(""),
    selectInput(ns("sel_cluster_method"), label = "Algorithm", 
                choices = c("None"=0, 
                            "Original Louvain algorithm"=1,
                            "Louvain algorithm with multilevel refinement"=2,
                            "SLM algorithm"=3),
                selected = 1),
    conditionalPanel(condition = paste0("input['", ns("sel_cluster_method"), "']", " != '0'"),
                     numericInput(ns("num_resolution"), "Resolution", value = 0.8, min = 0),
                     numericInput(ns("num_pca"), "Number of PCs to use", value = 20, min = 0)
    ),                     
    checkboxInput(ns("check_rename"), label = "Rename/ Merge Clusters", value=FALSE),
    conditionalPanel(condition = paste0("input['", ns("check_rename"), "']", " == true"), uiOutput(ns("ui_rename")))
  )

  results.box <- tabBox(
    width = NULL,
    title = "Visualize",
    tabPanel("Variable Genes", vargenes.ui),
    tabPanel("PCA", sc_pcavizUI(ns("sc_pcaviz"))),
    tabPanel("t-SNE", sc_tsnevizUI(ns("sc_tsneviz"))))
  
  fluidRow(
    column(width = 4,
           vg.box,
           pca.box,
           clust.box,
           tsne.box
    ),
    column(width=8,
           results.box
    )
  )
  
}

#' Server function for statistics module
#' 
#' @return A dataframe as a reactive value.
sc_dimredServer <- function(input, output, session, sessionData) {
  
  status <- sessionData$status
  
  #
  # Vargenes
  #
  
  vargenes.params <- reactive({
    list(
      method = input$sel_vargene_method,
      num.genes = input$num_vargenes
      # xmin = input$num_xmin,
      # xmax = input$num_xmax,
      # ycutoff = input$num_ycutoff
    )
  })
  
  sobj_vargenes <- reactive({
    req(sessionData$sobj_norm())
    
    print("Finding variable genes...")
    
    withProgress(message = 'Finding variable genes...', {
      sobj <- sc.vargenes(sessionData$sobj_norm(), vargenes.params())
    })
    
    status$vargenes_ready <- TRUE
    
    return(sobj)
  })
  
  output$plot_vargenes <- renderPlot({
    sc.vargenes.plot(sobj_vargenes(), vargenes.params())
  })
  
  # output$text_hvginfo <- renderText({
  #   sobj <- sobj_vargenes()
  #   
  #   paste("Number of variable genes: ", length(sobj@var.genes), "\n")
  # })
  
  #
  # PCA
  #
  
  observe({
    pars <- sessionData$filter.params()
    
    if (pars$mito.pattern != "") {
      enable(id = "check_reg_mito")
    } else {
      updateCheckboxInput(session, "check_reg_mito", value=FALSE)
      disable(id = "check_reg_mito")
    }
  })
  
  pca.params <- reactive({
    list(
      genes.use = input$sel_pca_genes,
      num.pcs = input$num_pcs,
      regress.umi = input$check_reg_umi,
      regress.mito = input$check_reg_mito
    )
  })

  sobj_scaled <- reactive({
    sobj <- sobj_vargenes()
    
    withProgress(message = 'Scaling data...', {
      sobj <- sc.scale(sobj)
    })
    
    return(sobj)
  })
  
  sobj_pca <- reactive({
    sobj <- sobj_scaled()
    
    print("Calculating PCA...")

    withProgress(message = 'Calculating PCA...', {
      sobj <- sc.pca(sobj, pca.params())
    })
    
    status$pca_ready <- TRUE
    
    return(sobj)
  })
  
  #
  # Clustering
  #
  
  #---------------------------------------------------------------------------------------
  #
  ## AGGS code:
  
  cluster.params <- reactive( {
    list(
      sel_cluster_method = input$sel_cluster_method,
      resolution = input$num_resolution,
      ndims = input$num_pca,
      algorithm = as.numeric(input$sel_cluster_method), 
      renamed = input$check_rename
      )
  } )
  
  sobj_cluster <- reactive({
    req(sobj_pca())
    sobj <- sobj_pca()
    
    print("Finding clusters...")
    
    withProgress(message = 'Finding clusters...', {
      sobj <- sc.cluster( sobj, cluster.params() )
    })
    
    status$clustering_ready <- TRUE
    
    return(sobj)
  })
  
  #cluster_names <- reactive({
  #  req(sobj_cluster())
  #  sobj <- sobj_cluster()
  #  df <- sc.cluster.names( sobj )
  #  return(df)
  #})
  
  cluster.renamed.params <- reactive( {
    req(sobj_cluster())
    sobj <- sobj_cluster()
    from <- sort( unique( sobj@meta.data$original.clusters ) )
    to <- sapply( paste0( "text_cluster_to_", from ), function(x) input[[ x ]] )
    list(
      from = from,
      to = to
    )
    #do.call( req, lapply( paste0("text_cluster_to_", from ), function(x) input[[ x ]] ) ) ## WARNING: this may need to be included in the future
    } )
  
  sobj_cluster_renamed <- reactive({
    req(sobj_cluster())
    sobj <- sobj_cluster()
    sobj <- sc.cluster.renamed( sobj, cluster.params(), cluster.renamed.params() )
    
    return(sobj)
  })
  
  #---------------------------------------------------------------------------------------
  
  output$ui_rename <- renderUI({
    req(sobj_cluster())
    
    print("Updating cluster rename inputs...")
    
    ns <- session$ns
    sobj <- sobj_cluster()
    
    clusters.original <- sort(unique(Idents(sobj)))
    
    boxes <- lapply(clusters.original, function(x) {
      textInput(ns(paste0("text_cluster_to_", x)), label = paste0("Cluster ", x), value = x)
    })
    boxes <- do.call(flowLayout, boxes)
    
    return(tagList(boxes))
  })
  
  #
  # t-SNE
  #

  #---------------------------------------------------------------------------------------
  #
  ## AGGS code:
  
  tsne.params <- reactive( { 
    list( 
      pcnum = input$num_tsne_pca,
      perplexity = input$num_perplexity,
      seed = input$num_tse_seed
    )
  } )
  
  sobj_tsne <- reactive({
    req(sobj_pca())
    sobj <- sobj_pca()
    print("Running t-SNE projection...")
    withProgress(message = 'Running t-SNE projection...', {
      sobj <- sc.tsne( sobj, tsne.params() )
    })
    status$tsne_ready <- TRUE
    return(sobj)
  })
  
  # add clustering information to t-SNE object
  sobj_tsne_cluster <- reactive({
    req(sobj_cluster(), sobj_tsne())
    sobj.cluster <- sobj_cluster_renamed()
    sobj <- sessionData$sobj_tsne()
    print("Adding clustering to t-SNE...")
    sobj <- sc.tsne.cluster( sobj, sobj.cluster )
    return(sobj)
  })

  #---------------------------------------------------------------------------------------  
  
  #
  # Final setup
  #
  
  pca.viz <- callModule(sc_pcavizServer, "sc_pcaviz", sobj_cluster)
  tsne.viz <- callModule(sc_tsnevizServer, "sc_tsneviz", status, sobj_tsne_cluster)
  
  sessionData$sobj_pca <- sobj_pca
  sessionData$sobj_tsne <- sobj_tsne
  sessionData$cluster.params <- cluster.params #aggs
  sessionData$cluster.fun <- reactive(sc.cluster) #aggs
  sessionData$cluster.renamed.params <- cluster.renamed.params # aggs
  sessionData$cluster.renamed.fun <- reactive(sc.cluster.renamed) #aggs 
  sessionData$sobj_cluster <- sobj_cluster_renamed 
  sessionData$sobj_tsne_cluster <- sobj_tsne_cluster
  sessionData$tsne.params <- tsne.params #aggs
  sessionData$tsne.fun <- reactive(sc.tsne) #aggs
  sessionData$tsne.cluster.fun <- reactive(sc.tsne.cluster) #aggs
  sessionData$plot.tsne.fun <- reactive(tsne.viz$plot.tsne) #aggs: add plot t-SNE to Rmd - ABORT IT - does not work as expected.
  
  sessionData$vargenes.params <- vargenes.params
  sessionData$vargenes.fun <- reactive(sc.vargenes)
  sessionData$vargenes.plot.fun <- reactive(sc.vargenes.plot)

  sessionData$scale.fun <- reactive(sc.scale)
  
  sessionData$pca.params <- pca.params
  sessionData$pca.fun <- reactive(sc.pca)
  sessionData$pca.scree.plot.fun <- reactive(pca.viz$scree.plot)
  sessionData$pca.plot.fun <- reactive(pca.viz$pca.plot)
  
  return(sessionData)
}



