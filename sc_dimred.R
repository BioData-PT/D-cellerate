
source("sc_modules.R")



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
    helpText("Select the dispersion cutoff and minimum and maximum",
             "average expression to use for selection of variable genes.", 
             "Recommendation: aim for 1,000 to 2,000 genes."),
    numericInput(ns("num_ycutoff"), "Dispersion cutoff", value = 0.5, min = 0),
    numericInput(ns("num_xmin"), "Expression min", value = 0.05, min = 0),
    numericInput(ns("num_xmax"), "Expression max", value = 6, min = 0)
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
    numericInput(ns("num_pca"), "Number of PCs to use", value = 20, min = 0),
    numericInput(ns("num_perplexity"), "Perplexity", value = 30, min = 1),
    numericInput(ns("num_tse_seed"), "RNG Seed", value = 42, min = 1)
  )
  
  clust.box <- box(
    collapsible=TRUE,
    width = NULL,
    title = "Clustering Options",
    helpText(""),
    selectInput(ns("sel_cluster_method"), label = "Algorithm", choices = c("None", "Louvain (original)")),
    conditionalPanel(condition = paste0("input['", ns("sel_cluster_method"), "']", " == 'Louvain (original)'"),
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
  
  sobj_vargenes <- reactive({
    sobj <- sessionData$sobj_norm()
    xmin <- input$num_xmin
    xmax <- input$num_xmax
    ycutoff <- input$num_ycutoff
    
    print("Finding variable genes...")
    
    withProgress(message = 'Finding variable genes...', {
      sobj <- FindVariableGenes(sobj, mean.function = ExpMean, dispersion.function = LogVMR,  
                                x.low.cutoff = xmin, x.high.cutoff = xmax, y.cutoff = ycutoff)
    })
    
    status$vargenes_ready <- TRUE
    
    return(sobj)
  })
  
  output$plot_vargenes <- renderPlot({
    sobj <- sessionData$sobj_norm()
    xmin <- input$num_xmin
    xmax <- input$num_xmax
    ycutoff <- input$num_ycutoff
    
    FindVariableGenes(sobj, mean.function = ExpMean, dispersion.function = LogVMR,  
                      x.low.cutoff = xmin, x.high.cutoff = xmax, y.cutoff = ycutoff)
    abline(v=c(xmin, xmax), h=ycutoff, col="darkred", lty="dashed")
  })
  
  output$text_hvginfo <- renderText({
    sobj <- sobj_vargenes()
    
    paste("Number of variable genes: ", length(sobj@var.genes), "\n")
  })
  
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

  sobj_pca <- reactive({
    sobj <- sobj_vargenes()
    genes.use <- input$sel_pca_genes
    num.pcs <- input$num_pcs
    
    if (genes.use == "All genes") {
      pc.genes <- rownames(sobj@data)
    } else {
      pc.genes <- sobj@var.genes
    }
    
    print("Calculating PCA...")

    # regress vars. TODO: make this ugly piece of code a little prettier
    rv <- character()
    if (input$check_reg_umi) { 
      rv <- append(rv, "nUMI")
    } 
    if (input$check_reg_mito) { 
      rv <- append(rv, "percent.mito")
    } 
    
    if (length(rv) == 0) {
      rv <- NULL
      detail.text <- NULL
    } else {
      detail.text <- paste0("Regressing out ", paste(rv, collapse=", "), ".")
    }

    withProgress(message = 'Scaling data...', detail = detail.text, {
      sobj <- ScaleData(object = sobj, vars.to.regress = rv)
    })
    

    withProgress(message = 'Calculating PCA...', {
      sobj <- RunPCA(object = sobj, pc.genes = pc.genes, pcs.compute = num.pcs, do.print=FALSE)
    })
    
    status$pca_ready <- TRUE
    
    return(sobj)
  })
  
  #
  # Clustering
  #
  
  sobj_cluster <- reactive({
    if (input$sel_cluster_method == "None") {
      req(sobj_pca())
      
      sobj <- sobj_pca()
    } else {
      req(sobj_pca())
      
      sobj <- sessionData$sobj_pca()
      resolution <- input$num_resolution    
      ndims <- input$num_pca
      
      print("Finding clusters...")
      
      withProgress(message = 'Finding clusters...', {
        sobj <- FindClusters(sobj, reduction.type = "pca", dims.use = 1:ndims, 
                             resolution = resolution, print.output = 0, save.SNN = FALSE)
      })
      
      sobj@meta.data$original.clusters <- sobj@ident
      
      status$clustering_ready <- TRUE
    }
    
    return(sobj)
  })
  
  sobj_cluster_renamed <- reactive({
    req(sobj_cluster())
    
    sobj <- sobj_cluster()
    
    if (input$check_rename == TRUE) {
      clusters <- cluster_names()
      
      sobj@meta.data$new.clusters <- plyr::mapvalues(sobj@meta.data$original.clusters, from=clusters$from, to=clusters$to)
    } else {    
      sobj@meta.data$new.clusters <- sobj@meta.data$original.clusters
    }
    
    sobj <- SetAllIdent(sobj, "new.clusters")
    
    return(sobj)
  })
  
  cluster_names <- reactive({
    req(sobj_cluster())
    
    sobj <- sobj_cluster()
    
    from <- sort(unique(sobj@ident))
    to <- sapply(paste0("text_cluster_to_", from), function(x) input[[ x ]])
    
    do.call(req, lapply(paste0("text_cluster_to_", from), function(x) input[[ x ]]))
    
    df <- data.frame(from=from, to=to, stringsAsFactors = FALSE)
    
    return(df)
  })

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
  
  #
  # t-SNE
  #
  
  sobj_tsne <- reactive({
    req(sobj_pca())
    
    sobj <- sobj_pca()
    pcnum <- input$num_pca
    perplexity <- input$num_perplexity
    seed <- input$num_tse_seed
    
    print("Running t-SNE projection...")
    
    withProgress(message = 'Running t-SNE projection...', {
      sobj <- RunTSNE(sobj, dims.use = 1:pcnum, do.fast = TRUE, perplexity=perplexity, seed.use = seed)
    })
    
    status$tsne_ready <- TRUE
    
    return(sobj)
  })
  
  # add clustering information to t-SNE object
  sobj_tsne_cluster <- reactive({
    #if (status$ready == TRUE) {
      req(sobj_cluster(), sobj_tsne())
      
      sobj.cluster <- sobj_cluster_renamed()
      sobj <- sessionData$sobj_tsne()
      
      print("Adding clustering to t-SNE...")
      
      sobj@meta.data$original.clusters <- sobj.cluster@meta.data$original.clusters
      sobj@meta.data$new.clusters <- sobj.cluster@meta.data$new.clusters
      
      sobj <- SetAllIdent(sobj, "new.clusters")
    # } else {
    #   req(sobj_tsne())
    #   
    #   sobj <- sobj_tsne()
    # }
    
    return(sobj)
  })
  
  
  #
  # Final setup
  #
  
  callModule(sc_pcavizServer, "sc_pcaviz", sobj_pca)
  callModule(sc_tsnevizServer, "sc_tsneviz", status, sobj_tsne_cluster)
  
  sessionData$sobj_pca <- sobj_pca
  sessionData$sobj_tsne <- sobj_tsne
  sessionData$sobj_cluster <- sobj_cluster_renamed
  sessionData$sobj_tsne_cluster <- sobj_tsne_cluster
  
  return(sessionData)
}


