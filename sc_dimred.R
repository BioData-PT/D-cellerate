
source("sc_modules.R")



#' UI function for statistics module
sc_dimredUI <- function(id) {
  ns <- NS(id)
  
  vargenes.ui <- tagList(
    plotOutput(ns("plot_vargenes")),
    verbatimTextOutput(ns("text_hvginfo"))
    )
  
  panels.ui <- tabsetPanel(type="pills",
                           tabPanel("Variable Genes", vargenes.ui),
                           tabPanel("PCA", sc_pcavizUI(ns("sc_pcaviz"))),
                           tabPanel("t-SNE", sc_tsnevizUI(ns("sc_tsneviz"))))
  
  sidepanel.ui <- tagList(
    h4("Variable Genes Options"),
    helpText("Select the dispersion cutoff and minimum and maximum",
             "average expression to use for selection of variable genes.", 
             "Recommendation: aim for 1,000 to 2,000 genes."),
    numericInput(ns("num_ycutoff"), "Dispersion cutoff", value = 0.5, min = 0),
    numericInput(ns("num_xmin"), "Expression min", value = 0.05, min = 0),
    numericInput(ns("num_xmax"), "Expression max", value = 6, min = 0),
    h4("PCA Options"),
    helpText("Select the set of genes to use and the number of principal components to compute."),
    selectInput(ns("sel_pca_genes"), label = "Genes to use", 
                choices = c("All genes", "Variable genes"), selected = "Variable genes"),
    checkboxInput(ns("check_reg_umi"), label = "Regress nUMI", value = FALSE),
    checkboxInput(ns("check_reg_mito"), label = "Regress percent.mito", value = FALSE),
    numericInput(ns("num_pcs"), "Number of PCs to compute", value = 40, min = 0),
    h4("t-SNE Options"),
    #helpText("More options??! Cool!"),
    numericInput(ns("num_pca"), "Number of PCs to use", value = 20, min = 0),
    numericInput(ns("num_perplexity"), "Perplexity", value = 30, min = 1),
    numericInput(ns("num_tse_seed"), "RNG Seed", value = 42, min = 1)
  )
  
  # sidebarLayout(
  #   sidebarPanel(
  #   ),
  #   mainPanel(panels.ui)
  # )
  
  fluidRow(
    box(sidepanel.ui, width = 4),
    box(panels.ui, width = 8)
  )
  
}

#' Server function for statistics module
#' 
#' @return A dataframe as a reactive value.
sc_dimredServer <- function(input, output, session, sessionData) {
  
  observe({
    pars <- sessionData$filter_params()
    
    if (pars$use_mito == TRUE) {
      enable(id = "check_reg_mito")
    } else {
      updateCheckboxInput(session, "check_reg_mito", value=FALSE)
      disable(id = "check_reg_mito")
    }
  })
  
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
    
    return(sobj)
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
    
    return(sobj)
  })
  
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
  
  callModule(sc_pcavizServer, "sc_pcaviz", sobj_pca)
  callModule(sc_tsnevizServer, "sc_tsneviz", sobj_tsne)
  
  sessionData$sobj_pca <- sobj_pca
  sessionData$sobj_tsne <- sobj_tsne
  
  return(sessionData)
}


