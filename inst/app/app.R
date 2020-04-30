  # Single-cell analysis pipeline 

# devtools::install_github('rstudio/DT')

library(DT)

library(shiny)
library(shinyBS)
library(shinyjs)
library(shinydashboard)
library(shinyWidgets)
      
library(Seurat)
library(gridExtra)

library(Matrix)

library(formatR)
  
# options(shiny.reactlog=TRUE)  
# options(shiny.trace=TRUE)

#source("tableFilterModule.R")
#source("tableTransformModule.R")
source("markdown_tools.R")
source("sc_import.R")
source("sc_filter.R")
source("sc_normalize.R")
source("sc_dimred.R")
# source("sc_cluster.R") # aggs     
source("sc_markers.R")
source("sc_export.R")

options(shiny.maxRequestSize=1*1024^3)

print("Reading datasets...")

example.datasets <- readRDS("data/example_datasets.rds")

ui <- function(request) {
  dashboardPage(
    skin = "red",
    dashboardHeader(title = "D-cellerate",
                    dropdownMenuOutput(outputId = "status")),
    dashboardSidebar(
      sidebarMenuOutput("tabs"),
      tags$hr(),
      helpText("Whether to automatically (re)calculate analysis results when an input is changed."),
      materialSwitch("check_auto", label = "Automatic calculation", status="success", value=TRUE)
    ),
    dashboardBody(
      tabItems(
        tabItem("ImportData", sc_importUI("sc_import", example.datasets)),
        tabItem("Filter", sc_filterUI("sc_filter")),
        tabItem("Normalize", sc_normalizeUI("sc_normalize")),
        tabItem("ReduceDimensions", sc_dimredUI("sc_dimred")),
        #tabItem("Cluster", sc_clusterUI("sc_cluster")),
        tabItem("DifferentialExpression", sc_markersUI("sc_markers")),
        tabItem("ExportAnalysis", sc_exportUI("sc_export")))
    )
  )
}

#' Main application server function 
server <- function(input, output, session) {
  # hideTab("tabs", target = "Filter")
  
  sessionData <- list(
    #example.datasets = example.datasets,
    dataframe = NULL,
    status = reactiveValues(
      data_ready = FALSE,
      filter_ready = FALSE,
      normalize_ready = FALSE,
      vargenes_ready = FALSE, 
      pca_ready = FALSE, 
      clustering_ready = FALSE,
      tsne_ready = FALSE, 
      marker_ready = FALSE #aggs 
      )
  )

  output$status <- renderMenu({
    ctab <- isolate(input$tabs)
    
    stats <- reactiveValuesToList(sessionData$status)
    
    msgs <- lapply(names(stats), function(x) {
      val <- stats[[x]]
      
      notificationItem(text = x, 
                       status = ifelse(val, "success", "info"),
                       icon = icon(ifelse(val, "check", "times")))
    })

    updateTabItems(session, "tabs", ctab)
    
    dropdownMenu(type = "notification", .list = msgs, icon = icon("cog"))
  })
  
  # import tab
  sessionData <- callModule(sc_importServer, "sc_import", example.datasets, sessionData)

  # filtering tab
  sessionData <- callModule(sc_filterServer, "sc_filter", sessionData)
  # 
  # # normalization
  sessionData <- callModule(sc_normalizeServer, "sc_normalize", sessionData)
  # 
  # # dimensionality reduction
  sessionData <- callModule(sc_dimredServer, "sc_dimred", sessionData)
  
  # 
  # # clustering
  # #sessionData <- callModule(sc_clusterServer, "sc_cluster", sessionData)
  # 
  # # markers
  sessionData <- callModule(sc_markersServer, "sc_markers", sessionData)
  # 
  # # export
  sessionData <- callModule(sc_exportServer, "sc_export", sessionData)

  # observeEvent(sessionData$dataframe(), {
  #   req(sessionData$dataframe())
  #   
  #   showTab("tabs", target = "Filter")
  # })
  
  # onBookmark(function(state) {
  #   
  #   showModal(modalDialog(
  #     title = "Saving application state... Please be patient.",
  #     easyClose = FALSE,
  #     footer = NULL
  #   ))
  #   
  #   #state$values$dataframe <- head(sessionData$dataframe())
  # })
  
  # TODO: keep selected tab open
  output$tabs <- renderMenu({
    selected <- isolate(input$tabs)
    print(selected)
    
    red <- "#FF0000"
    green <- "#00FF00"
    colors <- c(red, green)
    
    if (sessionData$status$data_ready) {
      analysis.tabs <- tagList(
        menuItem("Filter", tabName="Filter", icon = icon("filter")),
        menuItem("Normalize", tabName="Normalize", icon = icon("vials")),
        menuItem("Reduce Dimensions and Cluster", tabName="ReduceDimensions", icon = icon("project-diagram")),
        menuItem("Marker Gene Identification", tabName="DifferentialExpression", icon = icon("fingerprint")), #aggs
        menuItem("Export Analysis", tabName = "ExportAnalysis", icon = icon("cloud-download")) #aggs  
      )
    } else {
      analysis.tabs <- tagList()
    }
      
    sidebarMenu(
      id = "tabs",
      tags$style(".fa-filter {color:", colors[ sessionData$status$filter_ready + 1 ], "}"),
      tags$style(".fa-vials {color:", colors[ sessionData$status$normalize_ready + 1 ], "}"),
      tags$style(".fa-project-diagram {color:", colors[ sessionData$status$clustering_ready + 1 ], "}"),
      tags$style(".fa-fingerprint {color:", colors[ sessionData$status$marker_ready + 1 ], "}"), #aggs
      menuItem("Import Data", tabName="ImportData", icon = icon("cloud-upload"), selected = TRUE),
      analysis.tabs
    )
  })
  
  # onBookmarked(function(url) {
  #   removeModal()
  #   
  #   showBookmarkUrlModal(url)
  # })
  # 
  # 
  # # Read values from state$values when we restore
  # onRestore(function(state) {
  #   #sessionData$dataframe <- reactive(state$values$dataframe)
  #   
  #   #print(sessionData$dataframe)
  # })  
}

print("Launching application...")

shinyApp(ui, server, enableBookmarking = "server")


