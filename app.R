# Single-cell analysis pipeline 

# devtools::install_github('rstudio/DT')

library(DT)

library(shiny)
library(shinyBS)
library(shinyjs)
library(shinydashboard)
  
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
source("sc_cluster.R")
source("sc_markers.R")
source("sc_export.R")

options(shiny.maxRequestSize=1*1024^3)

print("Reading datasets...")

example.datasets <- readRDS("data/example_datasets.rds")

ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "D-cellerate",
                  dropdownMenuOutput(outputId = "status")),
  dashboardSidebar(
    sidebarMenu(id = "tabs",
      menuItem("Import Data", tabName="ImportData", icon = icon("cloud-upload"), selected = TRUE),
      menuItem("Filter", tabName="Filter"),
      menuItem("Normalize", tabName="Normalize"),
      menuItem("Reduce Dimensions and Cluster", tabName="ReduceDimensions"),
      menuItem("Marker Gene Identification", tabName="DifferentialExpression"),
      menuItem("Export Analysis", tabName = "ExportAnalysis")
    )
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


#' Main application server function 
server <- function(input, output, session) {
  # hideTab("tabs", target = "Filter")
  
  sessionData <- list(
    dataframe = NULL,
    status = reactiveValues(
      vargenes_ready = FALSE, 
      pca_ready = FALSE, 
      clustering_ready = FALSE,
      tsne_ready = FALSE)
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
  
  # normalization
  sessionData <- callModule(sc_normalizeServer, "sc_normalize", sessionData)

  # dimensionality reduction
  sessionData <- callModule(sc_dimredServer, "sc_dimred", sessionData)

  # clustering
  #sessionData <- callModule(sc_clusterServer, "sc_cluster", sessionData)

  # markers
  sessionData <- callModule(sc_markersServer, "sc_markers", sessionData)
  
  # export
  sessionData <- callModule(sc_exportServer, "sc_export", sessionData)

  # observeEvent(sessionData$dataframe(), {
  #   req(sessionData$dataframe())
  #   
  #   showTab("tabs", target = "Filter")
  # })
}

print("Launching application...")

shinyApp(ui, server)

