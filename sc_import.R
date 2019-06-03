
source("tableImportModule.R")

#' UI function for table import module
sc_importUI <- function(id, datasets) {
  ns <- NS(id)
  
  preview.ui <- div(style=paste0("width: 100%; overflow: auto;"),
                    tableOutput(ns("table_preview")))

  
  panels.ui <- tabsetPanel(type="pills",
                           tabPanel("Summary", htmlOutput(ns("table_summary1"))),
                           tabPanel("Preview", preview.ui))
  
  import.ui <- tagList(
    radioButtons(ns("radioImport"), label = "Import from",
                 choices = list("File upload" = "file", "Built-in dataset" = "dataset"), 
                 selected = "dataset"),
    tags$hr(),
    conditionalPanel(condition = paste0("input['", ns("radioImport"), "']", " == 'file'"), tableImportUI(ns("table_import"))),
    conditionalPanel(condition = paste0("input['", ns("radioImport"), "']", " == 'dataset'"), datasetImportUI_2(ns("dataset_import"), datasets))
  )
  
  # tagList(
  #   bsCollapse(id=ns("table"), multiple=TRUE, open = c("Import"),
  #              bsCollapsePanel("Import", 
  #                              import.ui,
  #                              style = "default")))

  # sidebarLayout(
  #   sidebarPanel(import.ui),
  #   mainPanel(panels.ui))
  
  fluidRow(
    box(import.ui, width = 4),
    box(panels.ui, width = 8)
  )
}

#' Server function for table loader module
#' 
#' @return A dataframe as a reactive value.
sc_importServer <- function(input, output, session, datasets, sessionData) {
  
  fileImportData <- callModule(tableImportServer, "table_import", stringsAsFactors = FALSE)
  datasetImportData <- callModule(datasetImportServer_2, "dataset_import", datasets)
  
  dataframe <- reactive({
    req(input$radioImport)
    
    if (input$radioImport == "file") {
      req(fileImportData$dataframe())
    } else {
      req(datasetImportData$dataframe())
    }
    
    switch(input$radioImport,
           "file"=fileImportData$dataframe(),
           "dataset"=datasetImportData$dataframe())
  })
  
  #observeEvent(dataframe())
  
  # show the table
  output$table_preview <- renderTable({
    df <- dataframe()
    
    validate(need(df, "Please import a table."))
    
    return(as.matrix(df[ 1:10, 1:10 ]))
  })
  
  # show the summary
  output$table_summary1 <- renderUI({
    df <- dataframe()
    
    txt <- paste(paste("Number of columns:", ncol(df)),
                 paste("Number of rows:", nrow(df)),
                 sep="<br/>")
    
    return (HTML(txt))
  })

  sessionData$dataframe <- dataframe
  
  return(sessionData)
}


