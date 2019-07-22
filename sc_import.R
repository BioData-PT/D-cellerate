
#source("tableImportModule.R")


import.tabular <- function(import.params) {
  mat <- read.table(import.params$filepath,
                    header = import.params$header, 
                    sep = import.params$sep,
                    quote = import.params$quote,
                    stringsAsFactors = import.params$stringsAsFactors,
                    check.names = import.params$check.names)
  mat <- Matrix(as.matrix(mat), sparse=TRUE)
}

import.dropseq <- function(import.params) {
  mat <- read.table(gzfile(import.params$filepath), header=TRUE, row.names=1)
  mat <- Matrix(as.matrix(mat), sparse=TRUE)
}

import.cellranger.hdf5 <- function(import.params) {
  mat <- Read10X_h5(import.params$filepath)
}

import.example <- function(import.params) {
  mat <- example.datasets[[ import.params$dataset ]]$dataframe
}

#' UI function for table import module
datasetImportUI <- function(id, datasets) {
  ns <- NS(id)
  
  tagList(
    selectInput(ns("selDataset"), label = "Dataset", choices = names(datasets)),
    tags$b("Description"),
    textOutput(ns("txtDescription"))
  )
}

#' Server function for table loader module
#' 
#' @return A dataframe as a reactive value.
datasetImportServer <- function(input, output, session, datasets) {
  import.params <- reactive({
    list(
      type = "Built-in",
      dataset = input$selDataset
    )
  })
  
  output$txtDescription <- renderText({ 
    print(as.character(datasets[[ input$selDataset ]]$description))
  })
  
  dataframe <- reactive({
    import.example(import.params())
  })
  
  name <- reactive({
    input$selDataset
  })
  
  list(dataframe=dataframe,
       name=name,
       params=import.params,
       import.fun=reactive(import.example))
}




#' UI function for table import module
mod_import_tableUI <- function(id) {
  ns <- NS(id)
  
  txt.ui <- tagList(
    checkboxInput(ns("heading"), "Has heading", value = TRUE),
    fluidRow(
      column(6, selectInput(ns("sep"), 
                            "Separator", 
                            c("Space" = " ",
                              "Tab" = "\t",
                              "Comma" = ",",
                              "Semicolon" = ";"), 
                            selected = "\t")),
      column(6, selectInput(ns("quote"), 
                            "Quote", 
                            c("None" = "",
                              "Double quote" = "\"",
                              "Single quote" = "'"), 
                            selected="None"))
    )
  )
  
  tagList(
    selectInput(ns("sel_type"), label = "File type", choices=c("Tabular", "Drop-seq tools", "Cellranger HDF5")),
    fileInput(ns("file"), "", width = "100%"),
    conditionalPanel(condition = paste0("input['", ns("sel_type"), "']", " == 'Tabular'"),
                     txt.ui)
  )
}

#' Server function for table loader module
#' 
#' @return A dataframe as a reactive value.
mod_import_tableServer <- function(input, output, session, stringsAsFactors=FALSE) {

  import.params <- reactive({
    switch(input$sel_type,
           "Tabular" = list(
             type = "Tabular",
             filepath = input$file$datapath,
             header = input$heading,
             sep = input$sep,
             quote = input$quote,
             stringsAsFactors = FALSE,
             check.names = FALSE
           ),
           "Drop-seq tools" = list(
             type = "Drop-seq tools",
             filepath = input$file$datapath
           ),
           "Cellranger HDF5" = list(
             type = "Cellranger HDF5",
             filepath = input$file$datapath
           )
    )
  })
  
  import.fun <- reactive({
    switch(input$sel_type,
           "Tabular" = import.tabular,
           "Drop-seq tools" = import.dropseq,
           "Cellranger HDF5" = import.cellranger.hdf5)
  })
  
  # parse into a data.frame
  dataframe <- reactive({
    if (!is.null(input$file)) {
      params <- import.params()
      
      withProgress(message = "Loading dataset...", expr = {
        mat <- import.fun()(params)
      }) 
      
      return(mat)
    } else {
      return (NULL)
    }
  })
  
  name <- reactive({
    input$file$name
  })
  
  return(list(dataframe=dataframe,
              name=name,
              params=import.params,
              import.fun=import.fun
              ))
}

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
    conditionalPanel(condition = paste0("input['", ns("radioImport"), "']", " == 'file'"), mod_import_tableUI(ns("table_import"))),
    conditionalPanel(condition = paste0("input['", ns("radioImport"), "']", " == 'dataset'"), datasetImportUI(ns("dataset_import"), datasets))
  )
  
  fluidRow(
    box(import.ui, width = 4),
    box(panels.ui, width = 8)
  )
}

#' Server function for table loader module
#' 
#' @return A dataframe as a reactive value.
sc_importServer <- function(input, output, session, datasets, sessionData) {
  
  fileImportData <- callModule(mod_import_tableServer, "table_import", stringsAsFactors = FALSE)
  datasetImportData <- callModule(datasetImportServer, "dataset_import", datasets)
  
  import.data <- reactive({
    if (input$radioImport == "file") {
      req(fileImportData$dataframe())
    } else {
      req(datasetImportData$dataframe())
    }
    
    switch(input$radioImport,
           "file"=fileImportData,
           "dataset"=datasetImportData)
  })

  dataframe <- reactive({
    import.data()$dataframe()
  })
  
  #observeEvent(dataframe())
  
  # show the table
  output$table_preview <- renderTable({
    df <- dataframe()
    
    validate(need(df, "Please import a table."))
    
    return(as.matrix(df[ 1:10, 1:10 ]))
  }, rownames = TRUE)
  
  # show the summary
  output$table_summary1 <- renderUI({
    df <- dataframe()
    
    txt <- paste(paste("Number of columns:", ncol(df)),
                 paste("Number of rows:", nrow(df)),
                 sep="<br/>")
    
    return (HTML(txt))
  })

  sessionData$dataframe <- dataframe
  sessionData$name <- reactive(import.data()$name())
  sessionData$import.params <- reactive(import.data()$params())
  sessionData$import.fun <- reactive(import.data()$import.fun())
  
  return(sessionData)
}


