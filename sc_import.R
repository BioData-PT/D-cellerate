
#source("tableImportModule.R")




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
  output$txtDescription <- renderText({ 
    print(as.character(datasets[[ input$selDataset ]]$description))
  })
  
  dataframe <- reactive({
    datasets[[ input$selDataset ]]$dataframe
  })
  
  name <- reactive({
    input$selDataset
  })
  
  list(dataframe=dataframe,
       name=name)
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


import.tabular <- function(import.params) {
  mat <- read.table(import.params$datapath,
                    header = import.params$filepath, 
                    sep = import.params$sep,
                    quote = import.params$quote,
                    stringsAsFactors = import.params$stringsAsFactors,
                    check.names = import.params$check.names)
  Matrix(as.matrix(mat), sparse=TRUE)
}



e1 <- expression(
    mat <- read.table(input$file$datapath,
                      header = input$heading, 
                      sep = input$sep,
                      quote = input$quote,
                      stringsAsFactors = stringsAsFactors,
                      check.names = FALSE),
    mat <- Matrix(as.matrix(mat), sparse=TRUE)
)




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
           )
    )
  })
  
  # parse into a data.frame
  dataframe <- reactive({
    if (!is.null(input$file)) {
      withProgress(message = "Loading dataset...", expr = {
        if (input$sel_type == "Tabular") {
          mat <- read.table(input$file$datapath,
                            header = input$heading, 
                            sep = input$sep,
                            quote = input$quote,
                            stringsAsFactors = stringsAsFactors,
                            check.names = FALSE)
          mat <- Matrix(as.matrix(mat), sparse=TRUE)
        } else if (input$sel_type == "Drop-seq tools") {
          mat <- read.table(gzfile(input$file$datapath), header=TRUE)
          rownames(mat) <- mat$GENE
          mat <- Matrix(as.matrix(mat[, -1]), sparse=TRUE)
        } else if (input$sel_type == "Cellranger HDF5") {
          mat <- Read10X_h5(input$file$datapath)
        }
      }) 
      
      return(mat)
    } else {
      return (NULL)
    }
  })
  
  name <- reactive({
    fileHandle()$filename
  })
  
  return(list(dataframe=dataframe,
              name=name))
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
  
  fileImportData <- callModule(mod_import_tableServer, "table_import", stringsAsFactors = FALSE)
  datasetImportData <- callModule(datasetImportServer, "dataset_import", datasets)
  
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
  
  return(sessionData)
}


