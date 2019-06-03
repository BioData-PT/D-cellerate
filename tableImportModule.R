


get.dataset.names <- function() {
  items <- data(package = "datasets")$results[ , "Item" ]
  items <- gsub(" .*$", "", items)
  #grep("^[a-z]", items, value = TRUE)
}

#' UI function for table import module
datasetImportUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    selectInput(ns("selDataset"), label = "Dataset", choices = get.dataset.names()),
    tags$b("Description"),
    textOutput(ns("txtDescription"))
  )
}

#' Server function for table loader module
#' 
#' @return A dataframe as a reactive value.
datasetImportServer <- function(input, output, session) {
  
  output$txtDescription <- renderText({ 
    dsets <- data.frame(data(package = "datasets")$results)
    selected <- input$selDataset
    
    print(as.character(dsets$Title[ which(dsets$Item == selected) ]))
  })
  
  dataframe <- reactive({
    df <- get(input$selDataset)
    
    as.data.frame(df)
  })

  name <- reactive({
    input$selDataset
  })

  list(dataframe=dataframe,
       name=name)
}



#' UI function for table import module
datasetImportUI_2 <- function(id, datasets) {
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
datasetImportServer_2 <- function(input, output, session, datasets) {
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
tableImportUI <- function(id) {
  ns <- NS(id)
  
  tagList(
  #  wellPanel(
    fileInput(ns("file"), "", width = "100%"),
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
}

#' Server function for table loader module
#' 
#' @return A dataframe as a reactive value.
tableImportServer <- function(input, output, session, stringsAsFactors) {
  # get the file
  fileHandle <- reactive({
    validate(need(input$file, message = FALSE))
    
    input$file
  })
  
  # parse into a data.frame
  dataframe <- reactive({
    read.table(fileHandle()$datapath,
             header = input$heading, 
             sep = input$sep,
             quote = input$quote,
             stringsAsFactors = stringsAsFactors)
  })
  
  name <- reactive({
    fileHandle()$filename
  })
  
  return(list(dataframe=dataframe,
              name=name))
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
  # get the file
  # fileHandle <- reactive({
  #   validate(need(input$file, message = FALSE))
  #   
  #   input$file
  # })
  
  # parse into a data.frame
  dataframe <- reactive({
    if (!is.null(input$file)) {
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


