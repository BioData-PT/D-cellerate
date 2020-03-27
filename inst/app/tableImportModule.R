


#' UI function for table import module
datasetImportUI <- function(id) {
  ns <- NS(id)
  
  get.dataset.names <- function() {
    items <- data(package = "datasets")$results[ , "Item" ]
    items <- gsub(" .*$", "", items)
    #grep("^[a-z]", items, value = TRUE)
  }
  
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



