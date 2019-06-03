
#' UI function for statistics module
sc_exportUI <- function(id) {
  ns <- NS(id)
  
  summary.ui <- tagList(
    h3("Summary"),
    htmlOutput(ns("text_summary"))
  )
  
  tagList(
    useShinyjs(),
    summary.ui,
    #uiOutput(ns("uiSections")),
    downloadButton(ns("report"), "Generate report")
  )
}

#' Server function for statistics module
#' 
#' @param cmatrix Counts matrix.
#' 
#' @return A dataframe as a reactive value.
sc_exportServer <- function(input, output, session, sessionData) {
  
  output$text_summary <- renderUI({
    print(sessionData)
    
    fparams <- sessionData$filter_params()
    
    req(fparams)
    
    ftext <- tagList(
      h4("Filtering parameters"),
      tags$ul(
        tags$li("Minimum UMI per cell: ", fparams$min_umi),
        tags$li("Maximum UMI per cell: ", fparams$max_umi),
        tags$li("Minimum genes per cell: ", fparams$min_genes),
        tags$li("Maximum genes per cell: ", fparams$max_genes))
    )
    
    return(ftext)
  })
  
  disable("report")

  # observe({
  #   req(sessionData$all_markers())
  # 
  #   enable("report")
  # })

  output$report <- downloadHandler(
    filename = "report.html",
    content = function(file) {
      req(sessionData$all_markers())
      
      tempReport <- file.path(tempdir(), "sc_report.Rmd")
      file.copy("sc_report.Rmd", tempReport, overwrite = TRUE)

      # Set up parameters to pass to Rmd document
      params <- list(all_markers = sessionData$all_markers())

      rmarkdown::render(tempReport,
                        output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv()))
    }
  )
  
  return(sessionData)
  
}











