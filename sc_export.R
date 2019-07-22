
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
    downloadButton(ns("report"), "Generate report"),
    h4("Debug"),
    verbatimTextOutput(ns("txt_debug"))
  )
}

#' Server function for statistics module
#' 
#' @param cmatrix Counts matrix.
#' 
#' @return A dataframe as a reactive value.
sc_exportServer <- function(input, output, session, sessionData) {
  
  output$text_summary <- renderUI({
    fparams <- sessionData$filter.params()
    
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
  
  output$txt_debug <- renderText({
    # print(str(sessionData$import.params()))
    # 
    # print(sessionData$import.fun())
    # tidy_function_body(sessionData$import.fun())
    # cat(make_chunk_from_function_body(sessionData$import.fun()))
    paste(report.source(), collapse="\n")
    
    #str(list_to_dataframe(sessionData$import.params()))
    
  })
  
  report.source <- reactive({
    req(sessionData$import.params(), 
        sessionData$filter.params())
    
    report <- readLines("sc_report_base.Rmd")
    
    insert.function <- function(report, tag, fun, chunk.name = "", chunk.options = list()) {
      w <- which(report == tag)
      report[w] <- make_chunk_from_function_body(fun, chunk.name = chunk.name, chunk.options = chunk.options)
      
      return(report)
    }

    # Import    
    report <- insert.function(report, "<!-- import.fun -->", sessionData$import.fun(), chunk.name = "import")

    # Filter
    report <- insert.function(report, "<!-- filter.fun -->", sessionData$filter.fun(), chunk.name = "filter")
    report <- insert.function(report, "<!-- barcode.plot -->", sessionData$barcode.plot.fun(), chunk.name = "barcode-plot")
    report <- insert.function(report, "<!-- violin.plot -->", sessionData$violin.plot.fun(), chunk.name = "violin-plot")
    
    report[ which(report == "#-- import_params") ] <- list_to_code(sessionData$import.params(), "import.params")
    report[ which(report == "#-- filter_params") ] <- list_to_code(sessionData$filter.params(), "filter.params")
    
    return(report)
  })
  
  # disable("report")

  # observe({
  #   req(sessionData$all_markers())
  # 
  #   enable("report")
  # })

  output$report <- downloadHandler(
    filename = "report.html",
    content = function(file) {
      req(report.source())
      
      # Set up parameters to pass to Rmd document
      params <- list(import.params = sessionData$import.params(),
                     filter.params = sessionData$filter.params())
      
      report <- report.source()
      
      tempReport <- file.path(tempdir(), "sc_report.Rmd")
      writeLines(report, con = tempReport, sep="\n")

      rmarkdown::render(tempReport,
                        output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv()))
      
      # replace dataset server pathname with dataset filename in the rendered html
      if (sessionData$import.params()$type != "Built-in") {
        pathname <- sessionData$import.params()$filepath
        filename <- sessionData$name()
        
        final_report <- readLines(file)
        final_report <- sub(pathname, filename, final_report)
        
        writeLines(final_report, file)
      }
    }
  )
  
  # output$report <- downloadHandler(
  #   filename = "report.html",
  #   content = function(file) {
  #     req(sessionData$all_markers())
  #     
  #     tempReport <- file.path(tempdir(), "sc_report.Rmd")
  #     file.copy("sc_report.Rmd", tempReport, overwrite = TRUE)
  # 
  #     # Set up parameters to pass to Rmd document
  #     params <- list(all_markers = sessionData$all_markers())
  # 
  #     rmarkdown::render(tempReport,
  #                       output_file = file,
  #                       params = params,
  #                       envir = new.env(parent = globalenv()))
  #   }
  # )
  
  
  return(sessionData)
  
}











