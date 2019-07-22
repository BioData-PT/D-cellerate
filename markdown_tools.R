
tidy_function_body <- function(fun) {
  tidy_source(text = as.character(body(fun))[-1])$text.tidy
}

make_chunk_from_function_body <- function(fun, chunk.name="", chunk.options=list(), clean.return=TRUE) {
  opts <- paste(paste(names(chunk.options), chunk.options, sep="="), collapse=", ")
  header <- paste0("```{r ", chunk.name, " ", chunk.options, "}")
  body <- tidy_function_body(fun)
  
  if (clean.return == TRUE) {
    body <- body[ grep("return(.*)", body, invert=TRUE) ]
  }

  paste(header, paste(body, collapse = "\n") , "```", sep="\n")
}
make_chunk_from_function_body(sc.filter, clean.return = TRUE)

list_to_code <- function(l, varname) {
  values <- ifelse(sapply(l, typeof) == "character", paste0("'", l, "'"), l)
  assignments <- paste(names(l), values, sep=" = ", collapse=",\n")
  res <- tidy_source(text = paste0(varname, " <- list(", assignments, ")"))$text.tidy

  res 
}



list_to_dataframe <- function(l) {
  data.frame(Value=unlist(l), row.names=names(l))
}