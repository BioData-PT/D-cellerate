
tidy_function_body <- function(fun) {
  paste(tidy_source(text = as.character(body(fun))[-1])$text.tidy, collapse="\n")
}

make_chunk_from_function_body <- function(fun, chunk.name="") {
  header <- paste("```{r ", chunk.name, "}")
  paste(header, tidy_function_body(fun), "```", sep="\n")
}

list_to_dataframe <- function(l) {
  data.frame(Value=unlist(l), row.names=names(l))
}