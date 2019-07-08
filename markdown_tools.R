
tidy_function_body <- function(fun) {
  tidy_source(text = as.character(body(fun))[-1])
}

make_chunk_from_function_body <- function(fun) {
  paste("``` r", tidy_function_body(fun), "```", sep="\n")
}

