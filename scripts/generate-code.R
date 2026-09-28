# Run from the project directory or through the Quarto pre-render hook.
# Extract the examples without executing them.
generate_code <- function() {
  source <- "parts/workshop.qmd"
  output <- "code/workshop.R"
  dir.create(dirname(output), showWarnings = FALSE)

  temp_source <- tempfile(fileext = ".qmd")
  temp_output <- tempfile(fileext = ".R")
  on.exit(unlink(c(temp_source, temp_output)))

  content <- readLines(source, warn = FALSE)
  content <- gsub("eval: false", "eval: true", content, fixed = TRUE)
  content <- gsub("eval=FALSE", "eval=TRUE", content, fixed = TRUE)
  writeLines(content, temp_source)
  knitr::purl(temp_source, output = temp_output, quiet = TRUE)
  if (!file.copy(temp_output, output, overwrite = TRUE)) {
    stop("Could not write ", output)
  }
  message("Generated ", output, " from ", source)
}

generate_code()
