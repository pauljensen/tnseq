
current_log_file <- NULL
indents <- c()

set_log_file <- function(filename) {
  current_log_file <<- filename
}

start_log_file <- function(filename) {
  current_log_file <<- filename
  close(file(filename, open="w"))
}

append_log_file <- function(file) {
  current_log_file <<- file
}

clear_log_file <- function() {
  current_log_file <<- NULL
}

indent_report <- function(indent="   ") {
  indents <<- c(indents, indent)
}

unindent_report <- function(n=1) {
  if (length(indents) - n < 1) {
    indents <<- c()
  } else {
    indents <<- indents[1:(length(indents)-n)]
  }
}

report <- function(string, tofile=T) {
  string <- paste0(paste(indents, collapse=""), string, "\n")
  cat(string)
  if (tofile && !is.null(current_log_file)) {
    cat(string, file=current_log_file, append=T)
  }
}

report0 <- function(...) {
  report(paste0(...))
}

reportf <- function(fmt, ...) {
  report(sprintf(fmt, ...))
}

debug <- function(...) {
  report(..., tofile=F)
}


