
reportEnv <- new.env()

assign("current_log_file", NULL, envir=reportEnv)
assign("indents", c(), envir=reportEnv)

set_log_file <- function(filename) {
  assign("current_log_file", filename, envir=reportEnv)
}

start_log_file <- function(filename) {
  set_log_file(filename)
  close(file(filename, open="w"))
}

append_log_file <- function(file) {
  set_log_file(file)
}

clear_log_file <- function() {
  set_log_file(NULL)
}

indent_report <- function(indent="   ") {
  assign("indents", c(get("indents", envir=reportEnv), indent), envir=reportEnv)
}

unindent_report <- function(n=1) {
  indents <- get("indents", envir=reportEnv)
  if (length(indents) - n < 1) {
    assign("indents", c(), envir=reportEnv)
  } else {
    assign("indents", indents[1:(length(indents)-n)], envir=reportEnv)
  }
}

report <- function(string, tofile=T) {
  indents <- get("indents", envir=reportEnv)
  current_log_file <- get("current_log_file", envir=reportEnv)
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


