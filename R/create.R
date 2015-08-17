
#' @export
print.tnseq <- function(tnseq) {
  cat("A Tn-seq experiment.")
}

#' Return paths for experiment files and directories.
#'
#' @param tnseq A \code{tnseq} object.
#' @param file Optional filename to append to end of path.
#' @param dir Subdirectory (\code{input}, \code{log}, etc.).  If not given, 
#'    returns the path to the \code{tnseq} directory.
#'    
#' @return String representing the directory or file path.  Directory paths 
#'    will always include a trailing "/".
#' 
#' @examples
#' get_path(tnseq, dir="split")
#' get_path(tnseq, file="lane1.fastq", dir="input")
#' 
#' @seealso \code{\link{change_path}}
#' 
#' @export
get_path <- function(tnseq, ...) UseMethod("get_path")
get_path.tnseq <- function(tnseq, file=NA, dir=NA) {
  path <- tnseq$path$path
  if (!is.na(dir)) {
    path <- paste0(path, tnseq$path[[dir]])
  }
  if (length(file) > 1 || !is.na(file)) {
    path <- paste0(path, file)
  }
  return(path)
}

#' Change the path for a \code{tnseq} directory.
#' 
#' @param tnseq A \code{tnseq} object.
#' @param newpath String with location of new directory.
#' @return The modified \code{tnseq} object.
#' 
#' A \code{tnseq} object stores the path given in 
#' \code{\link{create_tnseq_experiment}} to locate experiment files.  When 
#' renaming a tnseq directory, use this command to update the \code{tnseq} 
#' object.
#' 
#' @seealso \code{\link{get_path}}
#' 
#' @export
change_path <- function(tnseq, ...) UseMethod("change_path")
change_path.tnseq <- function(tnseq, newpath) {
  tnseq$path$path <- check_path_ending(newpath)
  return(tnseq)
}

# ============ internal commands for setting up tnseqr dirs ===================

# ensure path ends with an "/"
check_path_ending <- function(path) {
  if (!stringr::str_detect(path, "/$")) {
    path <- paste0(path, "/")
  }
  return(path)
}

# create a directory inside the path, unless it already exists
# returns a string with the directory name (ending with "/")
create_dir <- function(dirname, path=NA) {
  dirname <- check_path_ending(dirname)
  if (is.na(path)) {
    path <- ""
  } else {
    path <- check_path_ending(path)
  }
  fullpath <- paste0(path, dirname)
  if (!file.exists(fullpath)) {
    system(paste("mkdir", fullpath))
    report0("Creating directory ", fullpath)
  } else {
    report0("Directory ", fullpath, " already exists.")
  }
  return(check_path_ending(dirname))
}

# ============================= initializing ==================================

#' Create empty Tn-seq directory with sample files.
#' 
#' @param path Directory to be created.
#' @seealso \code{\link{create_tnseq_experiment}}
#' 
#' @export
create_tnseq_directory <- function(path) {
  path <- check_path_ending(path)
  
  # create the home directory
  create_dir(path)
  
  # create subdirectories
  create_dir("input", path)
  create_dir("split", path)
  create_dir("map", path)
  create_dir("log", path)
  
  # add sample barcode file
  
  # add sample samples file
}

#' Load a \code{tnseq} object from existing directory.
#' 
#' @param path Experimental directory to load.
#' @return A \code{tnseq} object.
#' @seealso \code{\link{create_tnseq_directory}}
#' 
#' @export
load_tnseq_experiment <- function(path) {
  path <- check_path_ending(path)
  if (!file.exists(path)) {
    stop(paste0("Directory does not exist.  Use `create_tnseq_directory` ",
                "to create a sample directory."))
  } else {
    report0("Creating tnseq experiment from ", path, ".")
  }
  
  indent_report()
  tnseq <- list(
    path = list(
      path=path,
      input = create_dir("input", path),
      split = create_dir("split", path),
      map = create_dir("map", path),
      log = create_dir("log", path),
      genomes = create_dir("genomes", path)
    )
  )
  unindent_report()
  
  class(tnseq) <- "tnseq"
  
  # check that barcode file exists in main directory
  barcode_file <- get_path(tnseq, file="barcodes.txt")
  if (!file.exists(barcode_file)) {
    stop("Cannot find barcode file: ", barcode_file)
  } else {
    report0("Found barcode file ", barcode_file, ".")
  }
  tnseq$barcode_file <- barcode_file
  barcodes <- read.table(tnseq$barcode_file, 
                         col.names=c("name","code"), 
                         as.is=TRUE)
  tnseq$barcodes <- barcodes$code
  names(tnseq$barcodes) <- barcodes$name
  reportf("   Loaded %i barcodes.", length(tnseq$barcodes))
  
  # check that sample sheet exists in the main directory
  sample_file <- get_path(tnseq, file="samples.csv")
  if (!file.exists(sample_file)) {
    stop("Cannot find sample file: ", sample_file)
  } else {
    report0("Found samples file ", sample_file, ".")
  }
  tnseq$sample_file <- sample_file
  tnseq$samples <- read.csv(sample_file, stringsAsFactors=F)
  # allow case insensitive names in samples file
  names(tnseq$samples) <- tolower(names(tnseq$samples))
  reportf("   Loaded %i samples.", nrow(tnseq$samples))
  
  # read all input files and create default filelog
  # TODO: remove dependence on ls system call
  input_files <- system(paste("ls -1", get_path(tnseq, dir="input")), intern=T)
  filelog <- data.frame(input=input_files, 
                        file=get_path(tnseq, dir="input", file=input_files),
                        stringsAsFactors=F)
  tnseq$filelog <- list(input=filelog)
  reportf("Found %i files of input data:", nrow(tnseq$filelog$input))
  indent_report()
  lapply(tnseq$filelog$input$file, report)
  unindent_report()
  
  # genome loading and processing
  tnseq$genomes <- unique(tnseq$samples$genome)
  reportf("Sample sheet references %i genome(s).", length(tnseq$genomes))
  
  return(tnseq)
}
