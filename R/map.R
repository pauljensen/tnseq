
BOWTIE_CMD <- "~/bowtie/bowtie"
BOWTIE_OPTS <- "-f -m 1 -n 1 --best -y -p 2"
INDEX_PATH <- "~/seqdata/index/"


# ================ sample mapping ================

parse_bowtie_log <- function(logfile) {
  text <- paste(readLines(logfile), collapse="")
  reads <- as.integer(str_match(text, "reads processed: (\\d+)")[1,2])
  aligned <- as.integer(str_match(text, "alignment: (\\d+)")[1,2])
  return(c(reads=reads, aligned=aligned))
}

# TODO: make parallel
map_reads <- function(tnseq) {
  n_samples <- nrow(tnseq$samples)
  tnseq$samples$mapfile <- character(n_samples)
  tnseq$samples$reads <- integer(n_samples)
  tnseq$samples$aligned <- integer(n_samples)
  for (i in 1:n_samples) {
    lane <- tnseq$samples$lane[i]
    code <- tnseq$samples$barcode[i]
    split_idx <- which(tnseq$filelog$split$lane == lane & 
                         tnseq$filelog$split$barcode == code)
    if (length(split_idx) > 1) {
      stop("Multiple samples for lane ", lane, " barcode ", code)
    }
    splitfile <- tnseq$filelog$split$file[split_idx]
    mapfile <- paste0(tnseq$map_path, lane, "_", code, ".map")
    logfile <- paste0(tnseq$log_path, lane, "_", code, "_bowtie.log")
    indexfile <- paste0(INDEX_PATH, tnseq$samples$genome[i])
    system(paste(BOWTIE_CMD, BOWTIE_OPTS, indexfile, splitfile, mapfile,
                 "2>", logfile))
    
    tnseq$samples$mapfile[i] <- mapfile
    results <- parse_bowtie_log(logfile)
    tnseq$samples$reads[i] <- results["reads"]
    tnseq$samples$aligned[i] <- results["aligned"]
  }
  return(tnseq)
}

