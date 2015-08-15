
BOWTIE_CMD <- "~/bowtie/bowtie"
BOWTIE_OPTS <- "-f -m 1 -n 1 --best -y -p 2"
INDEX_PATH <- "~/seqdata/index/"


# ================ sample mapping ================


# load a full map file, output from bowtie
load_mapfile <- function(filename) {
  reads <- readr::read_tsv(filename,
                           col_names=c("read_name", "strand", "ref_name", 
                                       "pos", "sequence", "quality", 
                                       "duplicity", "mismatch"),
                           col_types="cc_i____")
  reads$reads <- stringr::str_match(reads$read_name, "\\d+-(\\d+)")[,2] %>%
                    as.integer()
  reads$pos <- reads$pos + 1   # bowtie positions are zero-indexed
  return(reads)
}

compress_reads <- function(reads) {
  reads %>% 
    dplyr::select(pos, reads, strand) %>% 
    dplyr::group_by(pos, strand) %>% 
    dplyr::summarize(reads=sum(reads))
}

compress_mapfile <- function(infile, outfile) {
  reads <- load_mapfile(infile)
  readr::write_csv(compress_reads(reads), outfile, col_names=F)
}

load_compressed_mapfile <- function(filename) {
  reads <- readr::read_csv(filename,
                           col_names=c("pos", "strand", "reads"),
                           col_types="ici")
  return(reads)
}

parse_bowtie_log <- function(logfile) {
  text <- paste(readLines(logfile), collapse="")
  reads <- as.integer(stringr::str_match(text, "reads processed: (\\d+)")[1,2])
  aligned <- as.integer(stringr::str_match(text, "alignment: (\\d+)")[1,2])
  return(c(reads=reads, aligned=aligned))
}

# TODO: make parallel
map_reads <- function(tnseq) {
  # match samples to split files
  n_samples <- nrow(tnseq$samples)
  
  tnseq$samples$mapfile <- character(n_samples)
  tnseq$samples$reads <- integer(n_samples)
  tnseq$samples$aligned <- integer(n_samples)
  
  args <- data.frame(
    splitfile=character(n_samples),
    mapfile=character(n_samples),
    compressed_mapfile=character(n_samples),
    logfile=character(n_samples),
    indexfile=character(n_samples),
    bowtie_cmd=character(n_samples),
    stringsAsFactors=F
  )
  
  for (i in 1:n_samples) {
    input <- tnseq$samples$input[i]
    code <- tnseq$samples$barcode[i]
    split_idx <- which(tnseq$filelog$split$input == input & 
                         tnseq$filelog$split$barcode == code)
    if (length(split_idx) > 1) {
      stop("Multiple samples for file ", input, " barcode ", code)
    }
    args$splitfile[i] <- tnseq$filelog$split$file[split_idx]
    input <- tnseq$filelog$split$input[split_idx]
    barcode <- tnseq$filelog$split$barcode[split_idx]
    prefix <- paste0(input, "_", barcode)
    args$mapfile[i] <- get_path(tnseq, dir="map", file=paste0(prefix, ".map"))
    args$compressed_mapfile[i] <- paste0(args$mapfile[i], "c")
    args$logfile[i] <- get_path(tnseq, dir="log", 
                                file=paste0(prefix, "_bowtie.log"))
    args$indexfile[i] <- paste0(INDEX_PATH, tnseq$samples$genome[i])
    args$bowtie_cmd[i] <- paste(BOWTIE_CMD, BOWTIE_OPTS, args$indexfile[i], 
                                args$splitfile[i], args$mapfile[i], 
                                "2>", args$logfile[i])
  }
  
  for (i in 1:n_samples) {
    system(args$bowtie_cmd[i])
    compress_mapfile(args$mapfile[i], args$compressed_mapfile[i])
    tnseq$samples$mapfile[i] <- args$compressed_mapfile[i]
    results <- parse_bowtie_log(args$logfile[i])
    tnseq$samples$reads[i] <- results["reads"]
    tnseq$samples$aligned[i] <- results["aligned"]
  }
  
  tnseq$args <- args
  return(tnseq)
}

