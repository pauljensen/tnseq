
BOWTIE_CMD <- "~/bowtie/bowtie"
BOWTIE_OPTS <- "-f -n 3 -l 60 -m 1 -p 2"
INDEX_PATH <- "~/seqdata/index/"
BOWTIE_BUILD_CMD <- "~/bowtie/bowtie-build"


# ================ sample mapping ================


# load a full map file, output from bowtie
load_mapfile <- function(filename) {
  reads <- readr::read_tsv(filename,
                           col_names=c("read_name", "strand", "ref_name", 
                                       "pos", "sequence", "quality", 
                                       "duplicity", "mismatch"),
                           col_types="cc_ic_ic")
  reads$reads <- stringr::str_match(reads$read_name, "\\d+-(\\d+)")[ ,2] %>%
                    as.integer()
  
  reads$length <- nchar(reads$sequence)
  reads$pos <- as.integer(reads$pos + 1)   # bowtie positions are zero-indexed
  pos_strand <- reads$strand == "+"
  reads$pos[pos_strand] <- reads$pos[pos_strand] + reads$length[pos_strand] - 2
  
  reads$has_mismatch <- nchar(reads$mismatch) > 0
  reads$mismatch_pos <- as.integer(NA)
  has_mismatch <- reads$has_mismatch
  reads$mismatch_pos[has_mismatch] <- (reads$mismatch[has_mismatch] %>%
        stringr::str_match("(\\d+):"))[ ,2] %>%
        as.integer() + 1  # bowtie positions are zero-indexed
  
  reads$mismatch_in_TA <- with(reads, ifelse(has_mismatch,
                                             length - mismatch_pos < 2,
                                             FALSE))
  
  reads %<>%
    dplyr::group_by(read_name) %>%
    dplyr::mutate(duplicated=n() > 1)
  n_dup_orig <- nrow(reads[!reads$duplicated, ])
  
  reads <- reads[!reads$mismatch_in_TA, ]
  reads %<>% 
    dplyr::group_by(read_name) %>%
    dplyr::mutate(duplicated=n() > 1)
  n_dup <- nrow(reads[!reads$duplicated, ])
  
#  print((n_dup_orig - n_dup) / n_dup_orig * 100)
  
  return(reads[!reads$duplicated, ])
}

compress_reads <- function(reads) {
  reads %<>% 
    dplyr::select(pos, reads, strand) %>% 
    dplyr::group_by(pos) %>% 
    dplyr::summarize(reads=sum(reads),
                     strand=paste(sort(unique(strand), decreasing=T), 
                                  collapse="/"))
  reads$pos %<>% as.integer()
  return(reads)
}

compress_mapfile <- function(infile, outfile) {
  if (file.info(infile)$size == 0) {
    # empty file; create empty compressed file
    writeLines('"pos","reads","strand"\n', outfile)
  } else {
    infile %>%
      load_mapfile() %>%
      compress_reads() %>%
      readr::write_csv(outfile)
  }
}

load_compressed_mapfile <- function(filename) {
  reads <- readr::read_csv(filename, col_names=T, col_types="iic")
  if (nrow(reads) == 1 && is.na(reads[1,"pos"])) {
    # this was an empty file; return and empty data frame
    return(reads[-1, ])
  } else {
    return(reads)
  }
}

#' @export
parse_bowtie_log <- function(logfile) {
  text <- paste(readLines(logfile), collapse="")
  reads <- as.integer(stringr::str_match(text, "reads processed: (\\d+)")[1,2])
  aligned <- as.integer(stringr::str_match(text, "alignment: (\\d+)")[1,2])
  return(c(reads=reads, aligned=aligned))
}

# TODO: make parallel
#' @export
map_reads <- function(tnseq, show_bowtie_cmd=F) {
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
    } else if (length(split_idx) == 0) {
      stop("No samples for file ", input, " barcode ", code)
    }
    args$splitfile[i] <- tnseq$filelog$split$file[split_idx]
    input <- tnseq$filelog$split$input[split_idx]
    barcode <- tnseq$filelog$split$barcode[split_idx]
    prefix <- paste0(input, "_", barcode)
    args$mapfile[i] <- get_path(tnseq, dir="map", file=paste0(prefix, ".map"))
    args$compressed_mapfile[i] <- paste0(args$mapfile[i], "c")
    args$logfile[i] <- get_path(tnseq, dir="log", 
                                file=paste0(prefix, "_bowtie.log"))
    #args$indexfile[i] <- paste0(INDEX_PATH, tnseq$samples$genome[i])
    args$indexfile[i] <- get_path(tnseq, dir="genomes", 
                                  file=tnseq$samples$genome[i])
    args$bowtie_cmd[i] <- paste(BOWTIE_CMD, BOWTIE_OPTS, args$indexfile[i], 
                                args$splitfile[i], args$mapfile[i], 
                                "2>", args$logfile[i])
  }
  
  for (i in 1:n_samples) {
    reportf("Mapping file %s.", args$splitfile[i])
    indent_report()
    reportf("Index file: %s.", args$indexfile[i])
    if (show_bowtie_cmd) reportf("Bowtie call: %s.", args$bowtie_cmd[i])
    
    system(args$bowtie_cmd[i])
    results <- parse_bowtie_log(args$logfile[i])
    tnseq$samples$reads[i] <- results["reads"]
    tnseq$samples$aligned[i] <- results["aligned"]
    reportf("%i/%i (%5.2f%%) of reads aligned.", 
            results["aligned"],
            results["reads"],
            results["aligned"]/results["reads"] * 100)
    
    reportf("Compressing mapfile to %s.", args$compressed_mapfile[i])
    compress_mapfile(args$mapfile[i], args$compressed_mapfile[i])
    tnseq$samples$mapfile[i] <- args$compressed_mapfile[i]
    unindent_report()
  }
  
  tnseq$args <- args
  return(tnseq)
}


# ================ compiling mapped reads ================

#' @export
pair_samples <- function(tnseq) {
  t1 <- tnseq$samples[tnseq$samples$time == 1,]
  t2 <- tnseq$samples[tnseq$samples$time == 2,]
  n2 <- nrow(t2)
  
  pairs <- t2
  pairs$input <- NULL
  pairs$barcode <- NULL
  pairs$time <- NULL
  pairs$t2_mapfile <- pairs$mapfile
  pairs$t1_mapfile <- character(n2)
  pairs$mapfile <- NULL
  pairs$t2_reads <- pairs$reads
  pairs$t1_reads <- integer(n2)
  pairs$reads <- NULL
  pairs$t2_aligned <- pairs$aligned
  pairs$t1_aligned <- integer(n2)
  pairs$aligned <- NULL
  
  t1_sl <- paste(t1$strain, t1$library)
  t2_sl <- paste(t2$strain, t2$library)
  t1_slc <- paste(t1_sl, t1$condition)
  t2_slc <- paste(t2_sl, t2$condition)
  
  for (i in 1:n2) {
    idx <- which(t2_slc[i] == t1_slc)
    if (length(idx) == 0) {
      idx <- which(t2_sl[i] == t1_sl)
    }
    if (length(idx) == 0) {
      stop("No matching time 1 for sample: ", t2_slc[i])
    }
    if (length(idx) > 1) {
      stop("Multiple matching time 1 for sample: ", t2_slc[i])
    }
    pairs$t1_mapfile[i] <- t1$mapfile[idx]
    pairs$t1_reads[i] <- t1$reads[idx]
    pairs$t1_aligned[i] <- t1$aligned[idx]
  }
  
  tnseq$paired_samples <- pairs
  return(tnseq)
}

get_read_pairs <- function(mapfile1, mapfile2) {
  reads1 <- load_compressed_mapfile(mapfile1)
  reads2 <- load_compressed_mapfile(mapfile2)
  
  both <- dplyr::full_join(reads1, reads2, by=c("pos", "strand"))
  both$reads1 <- both$reads.x
  both$reads.x <- NULL
  both$reads2 <- both$reads.y
  both$reads.y <- NULL
  
  both$reads1[is.na(both$reads1)] <- 0
  both$reads2[is.na(both$reads2)] <- 0
  
  return(both)
}

#' @export
load_insertions <- function(tnseq) {
  n_paired <- nrow(tnseq$paired_samples)
  cols_to_keep <- c("strain", "library", "condition", "genome", "expansion")
  load_aux <- function(i) {
    cat(i, "\n")
    inserts <- get_read_pairs(tnseq$paired_samples$t1_mapfile[i],
                              tnseq$paired_samples$t2_mapfile[i])
    df <- do.call(data.frame, tnseq$paired_samples[i,cols_to_keep])
    cbind(df, inserts)
  }
  tnseq$insertions <- do.call(rbind, 
                              parallel::mclapply(1:nrow(tnseq$paired_samples), 
                                                 load_aux))
  return(tnseq)
}

#' @export
map_insertions_to_genes <- function(tnseq) {
  aux_map <- function(df) {
    genome <- df$genome[1]
    hits <- IRanges::findOverlaps(IRanges::IRanges(start=df$pos, width=1),
                                  GenomicRanges::ranges(tnseq$features[[genome]]))
    df$gene <- character(nrow(df))
    df$name <- character(nrow(df))
    df$product <- character(nrow(df))
    queries <- IRanges::queryHits(hits)
    subjects <- IRanges::subjectHits(hits)
    df$gene[queries] <- names(tnseq$features[[genome]])[subjects]
    df$name[queries] <- tnseq$features[[genome]]$name[subjects]
    df$product[queries] <- tnseq$features[[genome]]$product[subjects]
    return(df)
  }
  tnseq$insertions %<>%
    dplyr::group_by_("genome") %>%
    dplyr::do(aux_map(.))
  tnseq$insertions$gene[nchar(tnseq$insertions$gene) == 0] <- NA
  
  aux_relpos <- function(df) {
    genome <- df$genome[1]
    has_gene <- !is.na(df$gene)
    relpos <- numeric(nrow(df))
    gr <- tnseq$features[[genome]][df$gene[has_gene]]
    relpos[has_gene] <- ifelse(as.character(GenomicRanges::strand(gr)) == "+",
                               (df$pos[has_gene] - GenomicRanges::start(gr)) / GenomicRanges::width(gr),
                               (GenomicRanges::end(gr) - df$pos[has_gene]) / GenomicRanges::width(gr))
    relpos[!has_gene] <- NA
    df$relpos <- relpos
    return(df)
  }
  tnseq$insertions %<>%
    dplyr::group_by_("genome") %>%
    dplyr::do(aux_relpos(.))
  
  return(tnseq)
}


