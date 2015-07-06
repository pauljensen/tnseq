# ==================== preprocessing & de-multiplexing =======================

# Match a set of barcodes against a set of sequences.
#   dnas - DNAStringSet of sequences to match against.
#   barcodes - Character vector are barcode sequences.
#   mismatches - Maximum number of mismatches allowed.
#
# Barcodes are matched against beginning of the DNA strings.
# *All barcodes must be the same length.*
#
# Returns a logical matrix identifying matches.  Rows correspond to sequences,
# columns to barcodes.  ColNames are set to the barcode strings.  No rowNames.
#
# Function will warn if mismatches create ambiguity, i.e. if a sequence could
# begin with multiple barcodes.
match_barcodes <- function(dnas, barcodes, mismatches=0) {
  ncodes <- length(barcodes)
  ndnas <- length(dnas)
  candidates <- Biostrings::subseq(dnas, start=1, end=nchar(barcodes[1]))
  hits <- matrix(FALSE, nrow=ndnas, ncol=ncodes, dimnames=list(NULL, barcodes))
  for (code in barcodes) {
    hits[ ,code] <- Biostrings::isMatchingStartingAt(code, candidates, 
                                                     starting.at=1, 
                                                     max.mismatch=mismatches)
  }
  if (any(rowSums(hits) > 1)) {
    warning("Mismatch level creates ambiguous matches.")
  }
  
  return(hits)
}

#' Trim raw reads and split by DNA barcodes.
#' 
#' Adapters are trimmed, and barcodes are used to split input FASTQ files into
#' collapsed FASTA files (in the split/ directory).
#' 
#' @param tnseq A tnseq object.
#' @return A tnseq object with updated filelog.
#' @export
split_inputs <- function(tnseq) {
  barcodes <- tnseq$barcodes
  
  split_fastq_file <- function(lane, fastq) {
    set_log_file(get_path(tnseq, dir="log", 
                          file=paste0("split_", lane, ".log")))
    ptm <- proc.time()  # start timing
    build_hashes(barcodes)
    
    total_reads <- 0
    total_matched <- 0
    
    streamer <- ShortRead::FastqStreamer(fastq)
    repeat {
      fq <- ShortRead::yield(streamer)
      if (length(fq) == 0) {
        break
      }
      total_reads <- total_reads + length(fq)
      
      trimmed <- fq %>%
        ShortRead::narrow(start=10, end=34) %>%
        ShortRead::trimLRPatterns(Rpattern=Biostrings::DNAString("TAACAG"), 
                                  subject=., 
                                  max.Rmismatch=c(-1,-1,1,1,1,1))
      matched <- ShortRead::width(trimmed) <= 21
      total_matched <- total_matched + sum(as.integer(matched))
      trimmed <- trimmed[matched]
      
      barcode_strings <- trimmed %>% 
        ShortRead::narrow(start=1, end=6) %>%
        ShortRead::sread()
      assignments <- match_barcodes(barcode_strings, barcodes)
      for (code in barcodes) {
        keys <- trimmed[assignments[ ,code]] %>% 
          ShortRead::sread() %>% 
          as.character()
        add_keys(code, keys)
      }
    }
    
    reportf("Processing lane %s (%s).", lane, fastq)
    indent_report()
    reportf("%i total reads, %i (%.2f%%) with valid sequence.",
            total_reads, total_matched, total_matched/total_reads*100)
    report("Distribution of valid reads:")
    
    fileends <- paste0(lane, "_", names(barcodes), ".fasta")
    files <- get_path(tnseq, dir="split", file=fileends)
    reads <- integer(length(barcodes))
    
    for (i in seq_along(barcodes)) {
      reads[i] <- write_hash_to_file(barcodes[i], path.expand(files[i]))
      reportf("%9i %5.2f%%   %10s  [%s]", reads[i], reads[i]/total_matched*100, 
              names(barcodes)[i], barcodes[i])
    }
    erase_hashes(barcodes)
    
    unsplit <- total_matched - sum(reads)
    reportf("%9i %5.2f%%    unmatched", unsplit, unsplit/total_matched*100)
    unindent_report()
    elapsed <- proc.time() - ptm
    reportf("Elapsed time for file splitting: %.2f seconds.", elapsed[1])
    set_log_file(NULL)
    
    df <- data.frame(lane=lane, barcode=names(barcodes), 
                     file=files, reads=reads)
    rownames(df) <- NULL
    return(df)
  }
  
  dfs <- lapply(seq_along(tnseq$filelog$input$file), 
                function(i) split_fastq_file(tnseq$filelog$input$lane[i],
                                             tnseq$filelog$input$file[i]))
  tnseq$filelog$split <- do.call(rbind, dfs)
  return(tnseq)
}
