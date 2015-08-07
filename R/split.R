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
split_inputs <- function(tnseq, max_cycles=NA, dump_fails=FALSE,
                         match_transposon_sequence=FALSE, quality=NA) {
  barcodes <- tnseq$barcodes
  
  split_fastq_file <- function(lane, fastq) {
    start_log_file(get_path(tnseq, dir="log", 
                            file=paste0("split_", lane, ".log")))
    ptm <- proc.time()  # start timing
    
    if (dump_fails) {
      build_hashes(c(barcodes, "FAILS"))
    } else {
      build_hashes(barcodes)
    }
    
    total_reads <- 0
    n_cycles <- 0
    
    nbc <- length(barcodes)
    stats <- data.frame(reads=integer(nbc),
                        matched_tn=integer(nbc),
                        quality=integer(nbc))
    rownames(stats) <- barcodes
    
    streamer <- ShortRead::FastqStreamer(fastq)
    repeat {
      n_cycles <- n_cycles + 1
      if (!is.na(max_cycles) && n_cycles > max_cycles) {
        break
      } 
      
      fq <- ShortRead::yield(streamer)
      if (length(fq) == 0) {
        break
      }
      total_reads <- total_reads + length(fq)
      
      trimmed <- fq %>% ShortRead::narrow(start=10, end=34)
      barcode_strings <- trimmed %>%
        ShortRead::narrow(start=1, end=6) %>%
        ShortRead::sread()
      assignments <- match_barcodes(barcode_strings, barcodes)
      for (code in barcodes) {
        reads <- trimmed[assignments[ ,code]] %>%
          ShortRead::narrow(start=7, end=25)  # remove barcode
        stats[code,"reads"] <- stats[code,"reads"] + length(reads)
        
        if (match_transposon_sequence) {
          # matching here
          tn_sequence <- Biostrings::DNAString("TAACAG")
          reads %<>% ShortRead::trimLRPatterns(Rpattern=tn_sequence, 
                                               subject=., 
                                               max.Rmismatch=c(-1,-1,1,1,1,1))
          reads <- reads[ShortRead::width(reads) < 19]
        } else {
          reads %<>% ShortRead::narrow(start=1, end=12)
        }
        stats[code,"matched_tn"] <- stats[code,"matched_tn"] + length(reads)
        
        if (!is.na(quality)) {
          # filter reads for quality
          mean_quality <- rowMeans(as(quality(reads), "matrix"), na.rm=T)
          reads <- reads[mean_quality > quality]
        }
        stats[code,"quality"] <- stats[code,"quality"] + length(reads)
        
        keys <- reads %>% ShortRead::sread() %>% as.character()
        add_keys(code, keys)
      }
      close(streamer)

      if (dump_fails) {
        fails <- colSums(assignments) == 0
        add_keys("FAILS", fq %>% ShortRead::sread() %>% as.character())
      }
      
      tnseq$reads <<- reads
    }
    
    fileends <- paste0(lane, "_", names(barcodes), ".fasta")
    files <- get_path(tnseq, dir="split", file=fileends)
    for (i in seq_along(barcodes)) {
      write_hash_to_file(barcodes[i], path.expand(files[i]))
    }
    if (dump_fails) {
      failfile <- get_path(tnseq, dir="split", 
                           file=paste0(lane, "_", "FAILS.fasta"))
      write_hash_to_file("FAILS", path.expand(failfile))
    }
    
    # ==================== reporting ============================
    
    reportf("Processing lane %s (%s).", lane, fastq)
    indent_report()
    reportf("%i total reads.", total_reads)
    if (!is.na(max_cycles))
      reportf("%i FASTQ streamer cycles processed.", max_cycles)
    if (match_transposon_sequence)
      reportf("Reads were filtered by matching against Tn sequence.")
    if (!is.na(quality))
      reportf("Reads were quality filtered (Q >= %i).", quality)
    
    demux_fails <- total_reads - sum(stats$reads)
    reportf("Demultiplexing failed for %i reads (%.2f%%).", 
            demux_fails, demux_fails / total_reads * 100)
    
    elapsed <- proc.time() - ptm
    reportf("Elapsed time for file splitting: %.2f seconds.", elapsed[1])
    
    if (dump_fails)
      reportf("Demultiplexing failures written to file %s.", failfile)
    
    report("Distribution of valid reads:")
    report("    Reads  %Reads        Name   Barcode   Tn-match   Quality      Final  %Final")
    for (i in seq_along(barcodes)) {
      reportf("%9i  %5.2f%%  %10s  [%s]    %6.2f%%   %6.2f%%  %9i  %5.2f%%", 
              stats[i,"reads"], 
              stats[i,"reads"]/sum(stats$reads)*100, 
              names(barcodes)[i], barcodes[i],
              stats[i,"matched_tn"]/stats[i,"reads"]*100,
              stats[i,"quality"]/stats[i,"matched_tn"]*100,
              stats[i,"quality"],
              stats[i,"quality"]/sum(stats$quality)*100)
    }
    report("")
    
    set_log_file(NULL)
    
    erase_hashes(barcodes)
    if (dump_fails) erase_hashes("FAILS")
    
    rownames(stats) <- NULL
    df <- data.frame(lane=lane, barcode=names(barcodes), 
                     file=files)
    rownames(df) <- NULL
    return(cbind(df, stats))
  }
  
  dfs <- parallel::mclapply(seq_along(tnseq$filelog$input$file), 
                function(i) split_fastq_file(tnseq$filelog$input$lane[i],
                                             tnseq$filelog$input$file[i]))
  tnseq$filelog$split <- do.call(rbind, dfs)
  
  return(tnseq)
}

