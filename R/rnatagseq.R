
split_rtseq_inputs <- function(rtseq, max_cycles=NA) {
  barcodes <- rtseq$barcodes
  
  split_fastq_file <- function(lane, fastq) {
    start_log_file(get_path(rtseq, dir="log", 
                            file=paste0("split_", lane, ".log")))
    ptm <- proc.time()  # start timing
    
    total_reads <- 0
    n_cycles <- 0
    
    nbc <- length(barcodes)
    stats <- data.frame(reads=integer(nbc),
                        matched_tn=integer(nbc),
                        quality=integer(nbc))
    rownames(stats) <- barcodes

    fileends <- paste0(lane, "_", names(barcodes), ".fasta")
    files <- get_path(rtseq, dir="split", file=fileends)
    names(files) <- barcodes
    sapply(paste("touch", files), system)
    
    streamer <- ShortRead::FastqStreamer(fastq, n=10000)
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
      
      barcode_strings <- fq %>%
        ShortRead::narrow(start=1, end=6) %>%
        ShortRead::sread()
      assignments <- match_barcodes(barcode_strings, barcodes)
      for (code in barcodes) {
        reads <- fq[assignments[ ,code]] %>%
          ShortRead::narrow(start=7, end=61)  # remove barcode
        stats[code,"reads"] <- stats[code,"reads"] + length(reads)
        
        stats[code,"matched_tn"] <- stats[code,"matched_tn"] + length(reads)
        
        stats[code,"quality"] <- stats[code,"quality"] + length(reads)
        
        writeFasta(reads, files[code], mode="a")
      }
    }
    close(streamer)
    
    # ==================== reporting ============================
    
    reportf("Processing lane %s (%s).", lane, fastq)
    indent_report()
    reportf("%i total reads.", total_reads)
    if (!is.na(max_cycles))
      reportf("%i FASTQ streamer cycles processed.", max_cycles)

    demux_fails <- total_reads - sum(stats$reads)
    reportf("Demultiplexing failed for %i reads (%.2f%%).", 
            demux_fails, demux_fails / total_reads * 100)
    
    elapsed <- proc.time() - ptm
    reportf("Elapsed time for file splitting: %.2f seconds.", elapsed[1])

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
    
    rownames(stats) <- NULL
    df <- data.frame(lane=lane, barcode=names(barcodes), 
                     file=files)
    rownames(df) <- NULL
    return(cbind(df, stats))
  }
  
  dfs <- parallel::mclapply(seq_along(rtseq$filelog$input$file), 
                            function(i) split_fastq_file(rtseq$filelog$input$lane[i],
                                                         rtseq$filelog$input$file[i]))
  rtseq$filelog$split <- do.call(rbind, dfs)
  
  return(rtseq)
}
