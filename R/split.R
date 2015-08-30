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

transposons <- list(
  magellan=list(
    seq=Biostrings::DNAString("TAACAGGTTGGATGATAAGT"),
    append_seq="TA",
    min_pos=29,
    max_pos=32
  )  
)


split_inputs_ <- function(input_file,
                          barcodes,
                          barcode_start=10,
                          barcode_end=15,
                          output_files,
                          max_reads=Inf,
                          fail_file=NA,
                          output_tnq_fails=FALSE,
                          
                          match_transposon=TRUE,
                          transposon="magellan",
                          tn_params=NA,
                          default_cut=NA,
                          
                          output_format="fasta",
                          output_untrimmed=FALSE,
                          collapse=T,
                          min_quality=NA) {
  
  # ============= transposon matching =============
  if (match_transposon) {
    if (is.na(transposon) && is.na(tn_params))
      stop("Either a transposon name or transposon_params must be given.")
    if (is.character(transposon) && !(transposon %in% names(transposons)))
      stop("Unknown transposon ", transposon)
  }
  if (!is.na(transposon)) {
    tn_params <- transposons[[transposon]]
  }
  if (!match_transposon && is.na(default_cut) && !output_untrimmed)
    stop("default_cut must be specified if transposon matching is off.")
  
  # ============= output file formats =============
  valid_formats <- c(fastq="fastq", fq="fastq", fasta="fasta", fa="fa")
  if (!(output_format %in% names(valid_formats)))
    stop("Unsupported output format ", output_format)
  output_format <- valid_formats[output_format]  # choose 'fastq' or 'fasta'
  if (collapse && output_format == "fastq")
    stop("FASTQ files cannot be collapsed.  Please output FASTA files.")
  
  # ============= output file names =============
  if (length(output_files) == 1 && is.null(names(output_files))) {
    # path prefix given; create names
    output_files <- paste0(check_path_ending(output_files),
                           input_file, "_",
                           names(barcodes),
                           ".", output_format)
    names(output_files) <- names(barcodes)
  } else {
    # check that all barcodes are represented
    if (!all(names(barcodes) %in% names(output_files)))
      stop("All barcodes must have an output file.")
  }
  if (output_tnq_fails) {
    fail_files <- paste0(output_files, "_", "FAILS", output_format)
    names(fail_files) <- names(barcodes)
  }
  
  # ============= output function =============
  # touch all files first
  outputf <- function(reads, filename) {
    if (collapse) {
      # add to hashes
    } else if (output_format == "fasta") {
      ShortRead::writeFasta(reads, filename, mode="a")
    } else {
      ShortRead::writeFastq(reads, filename, mode="a")
    }
  }
  
  # ============= hashtables for collapsing reads =============
  if (collapse) {
    hashes <- output_files
    if (output_fails) {
      hashes <- c(output_files, fail_files, fail_file)
    }
    build_hashes(hashes)
  }
  
  # ============= collect stats about splitting =============
  n_barcodes <- length(barcodes)
  stats <- data.frame(reads=integer(n_barcodes),
                      matched_tn=integer(n_barcodes),
                      quality=integer(n_barcodes))
  rownames(stats) <- names(barcodes)
  total_reads <- 0
  break_at_end <- FALSE
  
  streamer <- ShortRead::FastqStreamer(fastq)
  repeat {
    fq <- ShortRead::yield(streamer)
    
    # see if this yield contains enough reads to stop
    if (total_reads + length(fq) >= max_reads) {
      fq <- fq[1:(max_reads - total_reads)]
      break_at_end <- TRUE  # stop after this iteration
    }
    total_reads <- total_reads + length(fq)
    
    barcode_strings <- fq %>%
      ShortRead::narrow(start=barcode_start, end=barcode_end) %>%
      ShortRead::sread()
    assignments <- match_barcodes(barcode_strings, barcodes)
    for (code in names(barcodes)) {
      reads <- fq[assignments[ ,barcodes(code)]]
      stats[code,"reads"] <- stats[code,"reads"] + length(reads)
      
      if (match_transposon) {
        read_length <- width(reads[1])
        start_pos <- 1  # match across entire read
        n_sites <- read_length - start_pos - length(tn_params$seq) + 2
        site_counts <- integer(n_sites)
        starting.at <- start_pos:(read_length-length(tn_params$seq)+1)
        edits <- Biostrings::neditStartingAt(tn_params$seq,
                                             ShortRead::sread(reads),
                                             starting.at=starting.at)
        argmins <- apply(edits, 2, which.min)
        if (!output_trimmed) {
          cutat <- ifelse(argmins < barcode_end+2, barcode_end+1, argmins-1)
          reads %<>% ShortRead::narrow(start=barcode_end+1, end=cutat)
        }
        passing <- argmins >= tn_params$min_pos && argmins <= tn_params$max_pos
        if (output_tnq_fails) {
          outputf(reads[!passing], fail_files[code])
        }
        outputf(reads[passing], output_files[code])
        
      } else {
        reads %<>% ShortRead::narrow(start=16, end=28)
      }
      stats[code,"matched_tn"] <- stats[code,"matched_tn"] + length(reads)
      
      if (!is.na(quality)) {
        # filter reads for quality
        mean_quality <- rowMeans(as(quality(reads), "matrix"), na.rm=T)
        reads <- reads[mean_quality > quality]
      }
      stats[code,"quality"] <- stats[code,"quality"] + length(reads)
      
      keys <- reads %>% ShortRead::sread() %>% as.character() %>% paste0("TA")
      add_keys(code, keys)
    }
    
    if (dump_fails) {
      fails <- colSums(assignments) == 0
      add_keys("FAILS", fq %>% ShortRead::sread() %>% as.character())
    }
    
    # exit if we've read enough reads
    if (break_at_end) break
  }
  close(streamer)
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
                         match_transposon_sequence=TRUE, quality=NA) {
  barcodes <- tnseq$barcodes
  
  split_fastq_file <- function(input, fastq) {
    start_log_file(get_path(tnseq, dir="log", 
                            file=paste0("split_", input, ".log")))
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
      
      barcode_strings <- fq %>%
        ShortRead::narrow(start=10, end=15) %>%
        ShortRead::sread()
      assignments <- match_barcodes(barcode_strings, barcodes)
      for (code in barcodes) {
        reads <- fq[assignments[ ,code]]
        stats[code,"reads"] <- stats[code,"reads"] + length(reads)
        
        if (match_transposon_sequence) {
          # matching here
          pattern <- Biostrings::DNAString("TAACAGGTTGGATGATAAGT")
          read_length <- 51
          start_pos <- 1
          n_sites <- read_length - start_pos - length(pattern) + 2
          site_counts <- integer(n_sites)
          avg_dists <- numeric(n_sites)
          starting.at <- start_pos:(read_length-length(pattern)+1)
          edits <- Biostrings::neditStartingAt(pattern,
                                               ShortRead::sread(reads),
                                               starting.at=starting.at)
          argmins <- apply(edits, 2, which.min)
          cutat <- ifelse(argmins < 29, 17, argmins-1)
          reads %<>% ShortRead::narrow(start=16, end=cutat)
          reads <- reads[argmins >= 29 & argmins <= 32]
        } else {
          reads %<>% ShortRead::narrow(start=16, end=28)
        }
        stats[code,"matched_tn"] <- stats[code,"matched_tn"] + length(reads)
        
        if (!is.na(quality)) {
          # filter reads for quality
          mean_quality <- rowMeans(as(quality(reads), "matrix"), na.rm=T)
          reads <- reads[mean_quality > quality]
        }
        stats[code,"quality"] <- stats[code,"quality"] + length(reads)
        
        keys <- reads %>% ShortRead::sread() %>% as.character() %>% paste0("TA")
        add_keys(code, keys)
      }

      if (dump_fails) {
        fails <- colSums(assignments) == 0
        add_keys("FAILS", fq %>% ShortRead::sread() %>% as.character())
      }      
    }
    close(streamer)
    
    fileends <- paste0(input, "_", names(barcodes), ".fasta")
    files <- get_path(tnseq, dir="split", file=fileends)
    for (i in seq_along(barcodes)) {
      write_hash_to_file(barcodes[i], path.expand(files[i]))
    }
    if (dump_fails) {
      failfile <- get_path(tnseq, dir="split", 
                           file=paste0(input, "_", "FAILS.fasta"))
      write_hash_to_file("FAILS", path.expand(failfile))
    }
    
    # ==================== reporting ============================
    
    reportf("Processing file %s (%s).", input, fastq)
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
    
    unindent_report()
    set_log_file(NULL)
    
    erase_hashes(barcodes)
    if (dump_fails) erase_hashes("FAILS")
    
    rownames(stats) <- NULL
    df <- dplyr::data_frame(input=input, barcode=names(barcodes), file=files)
    return(cbind(df, stats))
  }
  
  dfs <- parallel::mclapply(seq_along(tnseq$filelog$input$file), 
                function(i) split_fastq_file(tnseq$filelog$input$input[i],
                                             tnseq$filelog$input$file[i]))
  tnseq$filelog$split <- do.call(rbind, dfs)
  
  return(tnseq)
}

