
get_startend <- function(location) {
  if (location$operator == "span") {
    l <- list(start=location$start, end=location$end, strand="+")
  } else if (location$operator == "complement") {
    l <- get_startend(location$args[[1]])
    l$strand <- "-"
  } else {
    l <- list(start=NA, end=NA, strand=NA)
  }
  return(l)
}

get_tn_count <- function(sequence, tn_site) {
  freq <- Biostrings::oligonucleotideFrequency(sequence, nchar(tn_site))
  rc_tn_site <- tn_site %>% 
    Biostrings::DNAString() %>% 
    ShortRead::reverseComplement() %>% 
    as.character()
  return(unname(freq[tn_site] + freq[rc_tn_site]))
}

create_features_table <- function(gbk, feature="gene", name_key="locus_tag", 
                                  tn_site="TA") {
  features <- gbk$features[[feature]]
  n <- length(features)
  locs <- lapply(features, f(x, get_startend(x$location)))
  for (i in seq_along(features)) {
    sequence <- genbank::extract_sequence(features[[i]]$location, gbk$sequence)
    locs[[i]]$tn_count <- get_tn_count(sequence, tn_site)
    if (name_key %in% names(features[[i]])) {
      locs[[i]]$name <- features[[i]][[name_key]]
    } else {
      locs[[i]]$name <- paste0(feature, "--", i)
    }
  }
  
  locs <- Filter(f(l, !any(is.na(l))), locs)
  iranges <- IRanges::IRanges(start=sapply(locs, f(x, x$start)),
                              end=sapply(locs, f(x, x$end)))
  gr <- GenomicRanges::GRanges(ranges=iranges, 
                               strand=sapply(locs, f(x, x$strand)),
                               sites=sapply(locs, f(x, x$tn_count)),
                               seqnames="chr1")
  names(gr) <- sapply(locs, f(x, x$name))
  seqlengths(gr) <- length(gbk$sequence)
  return(gr)
}

load_genomes <- function(tnseq, cached=T) {
  tnseq$genome_sequences <- list()
  tnseq$features <- list()
  report("Loading genome files.")
  indent_report()
  for (genome in tnseq$genomes) {
    fullpath <- get_path(tnseq, dir="genomes", file=genome)
    reportf("%s: %s", genome, fullpath)
    indent_report()
    if (!file.exists(fullpath)) {
      stop("Could not find genome file ", fullpath)
    }
    
    fasta_path <- paste0(fullpath, ".fasta")
    features_path <- paste0(fullpath, ".features")
    
    use_cache <- cached && file.exists(fasta_path) && file.exists(features_path)
    if (!use_cache) {
      reportf("Parsing GenBank file %s.", fullpath)
      gbk <- genbank::parse_genbank(fullpath)
      tnseq$features[[genome]] <- create_features_table(gbk)
      tnseq$genome_sequences[[genome]] <- gbk$sequence
      reportf("Found %i bases, %i genes, %i transposon sites",
              seqlengths(tnseq$features[[genome]])[1],
              length(tnseq$features[[genome]]),
              sum(tnseq$features[[genome]]$sites))
      
      reportf("Writing features table to %s.", features_path)
      features <- tnseq$features[[genome]]
      readr::write_csv(data.frame(names=names(features),
                                  seqnames=seqnames(features),
                                  start=start(features),
                                  end=end(features),
                                  strand=strand(features),
                                  sites=features$sites,
                                  stringsAsFactors=F),
                       path=features_path)
      
      reportf("Writing FASTA sequence file to %s.", fasta_path)
      writeLines(c(paste0(">", genome), 
                   as.character(tnseq$genome_sequences[[genome]])),
                 fasta_path)
    } else {
      # load genome sequence
      reportf("Loading genome sequence from %s.", fasta_path)
      seq <-readLines(fasta_path)[2]
      tnseq$genome_sequences[[genome]] <- Biostrings::DNAString(seq)
      
      # load features table
      reportf("Loading features table from %s.", features_path)
      feats <- readr::read_csv(features_path, col_names=T, col_types="cciici")
      iranges <- IRanges::IRanges(start=feats$start,
                                  end=feats$end)
      gr <- GenomicRanges::GRanges(ranges=iranges, 
                                   strand=feats$strand,
                                   sites=feats$sites,
                                   seqnames=feats$seqnames)
      names(gr) <- feats$names
      seqlengths(gr) <- length(tnseq$genome_sequences[[genome]])
      tnseq$features[[genome]] <- gr
    }
    unindent_report()
  }
  unindent_report()
  return(tnseq)
}

create_bowtie_indexes <- function(tnseq) {
  reportf("Creating Bowtie index files.")
  indent_report()
  for (genome in tnseq$genomes) {
    genome_path <- get_path(tnseq, dir="genomes", file=genome)
    # for bowtie1, running bowtie-build -c with the sequence on the command line
    # does not produce the .rev index files. Need to use the fasta file.
    fasta_path <- paste0(genome_path, ".fasta")
    logfile <- get_path(tnseq, dir="log", file=paste0("index_", genome, ".log"))
    reportf("Building index %s.", genome_path)
    cmd <- paste(BOWTIE_BUILD_CMD, "-f",
                 fasta_path, 
                 genome_path,
                 ">", logfile)
    system(cmd)
  }
  unindent_report()
}

