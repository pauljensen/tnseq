
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

load_genome <- function(genome, feature="gene", name_key="locus_tag", tn_site="TA") {
  gbk <- genbank::parse_genbank(genome)
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

load_genomes <- function(tnseq) {
  tnseq$genomes <- list()
  report("Loading genome files.")
  indent_report()
  for (genome in unique(tnseq$samples$genome)) {
    fullpath <- get_path(tnseq, dir="genomes", file=genome)
    reportf("%s: %s", genome, fullpath)
    if (!file.exists(fullpath)) {
      stop("Could not find genome file ", fullpath)
    }
    tnseq$genomes[[genome]] <- load_genome(fullpath)
    indent_report()
    reportf("%i bases, %i genes, %i transposon sites",
            seqlengths(tnseq$genomes[[genome]])[1],
            length(tnseq$genomes[[genome]]),
            sum(tnseq$genomes[[genome]]$sites))
    unindent_report()
  }
  unindent_report()
  return(tnseq)
}
