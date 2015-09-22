
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

mcols_defaults <- list(
  name="gene", 
  product="product", 
  gene_length=f(x, length(x$sequence)),
  tn_sites=f(x, get_tn_count(x$sequence, "TA"))
)

create_features_table <- function(gbk, feature="CDS", key="locus_tag",
                                  mcols=mcols_defaults) {
  features <- gbk$features[[feature]]
  n <- length(features)
  
  locs <- lapply(features, f(x, get_startend(x$location)))
  genes <- character(n)
  meta <- vector(mode="list", length=n)
  for (i in seq_along(features)) {
    # add the DNA sequence to the feature in case any of the mcols function
    # need it
    features[[i]]$sequence <- genbank::extract_sequence(features[[i]]$location,
                                                        gbk$sequence)
    
    # use key to find the gene name
    if (key %in% names(features[[i]])) {
      genes[i] <- features[[i]][[key]]
    } else {
      # use generic number feature (e.g. "CDS--1")
      genes[i] <- paste0(feature, "--", i)
    }
    
    # extract the meta-column data
    meta[[i]] <- list()
    for (name in names(mcols)) {
      if (is.character(mcols[[name]])) {
        meta[[i]][[name]] <- ifelse(is.null(features[[i]][[mcols[[name]]]]),
                                    NA,
                                    features[[i]][[mcols[[name]]]])
      } else {
        # mcols[[name]] is a function that returns the meta-column
        meta[[i]][[name]] <- mcols[[name]](features[[i]])
      }
    }
  }
  
  # make sure we have complete location data for each gene
  to_keep <- sapply(locs, f(x, !any(is.na(x))))
  locs <- locs[to_keep]
  genes <- genes[to_keep]
  meta <- meta[to_keep]
  
  # build the final GenomicRanges structure
  iranges <- IRanges::IRanges(start=sapply(locs, f(x, x$start)),
                              end=sapply(locs, f(x, x$end)))
  gr <- GenomicRanges::GRanges(ranges=iranges, 
                               strand=sapply(locs, f(x, x$strand)),
                               seqnames="chr1")
  names(gr) <- genes
  if (length(mcols) > 0) {
    df <- list()
    for (name in names(meta[[1]])) {
      df[[name]] <- sapply(meta, f(x, x[[name]]))
    }
    GenomicRanges::mcols(gr) <- as.data.frame(df, stringsAsFactors=F)
  }
  
  return(gr)
}

#' @export
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
    features_path <- paste0(fullpath, "_features.Rdata")
    
    use_cache <- cached && file.exists(fasta_path) && file.exists(features_path)
    if (!use_cache) {
      reportf("Parsing GenBank file %s.", fullpath)
      gbk <- genbank::parse_genbank(fullpath)
      features <- create_features_table(gbk)
      GenomeInfoDb::seqlevels(features) <- genome
      tnseq$features[[genome]] <- features
      tnseq$genome_sequences[[genome]] <- gbk$sequence
      reportf("Found %i bases, %i genes, %i transposon sites",
              length(tnseq$genomic_sequences[[genome]]),
              length(tnseq$features[[genome]]),
              sum(tnseq$features[[genome]]$sites))
      
      reportf("Writing features table to %s.", features_path)
      save(features, file=features_path)
      
      reportf("Writing FASTA sequence file to %s.", fasta_path)
      writeLines(c(paste0(">", genome), 
                   as.character(tnseq$genome_sequences[[genome]])),
                 fasta_path)
    } else {
      # load genome sequence
      reportf("Loading genome sequence from %s.", fasta_path)
      seq <- readLines(fasta_path)[2]
      tnseq$genome_sequences[[genome]] <- Biostrings::DNAString(seq)
      
      # load features table
      reportf("Loading features table from %s.", features_path)
      load(features_path)
      tnseq$features[[genome]] <- features
    }
    unindent_report()
  }
  unindent_report()
  return(tnseq)
}

#' @export
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

#' @export
add_gene_metadata_ <- function(df, features) {
  features <- GenomicRanges::mcols(do.call(c, unname(features)), use.names=T)
  dplyr::bind_cols(df,
                   features[df$gene, ] %>% 
                     as.list() %>% 
                     dplyr::as_data_frame())
}

#' @export
add_gene_metadata <- function(tnseq, slots=c("insertions", "aggregate")) {
  for (slot in slots) {
    if (!is.null(tnseq[[slot]])) {
      tnseq[[slot]] <- add_gene_metadata_(tnseq[[slot]], tnseq$features)
    }
  }
  return(tnseq)
}

#' @export
remove_nongenic_ <- function(df) {
  df[!is.na(df$gene), ]
}

#' @export
remove_nongenic <- function(tnseq, slots=c("insertions", "aggregate")) {
  for (slot in slots) {
    if (!is.null(tnseq[[slot]])) {
      tnseq[[slot]] <- remove_nongenic_(tnseq[[slot]])
    }
  }
  return(tnseq)
}