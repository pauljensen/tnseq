
#' @export
split_inputs_rnatagseq <- function(rtseq, ...) {
  rtseq %>% split_inputs(match_transposon=F, 
                         barcode_start=1, 
                         barcode_end=6, 
                         default_cut=40, 
                         collapse=F, 
                         output_format="fasta",
                         ...)
}

#' @export
reads_to_fpkm_ <- function(rtseq, mapfile) {
  cnames <- c("read", "strand", "genome", "pos", "seq", "score", "v7", "v8")
  reads <- read.table(mapfile, sep="\t", col.names=cnames, stringsAsFactors=F)
  #print(reads$genome[1])
  features <- rtseq$features[[reads$genome[1]]]
  hits <- IRanges::findOverlaps(IRanges::IRanges(start=reads$pos, width=1),
                                GenomicRanges::ranges(features))
  hitsdf <- as.data.frame(table(IRanges::subjectHits(hits)))
  counts <- numeric(length(features))
  names(counts) <- names(features)
  counts[as.integer(levels(hitsdf$Var1))] <- hitsdf$Freq
  fpkm <- counts / (features$gene_length / 1000)
  fpkm <- fpkm / (sum(fpkm) / 1000000)
  return(list(counts=counts, fpkm=fpkm))
}

#' @export
reads_to_fpkm <- function(rtseq) {
  rtseq$counts <- list()
  rtseq$fpkm <- list()
  for (strain in unique(rtseq$samples$strain)) {
    print(strain)
    samples <- rtseq$samples[rtseq$samples$strain == strain, ]
    sample_names <- paste(samples$strain, samples$condition, samples$library)
    genome <- samples$genome[1]
    genes <- names(rtseq$features[[genome]])
    M_count <- matrix(nrow=length(genes), ncol=length(sample_names),
                      dimnames=list(genes, sample_names))
    M_fpkm <- M_count
    for (i in seq_along(sample_names)) {
      values <- reads_to_fpkm_(rtseq, substr(samples$mapfile[i], 1, nchar(samples$mapfile[i])-1))
      M_count[ ,i] <- values$counts
      M_fpkm[ ,i] <- values$fpkm
    }
    rtseq$counts[[genome]] <- M_count
    rtseq$fpkm[[genome]] <- M_fpkm
  }
  
  return(rtseq)
}

#' @export
process_rnatagseq <- function(dirpath, ...) {
  load_tnseq_experiment(dirpath) %>% 
    load_genomes() %>% 
    split_inputs_rnatagseq(...) %>% 
    map_reads() %>% 
    reads_to_fpkm()
}

