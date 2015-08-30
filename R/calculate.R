
#' @export
filter_insertions <- function(tnseq, drop_t1_zeros=T, drop_t2_zeros=T, 
                              min_total_reads=15,
                              trim_head_frac=0.0, trim_tail_frac=0.1) {
  passing <- !logical(nrow(tnseq$insertions))
  if (drop_t1_zeros) passing <- passing & tnseq$insertions$reads1 > 0
  if (drop_t2_zeros) passing <- passing & tnseq$insertions$reads2 > 0
  if (!is.na(min_total_reads))
    passing <- passing & with(tnseq$insertions, reads1 + reads2 >= min_total_reads)
  if (!is.na(trim_head_frac))
    passing <- passing & ifelse(is.na(tnseq$insertions$relpos),
                                TRUE,
                                tnseq$insertions$relpos > trim_head_frac)
  if (!is.na(trim_tail_frac))
    passing <- passing & ifelse(is.na(tnseq$insertions$relpos),
                                TRUE,
                                tnseq$insertions$relpos < (1 - trim_tail_frac))
  
  tnseq$insertions <- tnseq$insertions[passing, ]
  return(tnseq)
}

#' @export
calculate_fitness <- function(tnseq) {
  tnseq$insertions %<>%
    dplyr::group_by_("strain", "library", "condition") %>%
    dplyr::mutate(f1=reads1/sum(reads1),
                  f2=reads2/sum(reads2),
                  W=ifelse(reads2 == 0,
                           0,
                           log(expansion*f2/f1) / log(expansion*(1-f2)/(1-f1))),
                  method="CALCULATED")
  return(tnseq)
}

weighted_sd <- function(x, w) {
  sqrt(weighted.mean((x - mean(x))^2, w))
}

#' @export
aggregate <- function(tnseq, grouping=c("strain", "library", "condition", "gene"),
                      max_weight=50) {
  
  weight <- function(reads1, reads2) ifelse(reads1 + reads2 > max_weight,
                                            max_weight,
                                            reads1 + reads2)
  
  tnseq$aggregate <- tnseq$insertions %>%
    dplyr::group_by_(.dots=grouping) %>%
    dplyr::summarize(
      insertions=n(),
      zero_t2=sum(reads2 == 0),
      reads_t1=sum(reads1),
      reads_t2=sum(reads2),
      
      fitness=mean(W),
      stdev=sd(W),
      stderr=sd(W)/sqrt(n()),
      
      weighted_fitness=weighted.mean(W, weight(reads1, reads2)),
      weighted_stdev=weighted_sd(W, weight(reads1, reads2)),
      weighted_stderr=weighted_stdev/sqrt(n()),
      
      genome=paste(unique(genome), collapse=", "),
      method=unique(method)
    ) 
  
  return(tnseq)
}

#' @export
add_missing_genes <- function(tnseq) {
  add_miss <- function(df) {
    all_genes <- names(tnseq$features[[unique(df$genome)]])
    missing <- setdiff(all_genes, df$gene)
    n <- length(missing)
    newdf <- df[rep(1, n), ]
    newdf$gene <- missing
    newdf$insertions <- 0
    newdf$zero_t2 <- 0
    newdf$reads_t1 <- NA
    newdf$reads_t2 <- NA
    newdf$fitness <- 0
    newdf$stdev <- NA
    newdf$stderr <- NA
    newdf$weighted_fitness <- 0
    newdf$weighted_stdev <- NA
    newdf$weighted_stderr <- NA
    newdf$method <- "MISSING"
    
    rbind(df, newdf)
  }
  
  tnseq$aggregate %<>% 
    dplyr::group_by_("strain", "condition") %>%
    dplyr::do(add_miss(.)) %>%
    dplyr::arrange(gene) %>%
    dplyr::ungroup()

  return(tnseq)
}
