
#' @export
filter_insertions <- function(tnseq, drop_t1_zeros=T, drop_t2_zeros=F, 
                              min_total_reads=15,
                              trim_head_frac=0.0, trim_tail_frac=0.1) {
  passing <- !logical(nrow(tnseq$insertions))
  if (drop_t1_zeros) passing <- passing & tnseq$insertions$reads1 > 0
  if (drop_t2_zeros) passing <- passing & tnseq$insertions$reads2 > 0
  if (!is.na(min_total_reads))
    passing <- passing & with(tnseq$insertions, reads1 + reads2 >= min_total_reads)
  if (!is.na(trim_head_frac))
    passing <- passing & tnseq$insertions$relpos > trim_head_frac
  if (!is.na(trim_tail_frac))
    passing <- passing & tnseq$insertions$relpos < (1 - trim_tail_frac)
  
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
                           log(expansion*f2/f1) / log(expansion*(1-f2)/(1-f1))))
  return(tnseq)
}
