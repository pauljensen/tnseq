
get_genome_feature_table <- function(tnseq) {
  dfs <- lapply(tnseq$features, function(x) data.frame(gene=names(x),
                                                       name=x$name,
                                                       product=x$product))
  df <- do.call(rbind, dfs)
  rownames(df) <- df$gene
  df$gene <- NULL
  return(df)
}

#' @export
generate_reports <- function(tnseq) {
  lookup <- get_genome_feature_table(tnseq)
  
  readr::write_csv(tnseq$insertions, 
                   path=get_path(tnseq, dir="reports", file="insertions.csv"))

  writer <- function(df) {
    file <- paste(df$strain, df$library, df$condition)
    file <- paste0(file, ".csv")
    df <- df[!is.na(df$gene), ]
    df <- cbind(df, lookup[df$gene, ])
    readr::write_csv(df, path=get_path(tnseq, dir="reports", file=file)[1])
    return(df)
  }
  
  grouping <- c("strain", "library", "condition", "gene")
  tnseq %<>% aggregate(grouping=grouping) %>% add_missing_genes()
  tnseq$aggregate %>%
    dplyr::group_by_(.dots=setdiff(grouping, "gene")) %>%
    dplyr::do(writer(.)) %>%
    dplyr::ungroup()
  
  grouping <- c("strain", "condition", "gene")
  tnseq %<>% aggregate(grouping=grouping) %>% add_missing_genes()
  tnseq$aggregate %>%
    dplyr::group_by_(.dots=setdiff(grouping, "gene")) %>%
    dplyr::do(writer(.)) %>%
    dplyr::ungroup()
  
  return(NULL)
}

#' @export
quick_summary <- function(tnseq) {
  tnseq$aggregate %>%
    dplyr::group_by(strain, condition) %>%
    dplyr::summarize(insertions=sum(insertions),
                     mean_stdev=mean(stdev, na.rm=T),
                     mean_stderr=mean(stderr, na.rm=T),
                     total_genes=sum(!is.na(gene)),
                     missing_genes=sum(method=="MISSING"))
}