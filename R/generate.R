
#' @export
generate_reports <- function(tnseq) {
  readr::write_csv(tnseq$insertions, 
                   path=get_path(tnseq, dir="reports", file="insertions.csv"))

  writer <- function(df) {
    file <- paste(df$strain, df$library, df$condition)
    file <- paste0(file, ".csv")
    readr::write_csv(df, path=get_path(tnseq, dir="reports", file=file)[1])
    return(df)
  }
  
  grouping <- c("strain", "library", "condition", "gene")
  tnseq %<>% aggregate(grouping=grouping)
  tnseq$aggregate %>%
    dplyr::group_by_(.dots=setdiff(grouping, "gene")) %>%
    dplyr::do(writer(.)) %>%
    dplyr::ungroup()
  
  grouping <- c("strain", "condition", "gene")
  tnseq %<>% aggregate(grouping=grouping)
  tnseq$aggregate %>%
    dplyr::group_by_(.dots=setdiff(grouping, "gene")) %>%
    dplyr::do(writer(.)) %>%
    dplyr::ungroup()
}
