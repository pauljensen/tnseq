
plot_fitness_aux <- function(insertions, features, start, end) {
  plot(c(start, end), c(0, 0), type="l", ylim=c(-0.2, 1.0),
       ylab="fitness", xlab="genome position")
  abline(h=1.0, lty="solid", col="gray")
  abline(h=0.5, lty="dashed", col="gray")
  
  starts <- GenomicRanges::start(features)
  ends <- GenomicRanges::end(features)
  gene_names <- names(features)
  for (i in seq_along(starts)) {
    if (as.character(strand(features)[i]) == "+") {
      xs <- c(starts[i], ends[i], max(starts[i], ends[i]-75), starts[i])
    } else {
      xs <- c(starts[i], ends[i], ends[i], min(starts[i]+75, ends[i]))
    }
    polygon(xs, c(0, 0, -0.2, -0.2), col="lightblue")
    text(x=mean(c(starts[i], ends[i])), y=-0.1, gene_names[i])
  }
  
  for (i in 1:nrow(insertions)) {
    lines(c(insertions$pos[i], insertions$pos[i]), 
           c(0, insertions$W[i]))
    points(insertions$pos[i], 0)
  }
}

plot_fitness <- function(tnseq, strain, condition, gene, 
                         lib="all", margin=1000) {
  strn <- strain
  cond <- condition
  insertions <- tnseq$insertions %>%
    dplyr::filter(strain==strn, condition==cond)
  if (lib != "all") {
    insertions %<>% dplyr::filter(library==lib)
  }
  
  genome <- insertions$genome[1]
  
  if (is.null(gene) || !(gene %in% names(tnseq$features[[genome]]))) {
    gene <- names(tnseq$features[[genome]])[1]
  }
  
  plot_start <- start(tnseq$features[[genome]][gene]) - margin
  plot_end <- end(tnseq$features[[genome]][gene]) + margin
  
  insertions %<>% dplyr::filter(pos >= plot_start, pos <= plot_end)
  
  features <- IRanges::IRanges(start=plot_start, end=plot_end) %>%
    IRanges::findOverlaps(ranges(tnseq$features[[genome]])) %>% 
    IRanges::subjectHits() %>%
    tnseq$features[[genome]][.]
  
  plot_fitness_aux(insertions, features, plot_start, plot_end)
}

plot_fitness(tnseq, "316", "SDMM", "SP_0771", margin=1000)
