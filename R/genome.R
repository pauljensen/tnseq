
posmatch <- function(str, start, end, pattern) {
  if (nchar(str) < end) {
    return(FALSE)
  } else {
    return(str_detect(substr(str, start, end), pattern))
  }
}

has_header <- function(lines) {
  return(posmatch(lines[1], 21, 46, "Genetic Sequence Data Bank"))
}

parse_header <- function(lines) {
  return(list(
    filename = str_trim(substr(lines[1], 1, 20)),  # format says 9, but some are longer
    date = str_trim(lines[2]),
    release = str_trim(substr(lines[4], 48, 52)),
    title = str_trim(lines[6]),
    loci = strtoi(substr(lines[8], 1, 8)),
    bases = strtoi(substr(lines[8], 16, 26)),
    reports = strtoi(substr(lines[8], 40, 47))
  ))
}

break_by_keywords <- function(lines) {
  keywords <- list()
  start <- 1
  i <- 1
  while (i < length(lines)) {
    if (str_detect(lines[i+1], "^\\S")) {
      keywords[[length(keywords)+1]] <- lines[start:i]
      names(keywords[[length(keywords)]]) <- NULL
      start <- i + 1
    }
    i <- i + 1
  }
  # ignore last entry; should be "//"
  return(keywords)
}

group_list <- function(l, keyfun, valfun=function(x) x) {
  keys <- sapply(l, keyfun)
  grouped <- list()
  for (key in unique(keys)) {
    grouped[[key]] <- lapply(l[key == keys], valfun)
  }
  return(grouped)
}

gen_group_parser <- function(group, getsubkey, getval, valsep) {
  issubkey <- function(...) nchar(getsubkey(...)) > 0
  
  ret <- list()
  i <- 1
  while (i <= length(group)) {
    # should be a subkey
    if (!issubkey(group[i])) stop("subkey expected at: ", group[i])
    sk <- getsubkey(group[i])
    if (!(sk %in% names(ret))) ret[[sk]] <- character(0)
    val <- getval(group[i])
    while (i < length(group) && !issubkey(group[i+1])) {
      i <- i + 1
      val <- paste(val, getval(group[i]), sep=valsep)
    }
    i <- i + 1
    ret[[sk]] <- c(ret[[sk]], val)
  }
  
  # set name and value fields using first subkey
  first_sk <- getsubkey(group[1])
  ret$value <- ret[[first_sk]]
  ret[[first_sk]] <- NULL
  ret$name <- first_sk
  
  return(ret)
}

parse_keyword_group <- function(...) {
  gen_group_parser(..., 
                   getsubkey=function(line) str_trim(substr(line, 1, 12)),
                   getval=function(line) substr(line, 13, 80),
                   valsep=" ")
}

parse_feature_group <- function(group) {
  getsubkey <- function(l) {
    l <- str_trim(l)
    if (str_detect(l, "^/")) {
      return(str_match(l, "^/(\\w+)=")[1,2])
    } else {
      return("")
    }
  }
  getval <- function(l) {
    l <- str_trim(l)
    if (str_detect(l, "^/")) {
      return(str_match(l, "^/(\\w+)=(.+)")[1,3])
    } else {
      return("")
    }
  }
  gen_group_parser(group[-1], getsubkey=getsubkey, getval=getval, valsep="")
}

pattern_idx <- "(\\d+)"
pattern_loc <- paste0("([<>])?", pattern_idx)
fullstr <- function(...) paste0("^", ..., "$")
pattern_single <- fullstr(pattern_idx)
pattern_site <- fullstr(pattern_loc, "\\^", pattern_loc)
pattern_span <- fullstr(pattern_loc, "\\.\\.", pattern_loc)
pattern_range <- fullstr(pattern_loc, "\\.", pattern_loc)
pattern_op <- fullstr("(\\w+)\\((.+)\\)")

simple_op <- function(operator, start="", end="", start_open="", end_open="") {
  list(operator=operator, 
       start=as.integer(start), 
       end=as.integer(end), 
       start_open=ifelse(nchar(start_open) > 0, start_open, NA), 
       end_open=ifelse(nchar(end_open) > 0, end_open,  NA))
}

parse_location <- function(l) {
  if (str_detect(l, pattern_single)) {
    # single base
    m <- str_match(l, pattern_single)
    return(simple_op("base", start=m[1,2], end=m[1,2]))
    
  } else if (str_detect(l, pattern_site)) {
    # site between bases
    m <- str_match(l, pattern_site)
    return(simple_op("site", start_open=m[1,2], start=m[1,3], 
                     end_open=m[1,4], end=m[1,5]))
    
  } else if (str_detect(l, pattern_span)) {
    # base span
    m <- str_match(l, pattern_span)
    return(simple_op("span", start_open=m[1,2], start=m[1,3], 
                     end_open=m[1,4], end=m[1,5]))
    
  } else if (str_detect(l, pattern_range)) {
    # base range
    m <- str_match(l, pattern_range)
    return(simple_op("range", start_open=m[1,2], start=m[1,3], 
                     end_open=m[1,4], end=m[1,5]))
    
  } else if (str_detect(l, pattern_op)) {
    # operator
    m <- str_match(l, pattern_op)
    arg_strs <- str_split(m[1,3], ",\\s*")[[1]]
    return(list(operator=m[1,2], args=lapply(arg_strs, parse_location)))
  }
}

parse_genbank <- function(file) {
  lines <- readLines(file)
  ret <- list()
  if (has_header(lines)) {
    ret$HEADER <- parse_header(lines)
    lines <- lines[-(1:9)]
  } else {
    ret$HEADER <- NULL
  }
  
  groups <- break_by_keywords(lines)
  get_keyword <- function(x) str_split(x[1], "\\s+")[[1]][1]
  keywords <- sapply(groups, get_keyword)
  features <- groups[keywords == "FEATURES"][[1]]
  origin <- groups[keywords == "ORIGIN"][[1]]
  groups <- groups[!(keywords %in% c("FEATURES", "ORIGIN"))]
  
  ret <- lapply(groups, parse_keyword_group)
  ret <- group_list(ret, keyfun=function(x) x$name, 
                    valfun=function(x) x[names(x) != "name"])
  
  # singlets appear only once per file and have no subkeys
  singlets <- c("KEYWORDS", "VERSION", "ACCESSION", "DEFINITION", "LOCUS")
  for (singlet in singlets) {
    ret[[singlet]] <- ret[[singlet]][[1]]$value
  }
  
  # onesies appear only once per file but have subkeys
  onesies <- c("SOURCE")
  for (onesie in onesies) {
    ret[[onesie]] <- ret[[onesie]][[1]]
  }

  ret$features <- c(sapply(features[-1], function(x) substr(x, 6, 80)), "//") %>%
    break_by_keywords() %>% lapply(parse_feature_group)
  
  return(ret)
}

lines <- c(
  "LOCUS       SCU49845     5028 bp    DNA             PLN       21-JUN-1999",
  "  ORGANISM  Saccharomyces cerevisiae",
  "            Eukaryota; Fungi; Ascomycota; Saccharomycotina; Saccharomycetes;"
)



