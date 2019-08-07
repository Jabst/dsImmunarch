parse_mixcr <- function (.dataframe) {
  
  IMMCOL = new.env()
  
  IMMCOL$count  = "Clones"
  IMMCOL$prop   = "Proportion"
  IMMCOL$cdr3nt = "CDR3.nt"
  IMMCOL$cdr3aa = "CDR3.aa"
  IMMCOL$v      = "V.name"
  IMMCOL$d      = "D.name"
  IMMCOL$j      = "J.name"
  IMMCOL$ve     = "V.end"
  IMMCOL$ds     = "D.start"
  IMMCOL$de     = "D.end"
  IMMCOL$js     = "J.start"
  IMMCOL$vnj    = "VJ.ins"
  IMMCOL$vnd    = "VD.ins"
  IMMCOL$dnj    = "DJ.ins"
  IMMCOL$seq    = "Sequence"
  IMMCOL$order  = c(IMMCOL$count, IMMCOL$prop, IMMCOL$cdr3nt, IMMCOL$cdr3aa,
                    IMMCOL$v, IMMCOL$d, IMMCOL$j,
                    IMMCOL$ve, IMMCOL$ds, IMMCOL$de, IMMCOL$js,
                    IMMCOL$vnj, IMMCOL$vnd, IMMCOL$dnj, IMMCOL$seq)
  IMMCOL$type   = c("numeric", "numeric", "character", "character",
                    "character", "character", "character",
                    "integer", "integer", "integer", "integer",
                    "integer", "integer", "integer", "character")
  
                    
  
  require(stringr)
  require(immunarch)
  fix.alleles <- function (.data) {
    .data[[IMMCOL$v]] <- gsub("[*][[:digit:]]*", "", .data[[IMMCOL$v]])
    .data[[IMMCOL$d]] <- gsub("[*][[:digit:]]*", "", .data[[IMMCOL$d]])
    .data[[IMMCOL$j]] <- gsub("[*][[:digit:]]*", "", .data[[IMMCOL$j]])
    .data
  }
  
  
  .count <- 'clonecount'
  .sep = '\t'
  .vend <- "allvalignments"
  .jstart <- "alljalignments"
  .dalignments <- "alldalignments"
  .vd.insertions <- "VD.insertions"
  .dj.insertions <- "DJ.insertions"
  .total.insertions <- "Total.insertions"
  
  table.colnames <- tolower(colnames(.dataframe))
  table.colnames <- gsub(".", "", table.colnames, fixed = T)
  
  # Columns of different MiXCR formats
  # Clone count - Clonal sequence(s) - N. Seq. CDR3
  # cloneCount - clonalSequence - nSeqCDR3
  # cloneCount - targetSequences - nSeqImputedCDR3
  # cloneCount - targetSequences - nSeqCDR3
  if ("targetsequences" %in% table.colnames) {
    if ('nseqimputedcdr3' %in% table.colnames) {
      .nuc.seq <- 'nseqimputedcdr3'
    } else {
      .nuc.seq <- 'nseqcdr3'
    }
    
    .big.seq <- 'targetsequences'
  } else {
    .nuc.seq <- 'nseqcdr3'
    
    if ("clonalsequences" %in% table.colnames) {
      .big.seq <- 'clonalsequences'
    } else if ("clonalsequence" %in% table.colnames) {
      .big.seq <- 'clonalsequence'
    } else {
      .big.seq = NA
    }
  }
  
  if (!("allvalignments" %in% table.colnames)) {
    if ("allvalignment" %in% table.colnames) {
      .vend = "allvalignment"
    } else {
      .vend = NA
    }
  }
  if (!("alldalignments" %in% table.colnames)) {
    if ("alldalignment" %in% table.colnames) {
      .dalignments = "alldalignment"
    } else {
      .dalignments = NA
    }
  }
  if (!("alljalignments" %in% table.colnames)) {
    if ("alljalignment" %in% table.colnames) {
      .jstart = "alljalignment"
    } else {
      .jstart = NA
    }
  }
  
  if ("bestvhit" %in% table.colnames) {
    .vgenes <- 'bestvhit'
  } else if ('allvhits' %in% table.colnames) {
    .vgenes <- 'allvhits'
  } else if ('vhits' %in% table.colnames) {
    .vgenes <- 'vhits'
  } else if ('allvhitswithscore' %in% table.colnames) {
    .vgenes <- 'allvhitswithscore'
  } else {
    cat("Error: can't find a column with V genes\n")
  }
  
  if ("bestjhit" %in% table.colnames) {
    .jgenes <- 'bestjhit'
  } else if ('alljhits' %in% table.colnames) {
    .jgenes <- 'alljhits'
  } else if ('jhits' %in% table.colnames) {
    .jgenes <- 'jhits'
  } else if ('alljhitswithscore' %in% table.colnames) {
    .jgenes <- 'alljhitswithscore'
  } else {
    cat("Error: can't find a column with J genes\n")
  }
  
  if ("bestdhit" %in% table.colnames) {
    .dgenes <- 'bestdhit'
  } else if ('alldhits' %in% table.colnames) {
    .dgenes <- 'alldhits'
  } else if ('dhits' %in% table.colnames) {
    .dgenes <- 'dhits'
  } else if ('alldhitswithscore' %in% table.colnames) {
    .dgenes <- 'alldhitswithscore'
  } else {
    cat("Error: can't find a column with D genes\n")
  }
  
  df <- .dataframe
  
  #
  # return NULL if there is no clonotypes in the data frame
  #
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  names(df) = make.names(names(df))
  names(df) <- tolower(gsub(".", "", names(df), fixed = T))
  names(df) <- str_replace_all(names(df), " ", "")
  
  # check for VJ or VDJ recombination
  # VJ / VDJ / Undeterm
  recomb_type = "Undeterm"
  if (sum(substr(head(df)[[.vgenes]], 1, 4) %in% c("TCRA", "TRAV", "TRGV", "IGKV", "IGLV"))) {
    recomb_type = "VJ"
  } else if (sum(substr(head(df)[[.vgenes]], 1, 4) %in% c("TCRB", "TRBV", "TRDV", "IGHV"))) {
    recomb_type = "VDJ"
  }
  
  if (!is.na(.vend) && !is.na(.jstart)) {
    .vd.insertions <- "VD.insertions"
    df$VD.insertions <- -1
    if (recomb_type == "VJ") {
      df$VD.insertions <- -1
    } else if (recomb_type == "VDJ") {
      logic <- sapply(strsplit(df[[.dalignments]], "|", T, F, T), length) >= 4 &
        sapply(strsplit(df[[.vend]], "|", T, F, T), length) >= 5
      df$VD.insertions[logic] <-
        as.numeric(sapply(strsplit(df[[.dalignments]][logic], "|", T, F, T), "[[", 4)) -
        as.numeric(sapply(strsplit(df[[.vend]][logic], "|", T, F, T), "[[", 5)) - 1
    }
    
    .dj.insertions <- "DJ.insertions"
    df$DJ.insertions <- -1
    if (recomb_type == "VJ") {
      df$DJ.insertions <- -1
    } else if (recomb_type == "VDJ") {
      logic <- sapply(strsplit(df[[.jstart]], "|", T, F, T), length) >= 4 &
        sapply(strsplit(df[[.dalignments]], "|", T, F, T), length) >= 5
      df$DJ.insertions[logic] <-
        as.numeric(sapply(strsplit(df[[.jstart]][logic], "|", T, F, T), "[[", 4)) -
        as.numeric(sapply(strsplit(df[[.dalignments]][logic], "|", T, F, T), "[[", 5)) - 1
    }
    
    # VJ.insertions
    logic <- (sapply(strsplit(df[[.vend]], "|", T, F, T), length) > 4) & (sapply(strsplit(df[[.jstart]], "|", T, F, T), length) >= 4)
    .total.insertions <- "Total.insertions"
    if (recomb_type == "VJ") {
      df$Total.insertions <- NA
      if (length(which(logic)) > 0) {
        df$Total.insertions[logic] <-
          as.numeric(sapply(strsplit(df[[.jstart]][logic], "|", T, F, T), "[[", 4)) - as.numeric(sapply(strsplit(df[[.vend]][logic], "|", T, F, T), "[[", 5)) - 1
      }
    } else if (recomb_type == "VDJ") {
      df$Total.insertions <- df[[.vd.insertions]] + df[[.dj.insertions]]
    } else {
      df$Total.insertions <- NA
    }
    df$Total.insertions[df$Total.insertions < 0] <- -1
    
    df$V.end <- -1
    df$J.start <- -1
    df[[.vend]] = gsub(";", "", df[[.vend]], fixed = T)
    logic = sapply(strsplit(df[[.vend]], "|", T, F, T), length) >= 5
    df$V.end[logic] <- sapply(strsplit(df[[.vend]][logic], "|", T, F, T), "[[", 5)
    logic = sapply(strsplit(df[[.jstart]], "|", T, F, T), length) >= 4
    df$J.start[logic] <- sapply(strsplit(df[[.jstart]][logic], "|", T, F, T), "[[", 4)
  } else {
    df$V.end <- -1
    df$J.start <- -1
    df$Total.insertions <- -1
    df$VD.insertions <- -1
    df$DJ.insertions <- -1
    
    .dj.insertions <- "DJ.insertions"
    .vd.insertions <- "VD.insertions"
  }
  
  .vend <- "V.end"
  .jstart <- "J.start"
  
  if (!is.na(.dalignments)) {
    logic <- sapply(str_split(df[[.dalignments]], "|"), length) >= 5
    df$D5.end <- -1
    df$D3.end <- -1
    df$D5.end[logic] <- sapply(str_split(df[[.dalignments]][logic], "|"), "[[", 4)
    df$D3.end[logic] <- sapply(str_split(df[[.dalignments]][logic], "|"), "[[", 5)
    .dalignments <- c('D5.end', 'D3.end')
  } else {
    df$D5.end <- -1
    df$D3.end <- -1
  }
  
  .dalignments <- c('D5.end', 'D3.end')
  
  if (!(.count %in% table.colnames)) {
    warn_msg = c("  [!] Warning: can't find a column with clonal counts. Setting all clonal counts to 1.")
    warn_msg = c(warn_msg, "\n      Did you apply repLoad to MiXCR file *_alignments.txt?")
    warn_msg = c(warn_msg, " If so please consider moving all *.clonotypes.*.txt MiXCR files to")
    warn_msg = c(warn_msg, " a separate folder and applying repLoad to it.")
    warn_msg = c(warn_msg, "\n      Note: The *_alignments.txt file IS NOT a repertoire file suitable for any analysis.\n")
    cat(warn_msg)
    
    df[[.count]] = 1
  }
  .freq = "Proportion"
  df$Proportion = as.numeric(df[[.count]]) / sum(as.numeric(df[[.count]]))
  
  .aa.seq = IMMCOL$cdr3aa
  df[[.aa.seq]] = bunch_translate(df[[.nuc.seq]])
  
  if (is.na(.big.seq)) {
    .big.seq = "BigSeq"
    df$BigSeq = df[[.nuc.seq]]
  }
  
  df <- df[, make.names(c(.count, .freq,
                          .nuc.seq, .aa.seq,
                          .vgenes, .dgenes, .jgenes,
                          .vend, .dalignments, .jstart,
                          .total.insertions, .vd.insertions, .dj.insertions, .big.seq))]
  
  colnames(df) <- IMMCOL$order
  
  df[[IMMCOL$v]] <- gsub("([*][[:digit:]]*)([(][[:digit:]]*[.]*[[:digit:]]*[)])", "", df[[IMMCOL$v]])
  df[[IMMCOL$v]] <- gsub(",", ", ", df[[IMMCOL$v]])
  df[[IMMCOL$v]] = str_replace_all(df[[IMMCOL$v]], '"', "")
  
  df[[IMMCOL$d]] <- gsub("([*][[:digit:]]*)([(][[:digit:]]*[.]*[[:digit:]]*[)])", "", df[[IMMCOL$d]])
  df[[IMMCOL$d]] <- gsub(",", ", ", df[[IMMCOL$d]])
  df[[IMMCOL$d]] = str_replace_all(df[[IMMCOL$d]], '"', "")
  
  df[[IMMCOL$j]] <- gsub("([*][[:digit:]]*)([(][[:digit:]]*[.]*[[:digit:]]*[)])", "", df[[IMMCOL$j]])
  df[[IMMCOL$j]] <- gsub(",", ", ", df[[IMMCOL$j]])
  df[[IMMCOL$j]] = str_replace_all(df[[IMMCOL$j]], '"', "").immunarch::postprocess(fix.alleles(df))
}



