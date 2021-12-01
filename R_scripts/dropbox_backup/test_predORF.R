function (x, n = 1, type = "grl", mode = "orf", strand = "sense", 
          longest_disjoint = FALSE, startcodon = "ATG", stopcodon = c("TAA", 
                                                                      "TAG", "TGA")) 
{
  
  x <- geni_skupa["ENSEEUG00000011205"]
  n = 5
  type = "df"
  mode = "orf"
  strand = "both"
  longest_disjoint=FALSE
  startcodon = "ATG"
  stopcodon = c("TAA", "TAG", "TGA")
  
  if (any(nchar(c(startcodon, stopcodon)) != 3)) 
    stop("startcodon and stopcodons can only contain 3-letter strings.")
  if (!toupper(mode) %in% c("ORF", "CDS")) 
    stop("'mode' can only be assigned one of: 'orf' or 'cds'")
  if (length(names(x)) == 0 | any(duplicated(names(x)))) 
    stop("Sequence name slot of x need be populated with unique names.")
  x <- x[width(x) >= 6]
  
  
  .predORF <- function(x, n, mode, strand, ...) {
    
    
    tmp_pos <- lapply(x, function(y) .predORF(y, n = n, 
                                              mode = mode, strand = "sense", longest_disjoint, 
                                              startcodon, stopcodon))
    
    x <- x[[1]]
    strand = "sense"
    
    if (tolower(strand) == "sense") {
      mystrand <- 1
      startcodon <- startcodon
      stopcodon <- stopcodon
    }
    else if (tolower(strand) == "antisense") {
      mystrand <- 2
      stopcodon_temp <- as.character(reverseComplement(DNAStringSet(stopcodon)))
      startcodon_temp <- as.character(reverseComplement(DNAStringSet(startcodon)))
      stopcodon <- startcodon_temp
      startcodon <- stopcodon_temp
    }
    else {
      stop("strand can only be assigned 'sense', 'antisense' or 'both'")
    }
    if (alphabetFrequency(x)["N"] > 0) {
      orfRanges <- cbind(subject_id = numeric(), start = numeric(), 
                         end = numeric(), width = numeric(), strand = numeric(), 
                         inframe2end = numeric())
      warning("Skipped sequence containing Ns.")
      return(orfRanges)
    }
    c1 <- as.character(suppressWarnings(codons(x)))
    c2 <- as.character(suppressWarnings(codons(x[2:length(x)])))
    c3 <- as.character(suppressWarnings(codons(x[3:length(x)])))
    startpos1 <- which(c1 %in% startcodon)
    stoppos1 <- which(c1 %in% stopcodon)
    startpos2 <- which(c2 %in% startcodon)
    stoppos2 <- which(c2 %in% stopcodon)
    startpos3 <- which(c3 %in% startcodon)
    stoppos3 <- which(c3 %in% stopcodon)
    if (mode == "cds") {
      stoppos1 <- unique(c(0, stoppos1, length(c1) + 1))
      stoppos2 <- unique(c(0, stoppos2, length(c2) + 1))
      stoppos3 <- unique(c(0, stoppos3, length(c3) + 1))
      startpos1 <- stoppos1
      startpos2 <- stoppos2
      startpos3 <- stoppos3
    }
    if (tolower(strand) == "sense") {
      orfpos1 <- t(sapply(seq(along = startpos1), function(x) c((startpos1[x] * 
                                                                   3) - 2, stoppos1[stoppos1 > startpos1[x]][1] * 
                                                                  3)))
      if (length(orfpos1) == 0) 
        orfpos1 <- matrix(nrow = 0, ncol = 2)
      orfpos2 <- t(sapply(seq(along = startpos2), function(x) c((startpos2[x] * 
                                                                   3) - 1, (stoppos2[stoppos2 > startpos2[x]][1] * 
                                                                              3) + 1)))
      if (length(orfpos2) == 0) 
        orfpos2 <- matrix(nrow = 0, ncol = 2)
      orfpos3 <- t(sapply(seq(along = startpos3), function(x) c((startpos3[x] * 
                                                                   3) - 0, (stoppos3[stoppos3 > startpos3[x]][1] * 
                                                                              3) + 2)))
      if (length(orfpos3) == 0) 
        orfpos3 <- matrix(nrow = 0, ncol = 2)
    }
    if (tolower(strand) == "antisense") {
      startpos1 <- sort(startpos1, decreasing = TRUE)
      startpos2 <- sort(startpos2, decreasing = TRUE)
      startpos3 <- sort(startpos3, decreasing = TRUE)
      stoppos1 <- sort(stoppos1, decreasing = TRUE)
      stoppos2 <- sort(stoppos2, decreasing = TRUE)
      stoppos3 <- sort(stoppos3, decreasing = TRUE)
      orfpos1 <- t(sapply(seq(along = stoppos1), function(x) c((startpos1[startpos1 < 
                                                                            stoppos1[x]][1] * 3) - 2, stoppos1[x] * 3)))
      if (length(orfpos1) == 0) 
        orfpos1 <- matrix(nrow = 0, ncol = 2)
      orfpos2 <- t(sapply(seq(along = stoppos2), function(x) c((startpos2[startpos2 < 
                                                                            stoppos2[x]][1] * 3) - 1, (stoppos2[x] * 3) + 
                                                                 1)))
      if (length(orfpos2) == 0) 
        orfpos2 <- matrix(nrow = 0, ncol = 2)
      orfpos3 <- t(sapply(seq(along = stoppos3), function(x) c((startpos3[startpos3 < 
                                                                            stoppos3[x]][1] * 3) - 0, (stoppos3[x] * 3) + 
                                                                 2)))
      if (length(orfpos3) == 0) 
        orfpos3 <- matrix(nrow = 0, ncol = 2)
    }
    
    orfRanges <- rbind(orfpos1, orfpos2, orfpos3)
    orfRanges <- na.omit(orfRanges)
    if (mode == "cds") {
      orfRanges[, 1] <- orfRanges[, 1] + 3
      orfRanges[orfRanges[, 2] > length(x), 2] <- length(x)
    }
    orfRanges <- IRanges(start = orfRanges[, 1], end = orfRanges[, 
                                                                 2])
    orfRanges <- orfRanges[rev(order(width(orfRanges)))]
    orfRanges <- as.data.frame(orfRanges)
    inframe <- (length(x) - orfRanges$end)/3
    inframe2 <- inframe
    inframe2[abs((inframe - as.integer(inframe))) == 0] <- 1
    inframe2[abs(round((inframe - as.integer(inframe)), 2)) == 
               round(1/3, 2)] <- 2
    inframe2[abs(round((inframe - as.integer(inframe)), 2)) == 
               round(2/3, 2)] <- 3
    if (nrow(orfRanges) > 0) {
      orfRanges <- cbind(subject_id = 1:length(orfRanges[, 
                                                         1]), orfRanges, strand = mystrand, inframe2end = inframe2)
    }
    else {
      orfRanges <- cbind(subject_id = numeric(), start = numeric(), 
                         end = numeric(), width = numeric(), strand = numeric(), 
                         inframe2end = numeric())
    }
    if (n == "all") {
      if (nrow(orfRanges) > 0 & longest_disjoint == TRUE) {
        tmpgr <- GRanges(seqnames = "dummy", IRanges(orfRanges[, 
                                                               2], orfRanges[, 3]), strand = "+")
        orfRanges <- orfRanges[disjointBins(tmpgr) == 
                                 1, , drop = FALSE]
        orfRanges[, 1] <- 1:nrow(orfRanges)
      }
      return(orfRanges)
    }
    else if (is.numeric(n)) {
      if (nrow(orfRanges) < n & nrow(orfRanges) != 0) {
        upperlimit <- length(orfRanges[, 1])
      }
      else if (nrow(orfRanges) > 0) {
        upperlimit <- n
      }
      else {
        upperlimit <- NULL
      }
      if (is.null(upperlimit)) {
        return(orfRanges[upperlimit, , drop = FALSE])
      }
      else {
        return(orfRanges[1:upperlimit, , drop = FALSE])
      }
    }
    else {
      stop("n needs to be assigned positive integer or 'all'")
    }
  }
  
  if (class(x) == "DNAString") {
    if (tolower(strand) == "sense" | tolower(strand) == 
        "antisense") {
      tmp <- .predORF(x, n, mode, strand, longest_disjoint, 
                      startcodon, stopcodon)
      tmp[, "strand"] <- ifelse(as.numeric(tmp[, 
                                               "strand"]) == 1, "+", "-")
    }
    else if (tolower(strand) == "both") {
      tmp_pos <- .predORF(x, n, mode, strand = "sense", 
                          longest_disjoint, startcodon, stopcodon)
      tmp_neg <- .predORF(x, n, mode, strand = "antisense", 
                          longest_disjoint, startcodon, stopcodon)
      tmp <- rbind(tmp_pos, tmp_neg)
      tmp[, "strand"] <- ifelse(as.numeric(tmp[, 
                                               "strand"]) == 1, "+", "-")
    }
    else {
      stop("strand can only be assigned 'sense', 'antisense' or 'both'")
    }
    if (type == "df") 
      return(tmp)
    if (type == "gr") 
      return(makeGRangesFromDataFrame(data.frame(seqnames = "unknown", 
                                                 tmp), keep.extra.columns = TRUE))
  }
  
  if (class(x) == "DNAStringSet") {
    if (tolower(strand) == "sense" | tolower(strand) == "antisense") {
      tmp <- lapply(x, function(y) .predORF(y, n = n, mode = mode, 
                                            strand, longest_disjoint, startcodon, stopcodon))
      names(tmp) <- names(x)
      tmpdf <- do.call("rbind", tmp)
      rownames(tmpdf) <- NULL
      tmpdf <- data.frame(seqnames = rep(names(tmp), sapply(tmp, 
                                                            nrow)), tmpdf)
      tmpdf[, "strand"] <- ifelse(as.numeric(tmpdf[, 
                                                   "strand"]) == 1, "+", "-")
    }
    else if (tolower(strand) == "both") {
      
      tmp_pos <- lapply(x, function(y) .predORF(y, n = n, 
                                                mode = mode, strand = "sense", longest_disjoint, 
                                                startcodon, stopcodon))
      
      names(tmp_pos) <- names(x)
      tmpdf_pos <- do.call("rbind", tmp_pos)
      rownames(tmpdf_pos) <- NULL
      tmpdf_pos <- data.frame(seqnames = rep(names(tmp_pos), sapply(tmp_pos, nrow)), tmpdf_pos)
      
      
      tmp_neg <- lapply(x, function(y) .predORF(y, n = n, 
                                                mode = mode, strand = "antisense", longest_disjoint, 
                                                startcodon, stopcodon))
      names(tmp_neg) <- names(x)
      tmpdf_neg <- do.call("rbind", tmp_neg)
      rownames(tmpdf_neg) <- NULL
      tmpdf_neg <- data.frame(seqnames = rep(names(tmp_neg), 
                                             sapply(tmp_neg, nrow)), tmpdf_neg)
      tmpdf <- rbind(tmpdf_pos, tmpdf_neg)
      tmpdf[, "strand"] <- ifelse(as.numeric(tmpdf[, 
                                                   "strand"]) == 1, "+", "-")
      tmpdf <- tmpdf[order(tmpdf$seqnames, tmpdf$subject_id), 
                     ]
    }
    else {
      stop("strand can only be assigned 'sense', 'antisense' or 'both'")
    }
    if (tolower(type) == "df") {
      rownames(tmpdf) <- NULL
      return(tmpdf)
    }
    else if (tolower(type) == "gr") {
      tmpdf <- makeGRangesFromDataFrame(tmpdf, keep.extra.columns = TRUE)
      return(tmpdf)
    }
    else if (tolower(type) == "grl") {
      tmpdf <- makeGRangesFromDataFrame(tmpdf, keep.extra.columns = TRUE)
      tmpdf <- split(tmpdf, as.character(seqnames(tmpdf)))
      return(tmpdf)
    }
    else {
      stop("type can only be assigned 'df' of 'gr'")
    }
  }
}