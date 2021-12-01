easyGregexpr <- function(pattern, charvector, ...) {
  require(plyr)
  
  if (storage.mode(charvector) != "character") stop("easyGregexpr expects charvector to be a character vector.")
  
  #identify all matches
  regexpMatches <- gregexpr(pattern, charvector, ...)
  
  convertMatches <- c()
  for (i in 1:length(regexpMatches)) {
    thisLine <- regexpMatches[[i]]
    #only append if there is at least one match on this line
    if (thisLine[1] != -1) {
      convertMatches <- rbind(convertMatches, data.frame(element=i, start=thisLine, end=thisLine + attr(thisLine, "match.length") - 1))
    }
  }
  
  #if no matches exist, return null (otherwise, will break adply)
  if (is.null(convertMatches)) return(NULL)
  
  #We now have a data frame with the line, starting position, and ending position of every match
  #Add the matched string to the data.frame
  #Use adply to iterate over rows and apply substr func
  convertMatches <- adply(convertMatches, 1, function(row) {
    row$match <- substr(charvector[row$element], row$start, row$end)
    return(as.data.frame(row))
  })
  
  #need to convert from factor to character because adply uses stringsAsFactors=TRUE even when overridden
  convertMatches$match <- as.character(convertMatches$match)
  return(convertMatches)
}