# author: vfranke
### Merges a list of data tables given their ids

library(data.table)

MergeDataTable = function(l, key = NULL, all = TRUE, fill = NULL){
  
  if(is.null(key))
    stop('You have to give a key')
  
  l = lapply(l, function(x) setkeyv(x, cols = key))
  r = Reduce(function(x, y) merge(x, y, all = all, by = key), l)
  if(!is.null(fill))
    r[is.na(r)] = fill
  return(r)
}

