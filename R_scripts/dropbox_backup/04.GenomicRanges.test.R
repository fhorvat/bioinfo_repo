library(GenomicRanges)
library(IRanges)
library(ggbio)
library(ggplot2)
library(tibble)
library(dplyr)
library(magrittr)

IRangesMethodToTable <- function(IRangesMethod, ...){
  
  # apply method to IRanges
  irange_result <- IRangesMethod(...)
  
  # create list of IRanges
  irange_list <- list(..., irange_result)
  
  # loop through the list: transform to table, add unique id, bind to one table
  irange_tb <- lapply(1:length(irange_list), function(index){
    
    if(class(irange_list[[index]]) == "IRanges"){
      
      # get disjoint bins
      disjoin_bins <- irange_list[[index]] %>% disjointBins()
      
      irange_list[[index]] %>% 
        as_tibble(.) %>% 
        dplyr::mutate(range_id = index, 
                      range_ymin = disjoin_bins)
    }
    
  }) %>% 
    dplyr::bind_rows(.) %>% 
    
    # replace last added range with NA
    dplyr::mutate(range_id = ifelse(test = (range_id == max(range_id)), 
                                    yes = NA, 
                                    no = paste("range", range_id, sep = ".")))
  
  # return table
  return(irange_tb)
  
}


plotIRanges <- function(irange_tb){
  
  # calculate max. x-limit
  xlim_max <- max(irange_tb$end)
  
  # add range_ymin to separate ranges on y-axis
  irange_tb %<>% 
    dplyr::mutate(range_ymin = 1:n()) 
  
  # add range ID if it doesn't exist
  if(!("range_id" %in% names(irange_tb))){
    irange_tb %<>% 
      dplyr::mutate(range_id = "range.1")
  }
    
  # plot
  ggplot(irange_tb, aes(xmin = start - 0.4, xmax = end + 0.4, 
                        ymin = range_ymin - 1, ymax = range_ymin - 0.1, 
                        fill = range_id)) + 
    geom_rect() +
    scale_x_continuous(breaks = 1:xlim_max) +
    scale_y_reverse() +
    scale_fill_grey(start = 0.8, end = 0.3, na.value = "red3") + 
    xlim(1, xlim_max) +
    theme_bw() + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_blank(), 
          panel.grid.minor.y = element_blank(), 
          legend.position = "none")
  
}




### set - normal
# (p)union() # union
# (p)intersect() # intersect
# (p)setdiff() # difference
# (p)gaps() # complement
range_1 <- IRanges(start = c(4, 18, 1), end = c(7, 22, 4))
range_2 <- IRanges(start = c(8, 3, 16), end = c(13, 13, 20))
IRangesMethodToTable(union, range_1, range_2) %>% 
  plotIRanges(.)
IRangesMethodToTable(punion, fill.gap = TRUE, range_1, range_2) %>% 
  plotIRanges(.)

range_1 <- IRanges(start = c(9, 18, 1), end = c(13, 22, 9))
range_2 <- IRanges(start = c(8, 3, 2), end = c(13, 20, 5))
IRangesMethodToTable(intersect, range_1, range_2) %>% 
  plotIRanges(.)
IRangesMethodToTable(pintersect, resolve.empty = "max.start", range_1, range_2) %>% 
  plotIRanges(.)

range_1 <- IRanges(start = c(6, 18, 1), end = c(13, 22, 4))
range_2 <- IRanges(start = c(12, 9, 2), end = c(19, 20, 5))
IRangesMethodToTable(setdiff, range_1, range_2) %>% plotIRanges(.)
IRangesMethodToTable(psetdiff, range_1, range_2) %>% plotIRanges(.)

# - **gaps()** takes one IRanges object as argument 
# - it is one of the inter range transformations
range_1 <- IRanges(start = c(4, 18, 1), end = c(7, 22, 4))
range_2 <- IRanges(start = c(10, 3, 16), end = c(13, 13, 20))
IRangesMethodToTable(gaps, c(range_1, range_2)) %>% plotIRanges(.)
IRangesMethodToTable(pgap, range_1, range_2) %>% plotIRanges(.)


### intra
# shift()
# narrow()
# resize()
# flank()
# promoters()
# reflect()
# restrict()
range_1 <- IRanges(start = c(4, 18, 5, 8, 3, 16), 
                   end = c(9, 22, 9, 13, 13, 20))

IRangesMethodToTable(shift, shift = 7, range_1) %>% 
  plotIRanges(.)
IRangesMethodToTable(shift, shift = -2, range_1) %>% 
  plotIRanges(.)

IRangesMethodToTable(narrow, start = 2, end = 5, range_1) %>% 
  plotIRanges(.)

IRangesMethodToTable(resize, width = 3, fix = "start", range_1) %>% 
  plotIRanges(.)
IRangesMethodToTable(resize, width = 3, fix = "end", range_1) %>% 
  plotIRanges(.)
IRangesMethodToTable(resize, width = 3, fix = "center", range_1) %>% 
  plotIRanges(.)

IRangesMethodToTable(flank, range_1, width = 2, start = TRUE) %>%  plotIRanges(.)
IRangesMethodToTable(flank, range_1, width = 2, start = FALSE) %>%plotIRanges(.)
IRangesMethodToTable(flank, range_1, width = 2, start = TRUE, both = TRUE) %>% plotIRanges(.)

# similar to **flank**, but takes start as reference point
IRangesMethodToTable(promoters, upstream = 2, downstream = 2, range_1) %>% plotIRanges(.)

IRangesMethodToTable(reflect, range_1, bounds = IRanges(start = 14, end = 16)) %>% plotIRanges(.)

IRangesMethodToTable(restrict, range_1, start = 3, end = 9) %>% plotIRanges(.)


### inter
# range()
# reduce()
# disjoin()
IRangesMethodToTable(range, range_1) %>% 
  plotIRanges(.)

IRangesMethodToTable(reduce, range_1) %>% 
  plotIRanges(.)

IRangesMethodToTable(disjoin, range_1) %>% 
  plotIRanges(.)

disjointBins(range_2)


### other 
c(range_1, range_2)
range_2 + 1
range_2 - 1


### GenomicRanges
# plot
grange_1 <- GRanges("chr1", range_1, strand = "+")
plotIRanges(as_tibble(grange_1))

ggbio(autoplot(c(GRanges("chr1", range_1, strand = "+"), 
                 GRanges("chr2", range_2, strand = "-"))))

(grange_1 <- GRanges("chr1", range_1, strand = "+"))
(grange_2 <- GRanges("chr1", range_2, strand = "-"))

seqnames(grange_1)
strand(grange_1)
seqlevels(grange_1)
genome(grange_1)




