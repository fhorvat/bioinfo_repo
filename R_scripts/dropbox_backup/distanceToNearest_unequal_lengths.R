# finding distance to other elements of same type
distance_to_same <- 
  distanceToNearest(all_elements_original_ranges, ignore.strand = T) %>% 
  as.data.frame() %>% 
  right_join(., data.frame(queryHits = 1:length(all_elements_original_ranges)), by = "queryHits")

# adding distance to other elements to orignal ranges 
all_elements_original_ranges$distance_to_same <- distance_to_same$distance
