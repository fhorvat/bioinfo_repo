element_fpkm_filtered <- 
  element_genomic_filtered_fpkm %>% 
  mutate_at(.cols = vars(starts_with("s_")), 
            .funs = funs(ifelse(. > fpkm_limit, NA, .)))
# The use of funs() here indicates that ifelse(. > fpkm_limit, NA, .) is an anonymous function 
# that is being defined within the call to mutate_at().
