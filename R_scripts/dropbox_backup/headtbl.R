# converts to data.frame and diplays head (usefull for tibble and data.table)
headtbl <- function(df){
    head(as.data.frame(df))
}