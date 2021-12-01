### INFO: converts to data.frame and diplays head (usefull for tibble and data.table)
### DATE: 28. 08. 2017.
### AUTHOR: Filip Horvat

######################################################## FUNCTIONS
headt <- function(df){
    head(as.data.frame(df))
}
