'%my.rbind%' <- function(df1, df2){
  if (is.data.frame(df1) & is.data.frame(df2) == TRUE){
    
# sluèaj 1: df1 i df2 jednake dužine, imena stupaca = df1
    
    razlika.col <- ncol(df1) - ncol(df2)
    if (razlika.col == 0){ 
      # provjera razlikuju li se imena stupaca df1 i df2:
      df.col.diff <- setdiff(colnames(df1), colnames(df2))
      if (length(df.col.diff) != 0){
        warning("Imena kolona koje se preklapaju nisu jednaka")
      }   
      colnames(df2) <- colnames(df1)
    }
  
    
# sluèaj 2: df1 duži od df2, imena stupaca = df1
    
    razlika.col <- ncol(df1) - ncol(df2)
    if (razlika.col > 0){
      # provjera razlikuju li se imena stupaca df1 i df2:
      df.col.diff <- setdiff(colnames(df1[1:ncol(df2)]), 
                             colnames(df2))
      if (length(df.col.diff) != 0){
        warning("Imena kolona koje se preklapaju nisu jednaka")
      } 
      imena.col <- colnames(df1)
      for (i in 1 : razlika.col){
        df2 <- cbind(df2, NA)
      }
      colnames(df2) <- imena.col
    }

# sluèaj 3: df2 duži od df1, imena stupaca = poèetak df1, ostatak df2 
    
    razlika.col <- ncol(df2) - ncol(df1)
    if (razlika.col > 0){
      # provjera razlikuju li se imena stupaca df1 i df2:
      df.col.diff <- setdiff(colnames(df1), 
                             colnames(df2[1:ncol(df1)]))
      if (length(df.col.diff) != 0){
        warning("Imena kolona koje se preklapaju nisu jednaka")
      }
      imena.col <- colnames(df1)
      imena.col2 <- colnames(df2[(ncol(df1) + 1) : ncol(df2)])
      for (i in 1 : razlika.col){
        df1 <- cbind(df1, NA)
      }
      colnames(df2) <- c(imena.col, imena.col2)
      colnames(df1) <- c(imena.col, imena.col2)
    }
    return(rbind(df1, df2))
  }
  
  else{
    stop("Error: %my.rbind% radi samo s data.frame argumentima")
  }
}

my.df1 <- data.frame(a1 = 1:3, a2 = 5:7, a3 = c("a", "b", "c"),
                     row.names = 1:3) 
my.df2 <- data.frame(a1 = (1:3)*10, a2 = (5:7)*10, 
                     a3 = c("xx", "yy", "zz"), a4 = c(2:4),
                     row.names = c("xxx", "yyy", "zzz"))

my.df1 %my.rbind% my.df2
