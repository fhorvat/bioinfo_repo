df2 <- data.frame(v1 = 1:3, v2 = 5:7, v3 = c("x", "y", "z"))
df1 <- data.frame(v1 = (1:3)*10, v2 = (5:7)*10,
                  row.names = c("xxx", "yyy", "zzz"))

razlika.col <- setdiff(colnames(df1), colnames(df2))
if (length(razlika.col) > 0){
  imena.col <- colnames(df2)
  for (i in 1:length(razlika.col)){
    df2 <- cbind(df2, NA)
  }
  colnames(df2) <- c(imena.col, razlika.col)
}

razlika.col <- setdiff(colnames(df2), colnames(df1))
if (length(razlika.col) > 0){
  imena.col <- colnames(df1)
  for (i in 1:length(razlika.col)){
    df1 <- cbind(df1, NA)
  }
  colnames(df1) <- c(imena.col, razlika.col)
}

