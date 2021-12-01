m1 <- matrix (1:9, nrow = 8, ncol = 3)

My.replacement <- structure(NA, class = "My.replacement")

"(<-.My.replacement" <- function(My.replacement, m1, value) {
  b <- ceiling(nrow(m1)/2)
  if (nrow(m1) %% 2 == 0){
    m1[c(b,b+1), ] <- value
  }
  else {
    m1[b, ] <- value
  }
 My.replacement
 m1 <- My.replacement
 return(m1)
}

My.replacement(m1)<- 100:102
m1
