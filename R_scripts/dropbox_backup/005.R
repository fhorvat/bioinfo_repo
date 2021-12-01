My.fib <- function(x, n){
  if(x == 1){
    fib <- 1
  }
  else{
    fib <- c(1, 1)
    for (i in 1:x){
      fib[i + 2] <- fib[i] + fib[i + 1]
    }    
  }
  fib <- fib[fib <= x]
  pos.x <- length(fib)
  pos.n <- pos.x - n + 1
  if (n <= pos.x){
    fib.n <- fib[pos.n:pos.x]
  } 
  else{
    fib.n <- rep_len(fib[1:pos.x], length.out = n)
    if (isTRUE(all.equal(n/pos.x, as.integer(n/pos.x))) == FALSE){
      warning("recikliranje nije potpuno")
    } 
  }
return(fib.n)
}
My.fib(13, 4)

Fibonacciate <- function(v, ncols){
  v <- as.matrix(v)
  mat1 <- t(apply(v, 1, My.fib, n = ncols))
  return(mat1)
}

v1 <- 13
Fibonacciate(v1, 4)

