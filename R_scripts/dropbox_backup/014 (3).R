fun1 <- function (x){  
  y1 <- sin(x)
  return(y1)
}
  
fun2 <- function (x){
  y2 <- x
  return(y2)
}

fun0 <- function (fun1,fun2,int){
  x <- seq(int[1],int[2],length.out=100)
  plot (x, fun1(x), type='l', col='green')
  lines (x, fun2(x), col="red")
}
fun0(fun1,fun2,c(0,1))  

int1 <- c(0,1)
a <- 0
a <- seq(int1[1],int1[2],length.out=100)
b <- fun1(a)
d <- fun2(a)

xyz <- (b != 0) & (d != 0)
b <- b[xyz]
d <- d[xyz]
diff <- abs(mean((b-d)/b)*100)
diff

