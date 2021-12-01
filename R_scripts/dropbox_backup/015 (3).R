fun0 <- function (fun1,fun2,int){
  x <- seq(int[1],int[2],length(int)/199)
  
  fun1 <- function (x){  
    y1 <- factorial(x)
    
    fun2 <- function (x){
      y2 <- (sqrt(2*pi*x)*((x/exp(1))^x))
      matplot(x, cbind(y1,y2),type="l",col=c("red","green"),lty=c(1,1))
      dif <- mean((y1+y2)/y1)
      return(dif)
    }
    fun2(x)
  }
  fun1(x)
}

fun0(fun1,fun2,c(2,4))

#na intervalu [2,4] funkcije izgledaju gotovo potpuno jednako
#Stirlingova aproksimacija postaje sve manje precizna kako se x poveæava
#poveæanjem intervala raste i prosjeèna razlika izmeðu funkcija
