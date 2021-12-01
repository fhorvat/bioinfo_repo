  # uzorak 1
s1 <- rpois (100000, lambda = 6.5) 
  # mean je jednak lambda
sd1 <- sd(s1)
sd1 
  # standardna devijacija uzorka 1 (mean = 6.5)
var1 <- var(s1)
var1 
  # varijanca uzorka 1 (mean = 6.5)
dif1 <- abs(var1-6.5)
dif1
  # vrijednost varijance je vrlo bliska srednjoj vrijednosti
  # uzorka 
  
  # uzorak 2
s2 <- rpois (100000, lambda = 4.7) 
var2 <- var(s2)
dif2 <- abs(var2-4.7)
dif2


  # uzorak 3
s3 <- rpois (100000, lambda = 3.2) 
var3 <- var(s3)
dif3 <- abs(var3-3.2)
dif3

  # uzorak 4
s4 <- rpois (100000, lambda = 4.11) 
var4 <- var(s4)
dif4 <- abs(var4-4.11)
dif4


  # uzorak 5
s5 <- rpois (100000, lambda = 1.23) 
var5 <- var(s5)
dif5 <- abs(var5-1.23)
dif5

  # vrijednosti varijance svih uzoraka su vrlo bliske srednjoj 
  # vrijednosti uzorka, razlika te dvije vrijednosti je vrlo mala


