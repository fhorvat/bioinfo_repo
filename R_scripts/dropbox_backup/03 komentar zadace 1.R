####1. KREATIVNIJE ZBRAJANJE####
suma1<-function(x,y){z<-seq(x,y) 
                     sum(z)}



suma2<-function(x,y){z<- x:y
                     sum(z)}

###

f1 = function(m,n){
  a = sum(m:n)
  return(a)
}

f2 = function(m,n){
  z = n-m
  w = m + z*m + sum(1:z)
#   if (n > m)
#     return(w)
#   else 
#     z1 = z*-1
#   w = n + z1*n + sum(1:z1)
#   return (w)
  
}

#### 2. DOMENA FUNKCIJE ####
# indentacija

f1 <- function(n){
  if (n == 1)
    return (1)
  else return (n * factorial(n-1))
}

#### 3. CODE STYLE ####
for (i in 1:length(Xn)) red(1,Xn[i])->yn1[i]
for (i in 1:length(Xn)) red(3,Xn[i])->yn2[i]
for (i in 1:length(Xn)) red(6,Xn[i])->yn3[i]

###
manh <- function(v1,v2){
  return(sum(abs(v1-v2)))
}

eucl <- function (v1,v2){
  return(sqrt(sum((v2-v1)^2)))
}

cheb <- function (v1,v2){
  return(max(abs(v2-v1)))
}

dist <- function (v1,v2,metric){
  if (metric == 'manhattan'){manh (v1,v2)} else 
    if (metric == 'euclidean') {eucl (v1,v2)} else 
      if (metric == 'chebyshev') {cheb (v1,v2)} else 
        print ('Error: distance metric must be manhattan, euclidean or chebyshev')
}

#### 4. TABLE VAĐENJE VRIJEDNOSTI ####
disc <- c(5, 3, 0, 2, 0, 3, 2, 3, 6, 1, 2, 1, 2, 1, 3, 3, 3, 5, 2, 4, 4, 0, 2, 3, 7, 12, 3, 10, 9, 2, 3, 7,
          7, 2, 3, 3, 6, 2, 4, 3, 5, 2, 2, 4, 0, 4, 2, 5, 2, 3, 3, 6, 5, 8, 3, 6, 6, 0, 5, 2, 2, 2, 6, 3, 4, 4, 2, 2, 4,
          7, 5, 3, 3, 0, 2, 2, 2, 1, 3, 4, 2, 2, 1, 1, 1, 2, 1, 4, 4, 3, 2, 1, 4, 1, 1, 1, 0, 0, 2, 0);

tbl <- table(disc)
sum(tbl*as.numeric(names(tbl)))/sum(tbl)

### 5. UKLANJANJE VEKTORA S NULOM ###
v1 <- c(0, 1, 2, 3, 4)
v2 <- c(0, 2, 3, 5, 7)
new.v1 <- v1[v1!= 0]
new.v2 <- v2[v2!= 0]
new.v1
new.v2


v3 <- c(1, 2, 3, 4, 5)
v4 <- c(0, 2, 3, 5, 7)
new.v3 <- v3[v3!= 0]
new.v4 <- v4[v4!= 0]
new.v3
new.v4


mask <- (v3 != 0) & (v4 != 0)
mask
my.new.v3 <- v3[mask]
my.new.v4 <- v4[mask]
my.new.v3
my.new.v4


v3
v4
v3 & v4

v5 <- c(1, 2, 3, 4, 5)
v6 <- c(0, -2, 3, 5, 7)
v5&v6

v7 <- c(1, 2.5, 3, 4, 5)
v8 <- c(0.1, 2, 3.7, 5, 7)
v7 & v8

