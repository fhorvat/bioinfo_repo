# matrice
m1 = matrix(1:6, nrow=3, ncol=2)
m2 = matrix(1:6, nrow=3, ncol=2, byrow=TRUE)
m1
m2
m1 + 10

m1 + m2
m1 + 10:12
m1 + 10:13
m1 + 10:11
v = 1:5
v[2:3]

m3 = matrix(1:30, nrow=5)
m3[3]
m3[11]
m3[1,]
m3[,3]
m3[,3,drop=FALSE]
m3[3,,drop=FALSE]
dim(m3)
a = m3[,3]
dim(a)
m3[,2:3]
m3[2:3,2:3]

m3[m3[,3] < 14,]
rowSums(m3)

a = 1:5
b = 10:14
matrix(c(a,b), ncol=2)
matrix(c(a,b))
cbind(a,b)

# liste
l1 = list(1:3, 5:8)
l2 = list(1:3, c(T,T,F))
l1
l2
class(l2[2])
l2[2]
l2[[2]]

l2[2] + 1
l2[[2]] + 1
   
l2
l3 = c(l2[2],1:5)  
l3
c(l2[2],list(1:5))
length(1:5)
length(list(1:5))

?sum

sum(1,2,3,10,11)
sum(c(1,2,3,10,11))

# character
a = '1'
a + 1
c(a,a,a,a,a)
c(1:5,a)
c(T,1,'a')
c('a')

a = '1'
c(T,1,a)
is.character(a)
as.numeric(a)
as.numeric(c(T,1,a))
as.numeric('TRUE')

# factor
a = c(1,1,1,2,2,3)
factor(a)
a = factor(c('3','3','3','2','2','1'), levels=3:1)
as.numeric(a)
a = factor(rep(c('M','F'), times=c(2,3)))
as.numeric(a)
a = factor(c(2,2,2,5,5,7))
as.numeric(a)

#data.frame
d = data.frame(a = 1:5, b = rep(T, 5), l = c('a','b','c','d','e'))
d[1:3,]
d[,2]
colnames(d)
d$l
d$l[d$a < 3]

print('nesto',1)
print(c('nesto',1))
y = 1
print(paste('nesto',y))

'%p%' = function(x,y){
  x+y
}
1 %p% 2

# control stuctures












