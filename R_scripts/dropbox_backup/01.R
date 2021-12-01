1+1

a = 1
a 

a + 2

1:3
a = 1:3
a + 1 

b = 2
a + b

b = 2:3
a = 1:4

a + b

b = 2:3
a = 1:3
d = a + b
= 1
d + 5

a = 1:3

a > 2
a

a = 2
a = 1:3
i = a > 2
i
i + 1

a = 4:6
a[2]
a[2:1]

i[3]

a = 1:3
i = a > 2
a = 4:6
a[i]

a = 1:10
b = a[ a >= 7 ]
b
a == 3

a = 10:1
b = 1:10
b[ a == 3 ]

b[ a = 3 ]
a
a = 3
b[a]

b = 1:10
b
b[-3]
b = b[-3]

a = c(1,3,7)
d = c(a,a)
d
c(a,a)


c(a,1:10)

a = c(1,3,7)
c(a, c(FALSE,TRUE,TRUE))
c(a, c(F,T,T))
c(a, a(F,T,T))
c(a, a(T,F,T))

c(a, a[c(F,T,T)])
c(a, c(3,7))
c(a, a[c(T,F,T)])

d = NULL
d + 1
d = NA
d + 1

sum(1:3)
1 + 2 + 3
?sum
?exp
exp(2,3)
3^2

a <- 3
a
4 -> a
a
# <<-
#a 
# isuse kako je ovo dosadno

# 1.
f1 = function(x,y){
  
  z = x + y
  
  #return(z)
  z
}

a = f1(3,4)

# 2.
f2 = function(x,y){
  z = x + y
}

b = NULL
b = f2(3,4)

# 3.
f3 = function(x,y){
   x + y
}
f3(3,5)

#4.
f4 = function(x=1,y=6){
  x + y
}
f3()
f4()

#f5.
a = 10
f5 = function(x, a = 2){
  x + a
}
f5(3)

# 6. 
a = 2
f6 = function(x){
  a = 5
  x + a
}

f6(3)
a

a = 2
a = a + 1

a = 2; a = a + 1

g2 = function(x,y){
  c(x,y)
}

g1 = function(x, y){ 
  g2(x,y) + 2
}
g1(1,3)
g2(1,3)

p = function(x){
  print(x+1)
}
p(3)

a = 1
a = p(3)
a 

































