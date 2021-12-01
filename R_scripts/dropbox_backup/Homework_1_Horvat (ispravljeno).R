#############
# 1. excerise
#############
zbr1 <- function(m,n){
  sum(m:n)
}
zbr1(1,5)

zbr2 <- function (m,n){
  (m+n)*(n-m+1)/2
}
zbr2(1,6)

#############
# 2. excerise
#############
fact1 <- function (n) {
  if (n==0) return (1) else 
    return (n*fact1(n-1))
}
fact1(6)

fact2 <-function (n) {
  return (sqrt(2*pi*n)*((n/exp(1))^n))
}
fact2(4)

#############
# 3. excerise
#############
manh <- function(v1,v2){
  sum(abs(v1-v2))
}
v1 <- c(1,6,19)
v2 <- c(3,5,8)
manh (v1,v2)

#############
# 4. excerise
#############
eucl <- function (v1,v2){
  sqrt(sum((v2-v1)^2))
}
v1 <- c(4,9,18)
v2 <- c(5,7,4)
eucl(v1,v2)

#############
# 5. excerise
#############
cheb <- function (v1,v2){
  max(abs(v2-v1))
}
v1 <- c(4,9,18)
v2 <- c(5,7,4)
cheb(v1,v2)

#############
# 6. excerise
#############
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
metric <- ('ManHattan')
v1 <- c(4,9,18)
v2 <- c(5,7,4)
metric <- tolower(metric)
dist (v1,v2,metric)

### Skidam bod jer funkcija koju si napravio (dist) ne radi bez obzira na case metric argumenta.
#############
# 7. excerise
#############
inn_prod <- function (v1,v2){
  sum(v1*v2)
}
v1 <- c(1,3,8)
v2 <- c(4,5,9)
inn_prod(v1,v2)

#############
# 8. excerise
#############
kut <- function (v1,v2){
  acos(sum(v1*v2)/(sqrt(sum(v1*v1))*sqrt(sum(v2*v2))))
}
v1 <- c(3,0)
v2 <- c(5,5)
rad <- kut (v1,v2)
deg <- rad*(180/pi)
deg

#############
# 9. excerise
#############
a <- 1:5
b <- 1:3
d <- c(T,F)

a+b
# rezultat je 2 4 6 5 7
# program zbroji vektore, ali uz upozorenje da duljina vektora a nije vi?ekratnik duljine vektora b
# postupak zbrajanja:
#    a  1 2 3 4 5 
#    b  1 2 3 1 2 <- vektor b je kra?i od vektora a pa se prva dva ?lana vektora b ponavljaju
#    =  2 4 6 5 7

b+d
# rezultat je 2 2 4
# program zbroji vektore, ali uz upozorenje da duljina vektora b nije vi?ekratnik duljine vektora d
# pri zbrajanju program logi?kim podacima iz vektora d pridru?uje 1 ako su TRUE ili 0 ako su FALSE
# postupak zbrajanja:
#    b  1 2 3
#    d  T F T <- vektor d je kra?i od vektora b pa se prvi ?lan vektora d ponavlja
#       1 0 1
#    =  2 2 4

a+d
# rezultat je 2 2 4 4 6
# program zbroji vektore, ali uz upozorenje da duljina vektora a nije vi?ekratnik duljine vektora d
# pri zbrajanju program logi?kim podacima iz vektora d pridru?uje 1 ako su TRUE ili 0 ako su FALSE
# postupak zbrajanja:
#    a  1 2 3 4 5
#    d  T F T F T <- vektor d je kra?i od vektora b pa se vektor d ponavlja 2.5 puta
#       1 0 1 0 1
#    =  2 2 4 4 6

##############
# 10. excerise
##############
# 0 to 5, 0.01 increment, w/o seq
n <- 1:501
a[n] <- 0+((n-1)*0.01)
a
### Kod ce ti bacati gresku ako slucajno nije vec inicijalizirana a (tebi je vjerojatno ostala iz 
### proslog zadatka). To ti vrijedi i za iduci primjer pa cu ti sveukupno za to skinuti 0,5 bodova

#-pi to pi, w/o seq
n <- 1:25
a [n] <- -pi+((n-1)*((2*pi)/24))
a

#0 to 5, 0.01 increment, seq
a <- seq (0,5,0.01)
a

#-pi to pi, seq
a <- seq (-pi,pi,length.out=25)
a

##############
# 11. excerise
##############
### Funkcija ti ne radi jer si zamijenio u prvoj liniji vrijednost koja je proslijedjena kroz 
### argument, a u drugoj liniji radis kao da nisi...
sinus <- function (x=c(0,2*pi)){
  x <- seq(x[1], x[2], 0.01)
  plot(x, sin(x), type="l", xlim=c(x[1], x[2]))
}
sinus(c(1,10))

##############
# 12. excerise
##############
f <- function (t,x){
  n <- 0:t
  zbr <- (((-1)^n)/(factorial((2*n)+1)))*(x^((2*n)+1))
  y <- sum(zbr)
  return(y)
}
f(1,4)

##############
# 13. excerise
##############

##############
# 14. excerise
##############
### Ovo ne radi kako treba za proizvoljne funkcije fun1 i fun2 nego samo za one iz primjera koje 
### si definirao u tijelu funkcije
fun0 <- function (fun1,fun2,int){
  x <- seq(int[1],int[2],length.out=100)
  
  fun1 <- function (x){  
    y1 <- x
    
    fun2 <- function (x){
      y2 <- sin(x)
      matplot(x, cbind(y1,y2),type="l",col=c("red","green"),lty=c(1,1))
      ### Dobra upotreba range()
      legend(range(x)[1], range(y2,y1)[2],c("fun1", "fun2"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","green"))
      dif <- (y1-y2)/y1
      ### Ovo ne uklanja elemente vektora koji poprimaju vrijednost nula nego samo postavlja na nula razliku, ali je 
      ### svejedno uzima u obzir pri racunanju srednje vrijednosti. A i ne radis ovo za opceniti slucaj nego samo za 
      ### funkcije koje sam dao za primjer.
      ### Da si koristio npr. dif <- dif[-1] dobio bi tocnu vrijednost koju sam stavio na forum za testne funkcije
      # dif[1] <- 0
      dif <- dif[-1]
      dif <- mean(dif)*100
      return(dif)
    }
    fun2(x)
  }
  fun1(x)
}

fun0(fun1,fun2,c(0,1))

##############
# 15. excerise
##############
fun0 <- function (fun1,fun2,int){
  x <- seq(int[1],int[2],length.out=100)
  
  fun1 <- function (x){  
    y1 <- factorial(x)
    
    fun2 <- function (x){
      y2 <- (sqrt(2*pi*x)*((x/exp(1))^x))
      matplot(x, cbind(y1,y2),type="l",col=c("red","green"),lty=c(1,1))
      dif <- mean((y1-y2)/y1)*100
      return(dif)
    }
    fun2(x)
  }
  fun1(x)
}

fun0(fun1,fun2,c(2,4))

### Ok, to kazes, ali nisi nikakvim primjerom pokazao...
#na intervalu [2,4] funkcije izgledaju gotovo potpuno jednako
#Stirlingova aproksimacija postaje sve preciznija kako se x pove?ava
#pove?anjem intervala pada prosje?na razlika izme?u funkcija

##############
# 16. excerise
##############
### Sve je ok, osim sto ti kruzici nisu ispunjeni zutom bojom pa -0,5 zbog toga. To si mogao postici recimo s 
### argumentima funkcije points: pch=21, col="red", bg="yellow"
#koordinatni sustav
x_os <- seq(-2, 2, 0.333)
y_os <- seq(-3, 3, 0.5)
plot(xy.coords (x_os, y_os), main="Exercise 16", type="n", xlab="x values",
     ylab="y values", bty="l")
grid(nx=NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"))
axis(2, at=-2)

#linija1
x_lin <- seq(-2, -0.5, 0.5)
y_lin <- seq(-2.5, -1, 0.5)
lines(range(x_lin), range(y_lin), lwd=3)

#to?ke1
x_pnt1 <- seq(-0.5, 0.5, 0.25)
y_pnt1 <- seq(-1, 1, 0.5)
points(xy.coords (x_pnt1,y_pnt1), col="red")

#to?ke2 
x_pnt2 <- seq(0.5, 2, 0.5)
y_pnt2 <- seq(1, 2.5, 0.5)
points(xy.coords(x_pnt2, y_pnt2), cex=0.75, pch=4)

#strelica
arrows(-1, 1, -1.25,-1.25, length=0.15)

#text
text(-1, 1.5, labels="Line width 3", cex=0.75)

##############
# 17. excerise
##############
# Bar-plot se sastoji od stupaca na grafu koji se nalaze iznad 
# oznaka koje predstavljaju kvalitativne varijable. Visina svakog 
# stupca predstavlja veli?inu skupine definirane oznakom. Funkcija
# za crtanje je "barplot".

# Histogram se tako?er sastoji od stupaca na grafu ispod kojih
# su oznake koje ozna?avaju kvantitativne varijable. Oznaka ispod 
# stupca mo?e predstavljati jednu vrijednost ili raspon (range) 
# vrijednosti. Visina stupca predstavlja veli?inu grupe definiranu
# oznakom ispod stupca. Funkcija za crtanje je "hist". 

# Razlika je ?to stupac u bar-plotu predstavlja kvalitativnu
# varijablu, a stupac u histogramu predstavlja kvantitativnu 
# varijablu. 

a <- c(1,3,3,2,4,5,5,6,7,10)
barplot(a, main="barplot", names.arg=c("a","b","c","d","e","f","g","h","i","j"))
hist(a, main="histogram")

##############
# 18. excerise
##############
x <- 1:10
sample(x, size=5)

# Funkcija "sample" uzima nasumi?an uzorak veli?ine n iz vektora 
# duljine m. Ako ne preciziramo "size" uzima nasumi?an uzorak duljine
# jednake duljini vektora (permitira vektor): 

x <- 1:10
sample(x)

# Parametar "replace" dozvoljava da jedan ?lan vektora bude izabran
# vi?e puta:

x <- 1:10
sample(x, replace=TRUE)

#Pomo?u tog vektora mogu?e je uzeti nasumi?an uzorak 
# vektora koji je ve?i od same veli?ine vektora:

x <- 1:10
sample(x, size=15, replace=TRUE)

##############
# 19. excerise
##############
disc <- c(5, 3, 0, 2, 0, 3, 2, 3, 6, 1, 2, 1, 2, 1, 3, 3, 3, 5, 2, 4,
          4, 0, 2, 3, 7, 12, 3, 10, 9, 2, 3, 7, 7, 2, 3, 3, 6, 2, 4, 3, 5, 2, 
          2, 4, 0, 4, 2, 5, 2, 3, 3, 6, 5, 8, 3, 6, 6, 0, 5, 2, 2, 2, 6, 3, 4, 
          4, 2, 2, 4, 7, 5, 3, 3, 0, 2, 2, 2, 1, 3, 4, 2, 2, 1, 1, 1, 2, 1, 4, 
          4, 3, 2, 1, 4, 1, 1, 1, 0, 0, 2, 0);
tbl <- table(disc)
nm <- names(tbl)

# funkcija "table" prebroji koliko se pojedini element (kategorija)
# pojavljuje u vektoru i napravi tablicu u kojoj je svakom
# jedinstvenom elementu pridru?en broj ponavljanja u vektoru. 

# tbl varijabla sadr?i tablicu koja pokazuje koliko je puta u 
# godinama od 1850. do 1959. zabilje?en odre?en broj zna?ajnih 
# otkri?a. Npr. u tom je razdoblju u 12 pojedina?nih godina
# zabilje?eno 1 otkri?e, u 26 pojedin?anih godina 2 otkri?a itd.  

# nm varijabla sadr?i "imena" podataka u tablici, u na?em slu?aju
# imena su brojevi otkri?a (0 otki?a, 1 otkri?e, 2 otkri?a itd.).
# nm je vektor koji sadr?i character data type pa se podaci u njemu
# ne mogu manipulirati broj?anim operatorima (+ - * /). Za dijeljenje
# svakog ?lana nm vektora sa dva podatke u njemu treba prije pretvoriti
# u numeri?ki tip podataka funkcijom as.numeric(nm):

nm_num <- as.numeric(nm)
nm_num <- nm_num/2
nm_num

### Ovo nije tocno, rezultat ti treba biti jednak: mean(disc) tj. 3.1
# ra?unanje prosje?nog broja bitnih otki?a u periodu 1860-1959:
tbl_numeric <- as.numeric(tbl)
tbl_numeric
mean(tbl_numeric)

##############
# 20. excerise
##############
disc2 <- disc[seq(2,length(disc),2)]
disc2
