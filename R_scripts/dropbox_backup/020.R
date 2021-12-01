v1 <- seq(10, 60, length = 6)
v1

m1 <- matrix(1:30, nrow = 5, ncol = 6, byrow = TRUE)
m1

umn1 <- m1 * v1
umn1

# m1 * v1 - množe se elementi matrice i vektora s obzirom na njihovu poziciju
# Množenje ide po redovima: element prvog reda u prvom stupcu matrice se 
# množi s prvim elementom vektora, element drugog reda u prvom stupcu s 
# drugim elementom vektora itd. S obzirom da matrica ima 5 redova, a vektor 
# 6 elemenata, šesti element vektora se množi s prvim redom drugog stupca.
# Zatim se opet prvi element vektora množi s elementom drugog reda i drugog 
# stupca matrice itd. sve dok se svi elementi matrice ne pomnože s vektorom
# (recikliranje)

umn2 <- m1 %*% v1
umn2
m1[2,]*v1
sum(m1[2,]*v1)


# m1 %*% v1 - promatra vektor kao matricu s jednim stupcem i 6 redova,
# zatim množi matricu m1 s matricom v1:

#      [,1] [,2] [,3] [,4] [,5] [,6]        [,1]
# [1,]    1    2    3    4    5    6   [1,]   10
# [2,]    7    8    9   10   11   12   [2,]   20
# [3,]   13   14   15   16   17   18 * [3,]   30
# [4,]   19   20   21   22   23   24   [4,]   40
# [5,]   25   26   27   28   29   30   [5,]   50
#                                      [6,]   60
# 
# Množe se redovi sa stupcem na ovaj naèin:
# 1*10 + 2*20 + 3*30 + 4*40 + 5*50 + 6*60 = 910
# 7*10 + 8*20 + 9*30 + 10*40 + 11*50 + 12*60 = 2170 
# tj. prvi element stupca s prvim elementom reda, drugi element s drugim itd.
# rezultat je matrica s jednim stupcem i 5 redova. 

umn3 <- v1 %*% m1

# nije moguæe obaviti operaciju. 
# v1 se promatra kao matrica s jednim redom i 6 stupaca (6 elemenata), 
# a svaki stupac matrice m1 ima 5 redova (5 elemenata). Algoritam treba 
# pomnožiti red matrice v1 sa svakim stupcem matrice m1, ali su razlièite
# duljine (non-conformable) pa to nije moguæe. Kad bi npr. vektor v1 imao 
# 5 elemenata tada bi v1 %*% m1 bilo moguæe:

v1 <- seq(10, 50, length=5)
v1
umn3 <- v1 %*% m1
umn3
