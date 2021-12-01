v1 <- seq(10, 60, length = 6)
v1

m1 <- matrix(1:30, nrow = 5, ncol = 6, byrow = TRUE)
m1

zbr1 <- m1 + v1
zbr1

# m1 * v1 - zbrajaju  se elementi matrice i vektora s obzirom na njihovu poziciju
# Zbrajanje ide po redovima: element prvog reda u prvom stupcu matrice se 
# zbraja s prvim elementom vektora, element drugog reda u prvom stupcu s 
# drugim elementom vektora itd. S obzirom da matrica ima 5 redova, a vektor 
# 6 elemenata, šesti element vektora se zbraja s prvim redom drugog stupca.
# Zatim se opet prvi element vektora zbraja s elementom drugog reda i drugog 
# stupca matrice itd. sve dok se svi elementi matrice ne zbroje s vektorom
# (recikliranje)
