x <- sum(rep(c(T, T, F), times = 5))
x

a <- rep(c(T, T, F), times = 5)
a

b <- sum(c(T, T, F))
b
  
  # - rezultat sum(rep(c(T, T, F), times = 5)) je 10
 
  # OBAŠNJENJE:

  # - rep sekvencu T T F replicira 5 puta
  
  # - nakon rep(c(T, T, F), times = 5) imamo T T F sekvencu 5 puta tj. 
  # 10 puta T i 5 puta F

  # - sum radi na numeric argumentima

  # - kad pomoæu sum pokušamo zbrojiti logic argumente R ih pretvori u 
  # numeric tako da je TRUE = 1, a FALSE = 0
  
  # - kako je T = 1, a F = 0, kada sekvencu (T T F) x 5 koja ima 10xT i 5xF
  # zbrojimo pomoæu sum imamo 10*1 + 5*0 što daje 10


