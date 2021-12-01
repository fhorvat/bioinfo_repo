  # Normal distribution:
par(mfrow = c(1, 1))
rm(list=ls())

Distribution <- function(n){
  rnd.x <- list()
  mean.rnd.x <- list()
  for (i in 1:length(n)){
      rnd.x[i] <- list(matrix(rnorm(n[i]*10000), nrow = n[i], ncol = 10000, 
                      byrow = FALSE))
      mean.rnd.x[i] <- list(colMeans(rnd.x[[i]]))
    }
  stn.dev <- 0
  for (j in 1:length(n)){
    stn.dev[j] <- sd(unlist(mean.rnd.x[j]))
  }
  mean.rnd.x$standard.deviation <- data.frame(n, stn.dev)
  plot(stn.dev, n, type = 'l')
  return(mean.rnd.x)
}
set.seed(1)
Distribution(c(5, 50, 500, 5000))

  # plot distribution:
a <- Distribution(c(5, 50, 500, 5000))
par(mfrow = c(2, 2))
hist(a[[1]])
hist(a[[2]])
hist(a[[3]])
hist(a[[4]])
par(mfrow = c(1, 1))

  # standard deviation
  # Standardna devijacija je prikazana u data.frame-u na zadnjem mjestu 
  # u listi koju vraæa funkcija Nor.dist. Što je veæi n, standardna
  # devijacija je manja (obrnuto proporcionalno)


  # Uniformna distribucija 
par(mfrow = c(1, 1))
rm(list=ls())

Distribution <- function(n){
  rnd.x <- list()
  mean.rnd.x <- list()
  for (i in 1:length(n)){
    rnd.x[i] <- list(matrix(runif(n[i]*10000), nrow = n[i], ncol = 10000, 
                            byrow = FALSE))
    mean.rnd.x[i] <- list(colMeans(rnd.x[[i]]))
  }
  stn.dev <- 0
  for (j in 1:length(n)){
    stn.dev[j] <- sd(unlist(mean.rnd.x[j]))
  }
  mean.rnd.x$standard.deviation <- data.frame(n, stn.dev)
  plot(stn.dev, n, type = 'l')
  return(mean.rnd.x)
}
set.seed(1)
Distribution(c(5, 50, 500, 5000))

# plot distribution:
a <- Distribution(c(5, 50, 500, 5000))
par(mfrow = c(2, 2))
hist(a[[1]])
hist(a[[2]])
hist(a[[3]])
hist(a[[4]])
par(mfrow = c(1, 1))


  # Poissonova distribucija 
par(mfrow = c(1, 1))
rm(list=ls())

Distribution <- function(n){
  rnd.x <- list()
  mean.rnd.x <- list()
  for (i in 1:length(n)){
    rnd.x[i] <- list(matrix(rpois(n[i]*10000, l = 2), nrow = n[i], 
                            ncol = 10000, byrow = FALSE))
    mean.rnd.x[i] <- list(colMeans(rnd.x[[i]]))
  }
  stn.dev <- 0
  for (j in 1:length(n)){
    stn.dev[j] <- sd(unlist(mean.rnd.x[j]))
  }
  mean.rnd.x$standard.deviation <- data.frame(n, stn.dev)
  plot(stn.dev, n, type = 'l')
  return(mean.rnd.x)
}
set.seed(1)
Distribution(c(5, 50, 500, 5000))

  # plot distribution:
a <- Distribution(c(5, 50, 500, 5000))
par(mfrow = c(2, 2))
hist(a[[1]])
hist(a[[2]])
hist(a[[3]])
hist(a[[4]])
par(mfrow = c(1, 1))


  # Moja distribucija 
par(mfrow = c(1, 1))
rm(list=ls())

  # vraæa n brojeva po distribuciji, brojevi idu [-2*pi, 2*pi], 
  # a vjerojatnosti od [(sin(-2*pi))^2, (sin(2*pi))^2]
My.random.gen <- function(n){
  rnd.gen <- seq(-2*pi, 2*pi, length.out = 10000)
  rnd.gen.probability <- (sin(rnd.gen))^2
  xyz <- sample(rnd.gen, n, prob = rnd.gen.probability, replace = TRUE)
  return(xyz)
}

  # plot moje distribucije, napravio sam ga u posebnoj funkciji da se ne
  # bi kasnije kad pozovem funkciju unutar funkcije Distribution i puta
  # iscrtavao plot:
My.random.gen.plot <- function(n){
  rnd.gen <- seq(-2*pi, 2*pi, length.out = 10000)
  rnd.gen.probability <- (sin(rnd.gen))^2
  plot(rnd.gen, rnd.gen.probability, type = "l")
}
My.random.gen.plot(1)

Distribution <- function(n){
  rnd.x <- list()
  mean.rnd.x <- list()
  for (i in 1:length(n)){
    rnd.x[i] <- list(matrix(My.random.gen(n[i]*10000), nrow = n[i], 
                            ncol = 10000, byrow = FALSE))
    mean.rnd.x[i] <- list(colMeans(rnd.x[[i]]))
  }
  stn.dev <- 0
  for (j in 1:length(n)){
    stn.dev[j] <- sd(unlist(mean.rnd.x[j]))
  }
  mean.rnd.x$standard.deviation <- data.frame(n, stn.dev)
  plot(stn.dev, n, type = 'l')
  return(mean.rnd.x)
}
set.seed(1)
Distribution(c(5, 50, 500, 5000))

  # plot distribution:
a <- Distribution(c(5, 50, 500, 5000))
par(mfrow = c(2, 2))
hist(a[[1]])
hist(a[[2]])
hist(a[[3]])
hist(a[[4]])
par(mfrow = c(1, 1))
