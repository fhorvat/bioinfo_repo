library(ALL)

# A)
data(ALL)
# Ucitava data set iz paketa ALL
# Kada ucitamo neki package koji u sebi ima data.set taj data.set
# se ne uèitava automatski u memoriju zajedno s ostatkom paketa
# nego ga moramo pozvati s data kad ga trebamo. 
# To je napravljeno tako da se velika kolièina podataka ne bi direktno
# uèitavala i nepotrebno zauzimala memoriju. 

exprs(ALL) 
y <- exprs(ALL)
# exprs pristupa podacima spremljenim u ExpressionSet
# ExpressionSet je tip pospremanja podataka dobivenih 
# "high-throughput" istrazivanjima.
# Podaci su spremljeni u matricu i to tako da redovi 
# predstavljaju setove proba, a stupci predstavljaju
# uzorke. 
# U ovom slucaju radi se o microarray podacima o ekspresiji 
# gena za akutnu limafaticnu leukemiju iz data seta ALL. 
# Imena redova imenuju "features" (tocku na microarrayu). 
# Imena stupaca predstavljaju imena uzoraka dobivenih iz pacijenata
# koji boluju od ALL. 

pData(ALL)
pheno.data.ALL <- pData(ALL)
# pData naredba pristupa i daje podatke o uzorcima u ExpressionSetu
# ALL. U ovom slucaju imena redova su imena uzoraka, a
# u stupcima su smjesteni podaci o tom uzorku (pacijentu), npr., 
# dob, spol, datum dijagnoze itd. 


# B)
library(dplyr)

age.sex.F <- pheno.data.ALL %>%
  select(sex, age) %>%
  filter(sex == "F") 
  
  age.sex.M <- pheno.data.ALL %>%
  select(sex, age) %>%
  filter(sex == "M") 

# density dobi pacijenata s obzirom na 
densM <- density(age.sex.M$age, na.rm = T) 
densF <- density(age.sex.F$age, na.rm = T)

# limiti x i y osi, boje
xlim <- range(densM$x, densF$x)
ylim <- range(0, densM$y, densF$y)
densFcol <- rgb(1, 0, 0, 0.2)
densMcol <- rgb(0, 0, 1, 0.2)

# plot
plot(densM, xlim = xlim, ylim = ylim, xlab = 'Age',
     main = 'Distribution of age of female and male ALL patients', 
     panel.first = grid())
polygon(densM, density = -1, col = densMcol)
polygon(densF, density = -1, col = densFcol)

# legenda
legend('topleft',c('M','F'), fill = c(densMcol, densFcol), 
       bty = 'n',border = NA)
