#############
# 1. excerise
#############
# infix operatori su funkcije koje se pozivaju tako da ime funkcije dolazi 
# izmeðu dva argumenta (za razliku od veæine funkcija koje su "prefix" operatori
# i kod kojih ime funkcije dolazi ispred argumenata). 

# %in% provjerava postoji li lijevi argument u desnom i vraæa
# logièku vrijednost (T ako postoji, F ako ne postoji). Primjer:

x <- 1:10
2 %in% x # vraæa TRUE jer 2 postoji u x
11 %in% x # vraæa FALSE jer 11 ne postoji u x
"b" %in% x # vraæa FALSE je "b" ne postoji u x

#############
# 2. excerise
#############
'%my.rbind%' <- function(df1, df2){
  if (is.data.frame(df1) & is.data.frame(df2) == TRUE){
    
    # sluèaj 1: df1 i df2 jednake dužine, imena stupaca = df1
    
    razlika.col <- ncol(df1) - ncol(df2)
    if (razlika.col == 0){ 
      # provjera razlikuju li se imena stupaca df1 i df2:
      df.col.diff <- setdiff(colnames(df1), colnames(df2))
      if (length(df.col.diff) != 0){
        warning("Imena kolona koje se preklapaju nisu jednaka")
      }   
      colnames(df2) <- colnames(df1)
    }
    
    
    # sluèaj 2: df1 duži od df2, imena stupaca = df1
    
    razlika.col <- ncol(df1) - ncol(df2)
    if (razlika.col > 0){
      # provjera razlikuju li se imena stupaca df1 i df2:
      df.col.diff <- setdiff(colnames(df1[1:ncol(df2)]), 
                             colnames(df2))
      if (length(df.col.diff) != 0){
        warning("Imena kolona koje se preklapaju nisu jednaka")
      } 
      imena.col <- colnames(df1)
      for (i in 1 : razlika.col){
        df2 <- cbind(df2, NA)
      }
      colnames(df2) <- imena.col
    }
    
    # sluèaj 3: df2 duži od df1, imena stupaca = poèetak df1, ostatak df2 
    
    razlika.col <- ncol(df2) - ncol(df1)
    if (razlika.col > 0){
      # provjera razlikuju li se imena stupaca df1 i df2:
      df.col.diff <- setdiff(colnames(df1), 
                             colnames(df2[1:ncol(df1)]))
      if (length(df.col.diff) != 0){
        warning("Imena kolona koje se preklapaju nisu jednaka")
      }
      imena.col <- colnames(df1)
      imena.col2 <- colnames(df2[(ncol(df1) + 1) : ncol(df2)])
      for (i in 1 : razlika.col){
        df1 <- cbind(df1, NA)
      }
      colnames(df2) <- c(imena.col, imena.col2)
      colnames(df1) <- c(imena.col, imena.col2)
    }
    return(rbind(df1, df2))
  }
  
  else{
    stop("Error: %my.rbind% radi samo s data.frame argumentima")
  }
}

my.df1 <- data.frame(a1 = 1:3, a2 = 5:7, a3 = c("a", "b", "c"),
                     row.names = 1:3) 
my.df2 <- data.frame(a1 = (1:3)*10, a2 = (5:7)*10, 
                     a3 = c("xx", "yy", "zz"), a4 = c(2:4),
                     row.names = c("xxx", "yyy", "zzz"))

my.df1 %my.rbind% my.df2

#############
# 3. excerise
#############
mat <- matrix(rnorm(30, 5, 2), ncol = 6, nrow = 5)
mat

sd.row <- apply(mat, 1, sd)
sd.row

sd.col <- apply(mat, 2, sd)
sd.col


#############
# 4. excerise
#############
x <- tapply(iris$Sepal.Length, iris$Species, max)
x

#############
# 5. excerise
#############
My.fib <- function(x, n){
  if(x == 1){
    fib <- 1
  }
  else{
    fib <- c(1, 1)
    for (i in 1:x){
      fib[i + 2] <- fib[i] + fib[i + 1]
    }    
  }
  fib <- fib[fib <= x]
  pos.x <- length(fib)
  pos.n <- pos.x - n + 1
  if (n <= pos.x){
    fib.n <- fib[pos.n:pos.x]
  } 
  else{
    fib.n <- rep_len(fib[1:pos.x], length.out = n)
    if (isTRUE(all.equal(n/pos.x, as.integer(n/pos.x))) == FALSE){
      warning("recikliranje nije potpuno")
    } 
  }
  return(fib.n)
}
My.fib(13, 4)

Fibonacciate <- function(v, ncols){
  v <- as.matrix(v)
  mat1 <- t(apply(v, 1, My.fib, n = ncols))
  return(mat1)
}

v1 <- c(13, 17, 132)
Fibonacciate(v1, 4)


#############
# 6. excerise
#############
library(dplyr)
# A)
chr21 <- read.table("002.gtf", stringsAsFactors = F, sep = "\t")

# razdvajanje gene_id, transcript_id stupca na dva odvojena stupca
colnames(chr21) <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
tmp <- matrix(unlist(strsplit(as.character(chr21$"i"), '; ')),
              ncol = 2, byrow = TRUE)
chr21 <- cbind(chr21, as.data.frame(tmp))
chr21 <- select(chr21, -9)
colnames(chr21)[9:10] <- c("gene_id", "transcript_id")

# brisanje "gene_id" i "transcript_id" iz svakog reda u stupcima gene_id i transcript_id
chr21 <- chr21 %>%
  mutate(gene_id = gsub("gene_id ", "", x = chr21$"gene_id")) %>%
  mutate(transcript_id = gsub("transcript_id ", "", x = chr21$"transcript_id"))

# B) 
gene.transc <- read.table("mart_export (1).txt", stringsAsFactors=F, sep="\t",
                          header = TRUE)

# C) All genes that are protein coding
proteins.cod.genes <- gene.transc %>%
  filter(Gene.Biotype == "protein_coding") %>%
  group_by(Ensembl.Gene.ID) %>%
  summarize()

#############
# 7. excerise
#############
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("ALL")

# Paket sadrzi podatke o T- i B-stanicama akutne limfaticne 
# leukemije (ALL) iz Ritz laboratorija na Dana Farber Cancer 
# Institute

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

#podaci za plot
hist.F <- hist(age.sex.F$age, plot = F)                     
hist.M <- hist(age.sex.M$age, plot = F)  

# plot - oba histograma na istom plotu s preklapanjem,
# žene su prikazane rozo, muškarci plavo, a na mjestu preklapanja
# je ljubièasta (semitransparentne boje)

plot(hist.F, col = rgb(1, 0, 1, 0.25), xlim = c(0, 80), 
     ylim = c(0, 25), xlab = "Age",    
     main = "Distribution of age of female and male ALL patients")
plot(hist.M, col = rgb(0, 0, 1, 0.25), xlim = c(0, 80), 
     ylim = c(0, 25), add = T) 
legend('topright', c('M','F'), bty = 'n', border = NA,
       fill = c(rgb(0, 0, 1, 0.25), rgb(1, 0, 0, 0.25)))


# C) 
pheno.data.ALL <- pData(ALL)

# Posložimo tablicu tako da imamo pacijente samo s BCL/ABL 
# leukemijom i gledamo samo relapse i transplant stupce: 
bcr.abl <- pheno.data.ALL %>%
  select(mol.biol, relapse, transplant) %>%
  filter(mol.biol == "BCR/ABL", transplant != "NA") %>% 
  arrange(desc(transplant)) %>% 
  select(-1)

# u posebne tablice odvajam pacijente koji su primili transplant
# od onih koji nisu, grupiram pod relapseu i zbrajam koliko ima
# onih koji jesu u relapseu, a koliko onih koji nisu: 

transplant.TRUE <- bcr.abl %>%
  filter(transplant == "TRUE") %>%
  group_by(relapse) %>%
  summarise(n())

transplant.FALSE <- bcr.abl %>%
  filter(transplant == "FALSE") %>%
  group_by(relapse) %>%
  summarise(n())

# Brojeve stavljam u zajednièku tablicu i to tako da su u 
# stupcima transplant+ i transplant-, a u redovima relapse-
# i relapse+, tablica prikazuje koliko kojih pacijenata ima

transplant.all <- cbind(transplant.TRUE$"n()", transplant.FALSE$"n()")
colnames(transplant.all) <- c("transplant +", "transplant -")
rownames(transplant.all) <- c("relapse -", "relapse +")

# Na toj tablici radim Fisherov test, null hipoteza je da 
# relapse ne ovisi o tome je li primljen transplant 
fisher.test(transplant.all)

# p-vrijednost je 0.045 što znaèi da sa 95.5% sigurnošæu odbacujemo
# null hipotezu, tj. 95.5% je šansa da transplant utjeèe na 
# relapse


# D) 

# Iz pheno.data uzimam samo stupce koji prestavljaju mol.biol tumora
# i cod uzoraka
# Filtriram tablicu tako je u mol.biol samo BCR/ABL leukemija ili
# negativni uzorci. Imenujem redove prema codu s tim da dodajem
# "0" na poèetak:

bcr.abl.neg <- subset(pheno.data.ALL, select = "mol.biol")
bcr.abl.neg <- subset(bcr.abl.neg, subset = (mol.biol == "BCR/ABL" | mol.biol == "NEG"))

# Iz exprs.data uzimam samo one stupce èija imena odgovaraju
# imenima redova u pheno.data za BCR/ABL i negativne uzorke:

exprs.data.ALL <- exprs(ALL)

exprs.data.bcr.abl.neg <- subset(exprs.data.ALL, 
                                 select = (colnames(exprs.data.ALL) %in% rownames(bcr.abl.neg))) 

# exprs.data.bcr.abl.neg ima 111 stupaca što je jednako broju
# pacijenata s BCR/ABL mutacijom ili negativnim uzorkom 


# E) 

# 1. t-test - null hipoteza je da su srednje vrijednosti uzoraka
# jednake, tj. da su uzorci jednaki. p-vrijednost daje vjerojatnost
# da je hipoteza toèna. Vjerojatnost da hipoteza nije toèna, 
# tj. vjerojatnost da su uzorci razlièiti je ((1-p.vrijednost)*100)%

# subset normalnih uzoraka iz exprs.data.ALL
norm.samp <- subset(pheno.data.ALL, select = "mol.biol", 
                    subset = (mol.biol == "NEG"))
norm.samp.exprs <- subset(exprs.data.ALL, 
                          select = (colnames(exprs.data.ALL) %in% rownames(norm.samp)))

# subset patogenih uzoraka iz exprs.data.ALL
gene.samp <- subset(pheno.data.ALL, select = "mol.biol", 
                    subset = (mol.biol == "BCR/ABL"))
gene.samp.exprs <- subset(exprs.data.ALL, 
                          select = (colnames(exprs.data.ALL) %in% rownames(gene.samp))) 

# t-test
# prvo sam napravio transpose data.frame dobivenih iznad kako bi
# mogao primijeniti mapply.
# Nakon mapply u t.test.results matrici imam sve vrijedosti t-testa
# u redovima, a u stupcima su pojedini geni. Zato radim transpose
# i pretvaram u data.frame da bi na njemu mogao raditi s dpylr
# funkcijama. Zatim izbacujem sve stupce iz data.framea osim one 
# s p-vrijednostima. I na kraju u stupac "probability" stavljam 
# vrijednosti vjerojatnosti da je null hipoteza netoèna, tj. 
# vjerojatnost da je ekspresija pojedinog gena razlièita u 
# normalnim i patogenim uzorcima

norm.samp.exprs.t <- as.data.frame(t(norm.samp.exprs))
gene.samp.exprs.t <- as.data.frame(t(gene.samp.exprs))
t.test.results <- mapply(t.test, x = gene.samp.exprs.t, 
                         y = norm.samp.exprs.t, SIMPLIFY = T)
t.test.results <- as.data.frame(t(t.test.results))
t.test.results <- subset(t.test.results, select = "p.value")
t.test.results$probability <- round(as.numeric(t.test.results$p.value), digits = 5)*100

# 2.
a <- as.numeric(t.test.results$"p.value")
p.adj <- p.adjust(a, "fdr")

# Odbacujemo null hipotezu ako je vjerojatnost da podaci odgovaraju 
# null hipotezi mala (p<0.05 npr). Problem je što kada poveæavamo
# broj hipoteza u testu, poveæava se i vjerojatnost da naiðemo na
# neki rijetki dogaðaj zbog kojeg æemo onda odbaciti null hipotezu
# kada je istinita. Kada usporeðujemo dvije skupine u puno razlièitih
# kategorija, vjerojatnost da æemo naiæi na razliku u barem jednoj 
# od tih kategorija se poveæava èistim sluèajem, tj. zato jer ih 
# usporeðujemo više puta. Primjer: promatramo kako lijek djeluje
# na simptome neke bolesti i gledamo koliko je uèinkovit. 
# Što više simptoma promatramo, veæa je šansa da æe lijek ispasti
# uèinkovitiji od postojeæeg lijeka u lijeèenju barem jednog od tih 
# simptoma. 
# U našem sluèaju mi promatramo razliku u ekpresiji gena u puno
# uzoraka, a promatramo ih sve skupa. Postoji šansa da je ekspresija
# gena u jednom od tih uzoraka sluèajno veæa ili manja èisto zato
# jer promatramo puno uzoraka što æe onda utjecati na p-vrijednost
# dobivene t-testom. 
# Zato koristimo ispravak tih p-vrijednosti FDR metodom pomoæu koje
# kontroliramo oèekivani postotak netoèno odbijenih null hipoteza
# ("lažna otkriæa"). Ako je p=0.05 tada oèekujemo da æemo u 5% 
# sluèajeva odbaciti null hipotezu koja je ustvari bila toèna. 
# Pomoæu FDR metode ispravaka p-vrijednosti ta razina netoènih
# odbacivanja null hipoteza se drži na prihvatljivoj razini. 

# 3. 
# moramo odvojiti 50 gena s najveæom razlikom u ekpresiji gena 
# izmeðu patogenih i normalnih uzoraka, tj. gene s najmanjom
# q-vrijednosti

t.test.results$"q.value" <- round(p.adj, digits = 5)
top.50.genes <- subset(t.test.results, select = "q.value")
top.50.genes$gene.names <- rownames(top.50.genes)
top.50.genes <- top.50.genes[, c(2, 1)]
top.50.genes <- arrange(top.50.genes, q.value)
top.50.genes <- head(top.50.genes, 50)

top.50.pheno <- subset(pheno.data.ALL, select = "mol.biol", 
                       subset = (mol.biol == "NEG" | mol.biol == "BCR/ABL"))
exprs.data.top.50 <- subset(exprs.data.ALL, 
                            select = (colnames(exprs.data.ALL) %in% rownames(top.50.pheno)))
top.50.rownames <- row.names(exprs.data.top.50)
exprs.data.top.50 <- as.data.frame(exprs.data.top.50, 
                                   row.names = top.50.rownames)
exprs.data.top.50$gene.names <- top.50.rownames

logical.names <- rownames(exprs.data.top.50) %in% top.50.genes$gene.names
exprs.data.top.50 <- subset(x = exprs.data.top.50,
                            subset = (logical.names))
exprs.data.top.50 <- exprs.data.top.50[-112]


# 4. 

# 5. Heatmap
matrix.top.50 <- data.matrix(exprs.data.top.50)
top.50.heatmap <- heatmap(matrix.top.50, Colv=NA, scale="column")

# heatmapa je naèin prikazivanja podataka pomoæu boja s dendrogramom
# Redovi su poredani tako da imaju srednju vrijednost 0 i 
# standardnu devijaciju 1. 
