mydata <- read.table(file = 'hgnc_complete_set.txt', sep = '\t', 
                     header = T, fill = TRUE, quote = "", 
                     stringsAsFactors = FALSE)
mydata <- mydata[, -4:-17]
mydata <- mydata[, -6:-26]
a <- grepl('withdrawn', mydata$Approved.Symbol)
mydata <- subset(mydata, subset = !a)


# Funkcija za usporeğivanje (a < z) 
comp <- function(x, y){
  if (x == y){
    return(0)
  }
  if (x < y) {
    return(1)
  }
  if (x > y) {
    return(2)
  }
}


# Dodavanje na binary search tree (desno = veæi, lijevo = manji)
My.b.s.tree.add <- function(tree, new.item){
  
  new.node <- function(new.item){
    node <- as.list(new.item)
    node <- c(list(left = NULL, right = NULL), node)
    return(node)
  }
  
  if (is.null(tree)){
    return(new.node(new.item))
  }
  
  child <- tree
  path <- numeric(0)
  
  while (TRUE){
    comp.value <- comp(new.item$Approved.Symbol, child$Approved.Symbol)
    if (comp.value == 0){
      return(tree)
    }
    path <- c(path, comp.value)
    child <- child[[comp.value]]
    if (is.null(child) == TRUE){
      tree[[path]] <- new.node(new.item)
      return(tree)
    }
  }
}

my.tree <- NULL
my.tree <- My.b.s.tree.add(my.tree, mydata[3, ])


# Brisanje s binary search tree
# Nisam stigao napraviti brisanje kada node ima potomke na obje strane, 
# ali svi ostali sluèajevi rade (kada node nema potomaka, ili kada node
# ima potomka (potomke) samo na jednoj strani)
My.b.s.tree.rm <- function(tree, appr.symb){
  
  path <- numeric(0)
  tree.tmp <- tree
  
  while (TRUE){
    
    if (length(tree.tmp$Approved.Symbol) < 1){
      return("Traeni simbol ne postoji u stablu")
    }
    
    comp.value <- comp(appr.symb, tree.tmp$Approved.Symbol)
    if (comp.value == 0){
      
      # Brisanje kada nema potomaka
      if (is.null(tree.tmp$left) & is.null(tree.tmp$right)){
        tree[[path]] = NULL
        return(tree)
      }
      
      # Brisanje kada je potomak (potomci) desno 
      if (!is.null(tree.tmp$right) & is.null(tree.tmp$left)){
        tree.right <- tree[[c(path, 2)]]
        if (length(path) > 0){
          tree[[path]] <- tree.right
        }
        else{
          tree <- tree.right
        }
        return(tree)
      }
      
      # Brisanje kada je potomak (potomci) lijevo 
      if (!is.null(tree.tmp$left)& is.null(tree.tmp$right)){
        tree.left <- tree[[c(path, 1)]]
        if (length(path) > 0){
          tree[[path]] <- tree.left
        }
        else{
          tree <- tree.left
        }
        return(tree)
      }
      
      warning("Ima potomke u obje grane")
      return(tree)
    }
    path <- c(path, comp.value)
    tree.tmp <- tree.tmp[[comp.value]] 
    if (is.atomic(tree.tmp)){
      warning("Traeni gen ne postoji u stablu")
      return(tree)
    }
  }  
}

my.tree <- My.b.s.tree.rm(my.tree, "A1BG-AS1")
my.tree


# Traenje na binary search tree
My.b.s.tree.get <- function(tree, appr.symb){
  while (TRUE){
    comp.value <- comp(appr.symb, tree$Approved.Symbol)
    if (comp.value == 0){
      value.appr.symb <- tree[3:7]
      return(value.appr.symb)
    }
    tree <- tree[[comp.value]]
    if (is.atomic(tree)){
      return(warning("Traeni gen"))
    }
  }
}

binary.get <- My.b.s.tree.get(my.tree, "A2M")
binary.get



# Sekvencijalno dodavanje na listu 
My.seq.add <- function(my.seq, new.item){
  
  if (length(my.seq) == 0){
    my.seq <- list()
    my.seq[[1]] <- as.list(new.item)
    return(my.seq)
  }
  
  for (i in 1:length(my.seq)){
    if (my.seq[[i]]$Approved.Symbol == new.item$Approved.Symbol){
      return(my.seq)
    }
  }
  my.seq[[length(my.seq)+1]] <- as.list(new.item)
  return(my.seq)
}

my.seq <- NULL
my.seq <- My.seq.add(my.seq, mydata[3, ])
my.seq


# Sekvencijalno brisanje 
My.seq.rm <- function(my.seq, appr.symb){
  for (i in 1:length(my.seq)){
    if (my.seq[[i]]$Approved.Symbol ==  appr.symb){
      my.seq[[i]] <- NULL
      return(my.seq)
    }   
  }
  warning("Traeni gen ne postoji na listi")
  return(my.seq)
}

my.seq <- My.seq.rm(my.seq, "A")
my.seq


# Sekvencijalno traenje 
My.seq.get <- function(my.seq, appr.symb){
  for (i in 1:length(my.seq)){
    if (my.seq[[i]]$Approved.Symbol ==  appr.symb){
      return(my.seq[[i]])
    }   
  }
  return(warning("Traeni gen ne postoji u listi"))
}

seq.get <- My.seq.get(my.seq, "A1BG")
seq.get


# Testiranje brzine izvoğenja
library(microbenchmark)

# Prvo randomiziram data.frame i u svoje drvo stavljam prvih 10000 podataka 
# iz tog data.framea
my.random.data <- mydata[sample(nrow(mydata)), ]

my.tree <- NULL
for (i in 1:1000){
  my.tree <- My.b.s.tree.add(my.tree, my.random.data[i, ])
}

# To isto napravim za sekvencijalnu listu:
my.seq <- NULL
for (i in 1:1000){
  my.seq <- My.seq.add(my.seq, my.random.data[i, ])  
}

# Testiranje brzine dodavanja My.b.s.tree.add i My.seq.add
microbenchmark(My.b.s.tree.add(my.tree, my.random.data[1001, ]),
               My.seq.add(my.seq, my.random.data[1001, ]))

# Testiranje brzine brisanja My.b.s.tree.rm i My.seq.rm
microbenchmark(My.b.s.tree.rm(my.tree, "DMRTA2"),
               My.seq.rm(my.seq, "DMRTA2"))

# Testiranje brzine dohvaæanja My.b.s.tree.get i My.seq.get
microbenchmark(My.b.s.tree.get(my.tree, "DMRTA2"),
               My.seq.get(my.seq, "DMRTA2"))

# My.b.s.tree.* funkcije su znaèajno bre od My.seq.*. 

# Najveæa razlika je kod add funkcija, My.b.s.tree.add median je oko 700 ms,
# a My.seq.add median je oko 19000 ms (My.b.s.tree.add je oko 25 puta bra).

# My.b.s.tree.rm median je oko 80 ms, a My.seq.rm median je oko 650 
# (My.b.s.tree.rm je oko 8 puta bra).

# Get funkcije su najbre, a opet je My.b.s.tree.get (median oko 60) bra od
# My.seq.get (median oko 650). My.b.s.tree.get je oko 10 puta bra. 


# Big O notation

# My.b.s.tree.* funkcije imaju prosjeènu vremensku kompleksnost O(log n) 
# za sve funkcije (dodavanje, brisanje i dohvaæanje). Kada funkcija u
# binary search tree trai element u svakoj grani ima dva izbora, kada
# funkcija izabere koju granu æe slijediti polovica stabla je odbaèena.
# Tj. u svakom koraku prepolovimo problem, dakle kompleksnost je log(n). 
# Kada konstruiramo stablo svaki put kada udvostruèimo broj èvorova poveæamo
# broj koraka do rješenja za 1.

# My.seq.* funkcije imaju prosjeènu vremensku kompleksnost O(n) za sve 
# funkcije. Da bismo dohvatili n-ti element u povezanoj listi prvo moramo 
# proæi svih n-1 elemenata pa je zato vremenska kompleksnost O(n). Isto tako
# za dodavanje i brisanje, da bismo došli do n-tog elementa trebamo prvo 
# proæi n-1 elemenata prije toga. 