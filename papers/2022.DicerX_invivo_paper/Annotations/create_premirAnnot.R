library(data.table)

MIRANNOTDT <- load("mirAnnot.dt.rda",verbose=T)
MIRANNOTDT <- get(MIRANNOTDT)

MIRANNOTDT[ ID %in% c("MIMAT0049839_1","MIMAT0049840_1") ][["preID"]] <- "MI0040626_1"
MIRANNOTDT[ ID %in% c("MIMAT0014934",  "MIMAT0014935")   ][["preID"]] <- "MI0014099.2"

premirAnnot.dt <-
 MIRANNOTDT[ , {

  list(
   ID.5p        =       ID[ miRNA.STR=="5p"      ],
   ID.3p        =       ID[ miRNA.STR=="3p"      ],
   ID.unknown   =       ID[ miRNA.STR=="unknown" ],
   Name.5p      =     Name[ miRNA.STR=="5p"      ],
   Name.3p      =     Name[ miRNA.STR=="3p"      ],
   Name.unknown =     Name[ miRNA.STR=="unknown" ],
   type.5p      = pep.type[ miRNA.STR=="5p"      ],
   type.3p      = pep.type[ miRNA.STR=="3p"      ],
   type.unknown = pep.type[ miRNA.STR=="unknown" ]
  )
                } , by="preID" ]

save( premirAnnot.dt,file="premirAnnot.dt.rda" )

