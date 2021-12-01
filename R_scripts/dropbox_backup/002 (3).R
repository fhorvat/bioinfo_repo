  # a)
chilled.Mc1 <- CO2[CO2$Plant == "Mc1" & CO2$Treatment == "chilled", ]
chilled.Mc1

  # b)
chilled.Q.less250CO2 <- CO2[CO2$Type == "Quebec" & CO2$Treatment == "chilled"
                            & CO2$conc < 250, ]
chilled.Q.less250CO2

  # c)
chilled.Q.less250CO2.odd <- CO2[c(TRUE, FALSE) & CO2$Type == "Quebec" 
                                & CO2$Treatment == "chilled"
                                & CO2$conc < 250,  ]
chilled.Q.less250CO2.odd

  # d) 
more350CO2.more35up <- CO2[CO2$conc > 350 & CO2$uptake > 35, ] 
more350CO2.more35up

  # e) 
more350CO2.more35up.plant.type <- CO2[CO2$conc > 350 & CO2$uptake > 35,
                                      c(1, 2)]
more350CO2.more35up.plant.type

  # f)
  # a)
sub.chilled.Mc1 <- subset(CO2, Plant == "Mc1" & Treatment == "chilled")
sub.chilled.Mc1

  # b)
sub.chilled.Q.less250CO2 <- subset (CO2, Type == "Quebec" 
                                    & Treatment == "chilled" & conc < 250)
sub.chilled.Q.less250CO2 

  # c) 
sub.chilled.Q.less250CO2.odd <- subset(CO2, c(TRUE, FALSE) & Type == "Quebec" 
                                       & Treatment == "chilled" & conc < 250)
sub.chilled.Q.less250CO2.odd 

  # d) 
sub.more350CO2.more35up <- subset(CO2, conc > 350 & uptake > 35)
sub.more350CO2.more35up

  # e)
sub.more350CO2.more35up.plant.type <- subset(CO2, conc > 350 & uptake > 35, 
                                             c(Plant, Type))
sub.more350CO2.more35up.plant.type
