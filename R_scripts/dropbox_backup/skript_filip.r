library(ggplot2)

setwd("C:/Users/toonka/Desktop/simple_graph")
data <- read.csv(file="Biomasa.csv")

# prvo transformiras brojeve u tablici (podjelis ih s 10000000 da ih dovedes na isti red velicine kao i drugu varijablu)
# dakle falio ti je ovaj red:
data$Number_of_cell_per_L <- data$Number_of_cell_per_L / 10000000

p <- ggplot(data, aes(x = sample, group=1))
p <- p + geom_line(aes(y = Biomass_mgL.1, colour = "Biomass_mgL/L"))

p <- p + geom_line(aes(y = Number_of_cell_per_L, colour = "Number_of_cells/L")) 

# tu mu kazes da doda drugu os i da oznake pomnozi s istim brojem s kojim se gore dijelila da bi se vratili na originalne vrijednosti
# ustvari su prikazani podaci i dalje podjeljeni s 10000000, samo su oznake na osi promijenjene, ali to je za tebe nebitno
p <- p + scale_y_continuous(sec.axis = sec_axis(~ .* 10000000, name = "Number_of_cells/L"))


p <- p + scale_colour_manual(values = c("blue", "red"))
p <- p + labs(y = "Biomass_mg/L",
              x = "Sample",
              colour = "Parameter")+
  scale_x_discrete(limits=c("S06","S07","S08","S09", "S10", "S11", "S12", "S01", "S02", "S03"))
p <- p + theme(legend.position = c(0.8, 0.9))
p  <- p+ theme_classic()
p  


ggsave("p.png", p, width = 8, height =6 )
