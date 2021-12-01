library("dplyr")
library("ggplot2")
library("ggthemes")
library('scales') # visualization


############################################################################# 
##### Setting working directory and reading data
############################################################################# 
setwd("C:/Users/fhorvat/Dropbox/Praksa bioinfo/Projekti/test/kaggle")
train <- read.csv(file = "train.csv")
test <- read.csv(file = "test.csv")
full <- bind_rows(train, test)


############################################################################# 
##### Analysing name and titles
############################################################################# 
# Separate title from "Name" column
full$Title <- gsub("(.*, )|(\\..*)", "", full$Name)
table(full$Sex, full$Title)

# Titles with very low cell counts to be combined to "rare" level
rare_title <- c('Dona', 'Lady', 'the Countess','Capt', 'Col', 'Don', 'Dr', 'Major', 'Rev', 'Sir', 'Jonkheer')

# Also reassign mlle, ms, and mme accordingly
full$Title[full$Title == 'Mlle'] <- 'Miss' 
full$Title[full$Title == 'Ms'] <- 'Miss'
full$Title[full$Title == 'Mme'] <- 'Mrs' 
full$Title[full$Title %in% rare_title] <- 'Rare Title'

# Show title counts by sex again
table(full$Sex, full$Title)

# Separate surname from "Name" column
full$Surname <- sapply(full$Name, function(x) strsplit(x, split = '[,.]')[[1]][1])


############################################################################# 
##### Relation between family and surviving
############################################################################# 
# Create a family size variable including the passenger themselves
full$Fsize <- full$SibSp + full$Parch + 1

# Create a family variable 
full$Family <- paste(full$Surname, full$Fsize, sep = '_')

# Use ggplot2 to visualize the relationship between family size & survival (training data)
ggplot(full[1:891,], aes(x = Fsize, fill = factor(Survived))) +
  geom_bar(stat = 'count', position = 'dodge') +
  scale_x_continuous(breaks = c(1:11)) +
  labs(x = 'Family Size') + 
  theme_few()

# Discretize family size
full$FsizeD[full$Fsize == 1] <- 'singleton'
full$FsizeD[full$Fsize < 5 & full$Fsize > 1] <- 'small'
full$FsizeD[full$Fsize > 4] <- 'large'

# Show family size by survival using a mosaic plot
mosaicplot(table(full$FsizeD, full$Survived), main = 'Family Size by Survival', shade = TRUE)

# Adding a deck variable
full$Deck <-factor(sapply(full$Cabin, function(x) strsplit(x, NULL)[[1]][1]))


############################################################################# 
##### Dealing with missing info
############################################################################# 
##### Replacing missing values with a sensible values given the distribution 
##### of the data, e.g., the mean, median or mode
############################################################################# 

# Passengers 62 and 830 are missing Embarkment
full[c(62, 830), 'Embarked']

# Get rid of our missing passenger IDs
embark_fare <- full %>%
  filter(PassengerId != 62 & PassengerId != 830)

# Use ggplot2 to visualize embarkment, passenger class, & median fare
ggplot(embark_fare, aes(x = Embarked, y = Fare, fill = factor(Pclass))) +
  geom_boxplot() +
  geom_hline(aes(yintercept = 80), colour = 'red', linetype = 'dashed', lwd = 2) +
  scale_y_continuous(labels = dollar_format()) +
  theme_few() 

# Since their fare was $80 for 1st class, they most likely embarked from 'C'
full$Embarked[c(62, 830)] <- 'C'

# Passenger on row 1044 has an NA Fare value.
full[1044, ]

# Let’s visualize Fares among all others sharing their class and embarkment (n = 494).
ggplot(full[full$Pclass == '3' & full$Embarked == 'S', ], aes(x = Fare)) +
  geom_density(fill = '#99d6ff', alpha = 0.4) + 
  geom_vline(aes(xintercept = median(Fare, na.rm = T)), colour = 'red', linetype = 'dashed', lwd = 1) +
  scale_x_continuous(labels = dollar_format()) +
  theme_few()

# Replace missing fare value with median fare for class/embarkment
full$Fare[1044] <- median(full[full$Pclass == '3' & full$Embarked == 'S', ]$Fare, na.rm = TRUE)

############################################################################# 
##### Predictive imputation
#############################################################################
# Show number of missing Age values
sum(is.na(full$Age))
