
####################################################
#              Load packages                       #
####################################################

library(lme4) # approach 1, for mixed effects anova
library(lmerTest) # approach 1, for mixed effects anova (to do significance test)
library(cocor) # approach 2, compare two correlations to each other 

####################################################
#              Load data file                      #
####################################################

# load pilot data
dataset = read.csv("https://raw.githubusercontent.com/kekecsz/Gendered-nature-of-prejudice/master/dataset_pilot.csv")

###########################################################
#              Statistical analysis                       #
###########################################################

### Approach 1: ANOVA looking for race * sex interaction

## data management

# creating a long format from a wide format
warmth = c(dataset$Therm_8, dataset$Therm_9, dataset$Therm_14, dataset$Therm_15)
race = rep(c("black", "white"), each = 2*length(dataset$Therm_8))
sex = c(rep(c("men", "women"), each = length(dataset$Therm_8)), rep(c("men", "women"), each = length(dataset$Therm_8)))
ID = rep(paste("ID_", 1:length(dataset$Therm_8), sep = ""), 4)
data = as.data.frame(warmth)
data = as.data.frame(cbind(data, ID, race, sex))

# specifying factors
data$ID = factor(data$ID)
data$race = factor(data$race)
data$sex = factor(data$sex)

# statistical test
mod = lmer(warmth ~ race * sex + (1|ID), data = data)
# results
summary(mod)



### Approach 2: comparing the correlations of AW_B and AWM_BM vs. the correlation of AW_B and AWF_BF. We assume that the corelation of AW_B and AWM_BM would be higher than the correlation of AW_B and AWF_BF.

# statistical test
cor_results = cocor(~AW_B + AWM_BM | AW_B + AWF_BF, alternative = "greater", data = dataset)
# results
cor_results


# the same analysis with paired.r from the psych package
dataset_nomising = as.data.frame(na.omit(cbind(dataset$AW_B, dataset$AWM_BM, dataset$AWF_BF)))
names(dataset_nomising) = c("AW_B", "AWM_BM", "AWF_BF")

AW_B_AWM_BM = cor(dataset_nomising$AW_B, dataset_nomising$AWM_BM)
AW_B_AWF_BF= cor(dataset_nomising$AW_B, dataset_nomising$AWF_BF)
AWM_BM_AWF_BF= cor(dataset_nomising$AWM_BM, dataset_nomising$AWF_BF)


paired.r(AW_B_AWM_BM, AW_B_AWF_BF, AWM_BM_AWF_BF, n = nrow(dataset_nomising), twotailed = F)
