
# This is a preliminary analysis plan for the PSA 003 study

####################################################
#              Load packages                       #
####################################################
library(cocor) # compare two correlations to each other 
library(MASS) # for mvrnorm

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


### Approach 2: comparing the correlations of AW_B and AWM_BM vs. the correlation of AW_B and AWF_BF. We assume that the corelation of AW_B and AWM_BM would be higher than the correlation of AW_B and AWF_BF.

# test whether there is a significant difference between the correlation of AW_B with AWM_BM and AWF_BF
# (this is a one-sided test expecting the correlation of AW_B with AWM_BM to be greater than that of AWF_BF)
# the same analysis could be done using paired.r from the psych package, but it would not give us confidence intervals for the difference in correlations


dataset_nomising = as.data.frame(na.omit(cbind(dataset$AW_B, dataset$AWM_BM, dataset$AWF_BF)))
# name columns
names(dataset_nomising) = c("AW_B", "AWM_BM", "AWF_BF")



### simulating data
N_sim_per_group = 89

# NORTHAM data
correlation_matrix_NORTHAM = cor(dataset_nomising)
means_NORTHAM = apply(dataset_nomising, 2, mean)
SDs_NORTHAM = apply(dataset_nomising, 2, sd)
country = "NORTHAM"


# creating standardized data with the right correlation structure determined in the above section
sim_data_pre_NORTHAM = mvrnorm(N_sim_per_group, mu = rep(0, ncol(correlation_matrix_NORTHAM)),
                               Sigma = correlation_matrix_NORTHAM,
                               empirical = F)

# duplicate sim_data_pre and fill it up with correct means and standard deviations matching the original dataset
sim_data_NORTHAM = as.data.frame(sim_data_pre_NORTHAM)
for(i in 1:ncol(correlation_matrix_NORTHAM)){
  sim_data_NORTHAM[,i] = sim_data_pre_NORTHAM[,i] * SDs_NORTHAM[i] + means_NORTHAM[i]
}


sim_data_NORTHAM = cbind(sim_data_NORTHAM, country)

# test whether there is a significant difference between the correlation of AW_B with AWM_BM and AWF_BF
# (this is a one-sided test expecting the correlation of AW_B with AWM_BM to be greater than that of AWF_BF)
# the same analysis could be done using paired.r from the psych package, but it would not give us confidence intervals for the difference in correlations
cor_results_NORTHAM = cocor(~AW_B + AWM_BM | AW_B + AWF_BF, alternative = "greater", data = sim_data_NORTHAM)
# results
cor_results_NORTHAM


# ASIA data (simulating no difference in correlation in this region)


correlation_matrix_ASIA = matrix(c(1, 0.77, 0.77,
                                   0.77, 1, 0.82,
                                   0.77, 0.82, 1),
                                 ncol = 3, byrow = TRUE)
rownames(correlation_matrix_ASIA) = rownames(correlation_matrix_NORTHAM)

means_ASIA = apply(dataset_nomising, 2, mean)
SDs_ASIA = apply(dataset_nomising, 2, sd)
country = "ASIA"


# creating standardized data with the right correlation structure
sim_data_pre_ASIA = mvrnorm(N_sim_per_group, mu = rep(0, ncol(correlation_matrix_ASIA)),
                            Sigma = correlation_matrix_ASIA,
                            empirical = F)

# duplicate sim_data_pre and fill it up with correct means and standard deviations matching the original dataset
sim_data_ASIA = as.data.frame(sim_data_pre_ASIA)
for(i in 1:ncol(correlation_matrix_ASIA)){
  sim_data_ASIA[,i] = sim_data_pre_ASIA[,i] * SDs_ASIA[i] + means_ASIA[i]
}

sim_data_ASIA = cbind(sim_data_ASIA, country)


### Approach 2: comparing the correlations of AW_B and AWM_BM vs. the correlation of AW_B and AWF_BF. We assume that the corelation of AW_B and AWM_BM would be higher than the correlation of AW_B and AWF_BF.

# statistical test
cor_results_ASIA = cocor(~AW_B + AWM_BM | AW_B + AWF_BF, alternative = "greater", data = sim_data_ASIA)
# results
cor_results_ASIA

####### Determining region effect (prejudice is gendered in North america, but not in Asia)

### Calculating residuals for the NORTHAM region

mod_male = lm(AW_B ~ AWM_BM, data = sim_data_NORTHAM)
mod_female = lm(AW_B ~ AWF_BF, data = sim_data_NORTHAM)

AWM_BM_error = (sim_data_NORTHAM$AW_B - predict(mod_male))^2
AWF_BF_error = (sim_data_NORTHAM$AW_B - predict(mod_female))^2

# differences in the residual error in the two models
sim_data_NORTHAM$Diff_error = AWF_BF_error - AWM_BM_error


### Calculating residuals for the ASIA region

mod_male = lm(AW_B ~ AWM_BM, data = sim_data_ASIA)
mod_female = lm(AW_B ~ AWF_BF, data = sim_data_ASIA)

AWM_BM_error = (sim_data_ASIA$AW_B - predict(mod_male))^2
AWF_BF_error = (sim_data_ASIA$AW_B - predict(mod_female))^2

# differences in the residual error in the two models
sim_data_ASIA$Diff_error = AWF_BF_error - AWM_BM_error



### Establishing Region effects
sim_data_COMBINED = rbind(sim_data_NORTHAM, sim_data_ASIA)

summary(lm(Diff_error ~ country, data = sim_data_COMBINED))
