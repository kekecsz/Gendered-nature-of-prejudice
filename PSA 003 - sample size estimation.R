library(MASS) # for mvrnorm
library(cocor) # for cocor

# function doing statistical inference
statistical_inference = function(data, minimal_effect_of_interest, CI_width){
  
  corr_table <- cor(data)
  
  results <- get.cocor.results(cocor.dep.groups.overlap(corr_table[2,1], 
                                                        corr_table[3,1], 
                                                        corr_table[3,2],
                                                        n = nrow(data),
                                                        conf.level = CI_width))
  
  CI_lb = results$zou2007$conf.int[1]
  CI_ub = results$zou2007$conf.int[2]
  
  inference = if(((CI_lb > 0)|(CI_ub < 0)) & (CI_lb > minimal_effect_of_interest*-1) & (CI_ub < minimal_effect_of_interest)){"null_or_too_small_effect"
  } else if((CI_lb > minimal_effect_of_interest*-1) & (CI_ub < minimal_effect_of_interest)){"null_or_too_small_effect"
  } else if(CI_lb > 0){"greater"
  } else if (CI_ub < 0){"less"
  } else {"inconclusive"}
  
  return(inference)
}



# function simulating new samples and checking inference
simulation_function = function(N_per_subgroup, true_effect, 
                               baseline_corr_matrix,
                               sides_of_test, minimal_effect_of_interest,
                               CI_width){
  data = list(NA)
  
  for(i in 1:length(true_effect)){
    correletion_modifier = t(matrix(c(0, true_effect[i]/2, true_effect[i]/2*-1,
                                      true_effect[i]/2, 0, 0,
                                      true_effect[i]/2*-1, 0, 0), nrow = 3))
    
    cor_matrix = baseline_corr_matrix + correletion_modifier
    
    #generate data
    data[[i]] = mvrnorm(n = N_per_subgroup, mu = c(0, 0, 0), Sigma = cor_matrix, empirical = F)
  }
  
  if(sides_of_test == "one.sided"){
    correct_inferences = sapply(true_effect, function(x) if(((x < 0)|(x > 0)) & (x < minimal_effect_of_interest) &  (x > (minimal_effect_of_interest*-1))){"null_or_too_small_effect"} else if((x < minimal_effect_of_interest) &  (x > (minimal_effect_of_interest*-1))){"null_or_too_small_effect"} else if(x < 0){"less"} else if(x > 0){"greater"})
    inferences = unlist(lapply(data, function(x) statistical_inference(data = x, minimal_effect_of_interest = minimal_effect_of_interest, CI_width = CI_width)))
  } else if(sides_of_test == "two.sided"){
    correct_inferences = sapply(true_effect, function(x) if(((x < 0)|(x > 0)) & (x < minimal_effect_of_interest) &  (x > (minimal_effect_of_interest*-1))){"null_or_too_small_effect"} else if((x < minimal_effect_of_interest) &  (x > (minimal_effect_of_interest*-1))){"null_or_too_small_effect"} else if((x < 0)|(x > 0)){"different"})
    inferences_pre = unlist(lapply(data, function(x) statistical_inference(data = x, minimal_effect_of_interest = minimal_effect_of_interest, CI_width = CI_width)))
    inferences = sapply(inferences_pre, function(x) if((x == "less") | (x == "greater")){"different"} else {x})
  }
  
  inferences_not_inconclusive = inferences[inferences != "inconclusive"]
  correct_inferences_where_inference_is_not_inconclusive = correct_inferences[inferences != "inconclusive"]
  
  inference_correct = if(all(inferences == correct_inferences)){"correct inference"} else if(any(inferences_not_inconclusive != correct_inferences_where_inference_is_not_inconclusive)){"false inference"} else {"inconclusive"}
  
  return(inference_correct)
}



# Number of simulations to run
iterations = 10000
N_per_subgroup = 2300
final_results = as.data.frame(matrix(NA, nrow = 4, ncol = 4))
names(final_results) = c("scenario", "all correct", "any false", "any inconclusive")
CI_width = 0.999

########### Testing both hypothesis 1, 2, and 3

## Scenario 1. simulating the originally expected pattern of effects

simulation_output = replicate(iterations, 
                              simulation_function(
                                N_per_subgroup = N_per_subgroup,
                                true_effect = c(# we can specify the true effect to be simulated here (in r difference) for each hypnothesis. This is expressed in difference in the correlation of the prejudice against the group in general and the males of the group, and the correlation of the prejudice against the group in general and the females of the group
                                  0.1, # Black prejudice correlation difference among non-blacks
                                  0, # Black prejudice correlation difference among black women
                                  0.1, # White prejudice correlation difference among non-whites
                                  0, # White prejudice correlation difference among white women
                                  0.1, # EastAsian prejudice correlation difference among non-east asians
                                  0, # EastAsian prejudice correlation difference among east asian women
                                  0.1, # Politician prejudice correlation difference among non-politicians
                                  0.1, # Criminal prejudice correlation difference among non-criminals
                                  0.1, # Police prejudice correlation difference among non-police
                                  0), # EastAsian prejudice correlation difference among US citizens who are non-east asians
                                # in the example data the effect sizes were: 0.55, 0, 0.45, 0, 0.1, 0, 0.3, 0.3, 0.5 respectively.
                                baseline_corr_matrix = t(matrix(c(1, 0.6, 0.6, # this is the correlation matrix that we will simulate if there is no correlation difference. The order of the variables in this matrix is: prejudice against the group in general, prejudice against the males of the group, prejudice against the females of the group 
                                                                  0.6, 1, 0.4,
                                                                  0.6, 0.4, 1), nrow = 3)),
                                sides_of_test = "one.sided",       # you can specify "one.sided" or "two.sided" here
                                minimal_effect_of_interest = 0.1, # expressed in r difference, just like the true effect above
                                CI_width = CI_width # width of the confidence interval for statistical inference. This is effectively the inverse of the p-value threshold for rejecting the null or the alternative hypothesis
                              ))


# power and inferential error rates
scenario_1 = table(simulation_output)/sum(table(simulation_output))



## Scenario 2. simulating no effect across the board
simulation_output = replicate(iterations, 
                              simulation_function(
                                N_per_subgroup = N_per_subgroup,
                                true_effect = c(# we can specify the true effect to be simulated here (in r difference) for each hypnothesis. This is expressed in difference in the correlation of the prejudice against the group in general and the males of the group, and the correlation of the prejudice against the group in general and the females of the group
                                  0, # Black prejudice correlation difference among non-blacks
                                  0, # Black prejudice correlation difference among black women
                                  0, # White prejudice correlation difference among non-whites
                                  0, # White prejudice correlation difference among white women
                                  0, # EastAsian prejudice correlation difference among non-east asians
                                  0, # EastAsian prejudice correlation difference among east asian women
                                  0, # Politician prejudice correlation difference among non-politicians
                                  0, # Criminal prejudice correlation difference among non-criminals
                                  0, # Police prejudice correlation difference among non-police
                                  0), # EastAsian prejudice correlation difference among US citizens who are non-east asians
                                # in the example data the effect sizes were: 0.55, 0, 0.45, 0, 0.1, 0, 0.3, 0.3, 0.5 respectively.
                                baseline_corr_matrix = t(matrix(c(1, 0.6, 0.6, # this is the correlation matrix that we will simulate if there is no correlation difference. The order of the variables in this matrix is: prejudice against the group in general, prejudice against the males of the group, prejudice against the females of the group 
                                                                  0.6, 1, 0.4,
                                                                  0.6, 0.4, 1), nrow = 3)),
                                sides_of_test = "one.sided",       # you can specify "one.sided" or "two.sided" here
                                minimal_effect_of_interest = 0.1, # expressed in r difference, just like the true effect above
                                CI_width = CI_width # width of the confidence interval for statistical inference. This is effectively the inverse of the p-value threshold for rejecting the null or the alternative hypothesis
                              ))


# power and inferential error rates
scenario_2 = table(simulation_output)/sum(table(simulation_output))


## Scenario 3. simulating correlation difference in all subgroups and prejudice types
simulation_output = replicate(iterations, 
                              simulation_function(
                                N_per_subgroup = N_per_subgroup,
                                true_effect = c(# we can specify the true effect to be simulated here (in r difference) for each hypnothesis. This is expressed in difference in the correlation of the prejudice against the group in general and the males of the group, and the correlation of the prejudice against the group in general and the females of the group
                                  0.1, # Black prejudice correlation difference among non-blacks
                                  0.1, # Black prejudice correlation difference among black women
                                  0.1, # White prejudice correlation difference among non-whites
                                  0.1, # White prejudice correlation difference among white women
                                  0.1, # EastAsian prejudice correlation difference among non-east asians
                                  0.1, # EastAsian prejudice correlation difference among east asian women
                                  0.1, # Politician prejudice correlation difference among non-politicians
                                  0.1, # Criminal prejudice correlation difference among non-criminals
                                  0.1, # Police prejudice correlation difference among non-police
                                  0.1), # EastAsian prejudice correlation difference among US citizens who are non-east asians
                                # in the example data the effect sizes were: 0.55, 0, 0.45, 0, 0.1, 0, 0.3, 0.3, 0.5 respectively.
                                baseline_corr_matrix = t(matrix(c(1, 0.6, 0.6, # this is the correlation matrix that we will simulate if there is no correlation difference. The order of the variables in this matrix is: prejudice against the group in general, prejudice against the males of the group, prejudice against the females of the group 
                                                                  0.6, 1, 0.4,
                                                                  0.6, 0.4, 1), nrow = 3)),
                                sides_of_test = "one.sided",       # you can specify "one.sided" or "two.sided" here
                                minimal_effect_of_interest = 0.1, # expressed in r difference, just like the true effect above
                                CI_width = CI_width # width of the confidence interval for statistical inference. This is effectively the inverse of the p-value threshold for rejecting the null or the alternative hypothesis
                              ))


# power and inferential error rates
scenario_3 = table(simulation_output)/sum(table(simulation_output))



## Scenario 4. simulating some very small effects (below the minimal effect of interest threshold)
simulation_output = replicate(iterations, 
                              simulation_function(
                                N_per_subgroup = N_per_subgroup,
                                true_effect = c(# we can specify the true effect to be simulated here (in r difference) for each hypnothesis. This is expressed in difference in the correlation of the prejudice against the group in general and the males of the group, and the correlation of the prejudice against the group in general and the females of the group
                                  0.03, # Black prejudice correlation difference among non-blacks
                                  0, # Black prejudice correlation difference among black women
                                  0.03, # White prejudice correlation difference among non-whites
                                  0, # White prejudice correlation difference among white women
                                  0.03, # EastAsian prejudice correlation difference among non-east asians
                                  0, # EastAsian prejudice correlation difference among east asian women
                                  0.03, # Politician prejudice correlation difference among non-politicians
                                  0.03, # Criminal prejudice correlation difference among non-criminals
                                  0.03, # Police prejudice correlation difference among non-police
                                  0), # EastAsian prejudice correlation difference among US citizens who are non-east asians
                                # in the example data the effect sizes were: 0.55, 0, 0.45, 0, 0.1, 0, 0.3, 0.3, 0.5 respectively.
                                baseline_corr_matrix = t(matrix(c(1, 0.6, 0.6, # this is the correlation matrix that we will simulate if there is no correlation difference. The order of the variables in this matrix is: prejudice against the group in general, prejudice against the males of the group, prejudice against the females of the group 
                                                                  0.6, 1, 0.4,
                                                                  0.6, 0.4, 1), nrow = 3)),
                                sides_of_test = "one.sided",       # you can specify "one.sided" or "two.sided" here
                                minimal_effect_of_interest = 0.1, # expressed in r difference, just like the true effect above
                                CI_width = CI_width # width of the confidence interval for statistical inference. This is effectively the inverse of the p-value threshold for rejecting the null or the alternative hypothesis
                              ))


# power and inferential error rates
scenario_4 = table(simulation_output)/sum(table(simulation_output))

final_results["scenario"] = c("expected pattern", "no effect", "effect everywhere", "small effects only")
final_results["all correct"] = c(scenario_1["correct inference"], scenario_2["correct inference"], scenario_3["correct inference"], scenario_4["correct inference"])
final_results["any false"] = c(scenario_1["false inference"], scenario_2["false inference"], scenario_3["false inference"], scenario_4["false inference"])
final_results["any inconclusive"] = c(scenario_1["inconclusive"], scenario_2["inconclusive"], scenario_3["inconclusive"], scenario_4["inconclusive"])

final_results
