##########################################################
## Title: 1st dose COVID-19 vaccine waning
## Code author(s): Rachel Mulholland <rachel.mulholland@ed.ac.uk> 
##                 Chris Robertson <chrisrobertson@nhs.net>
## Description: 00_functions - Unique functions for analysis
##########################################################

#### Libraries ####
library("spatstat")


#### Summary table using weights ####
# Creates a table of cohort summaries using weights
# For categorical variables, the sum of the weights and the % is calculated
# For numerical variables, the weighed mean and weighted sd are calculated, as well as
# the weighted median and weighted IQR

# Input:
# - data = the dataset (must have weights named eave_weight)
# - dependent = a character of the dependent variables name
# - explanatory = a string of characters of the explanatory variables

# Output:
# A table with each explanatory variable as a row (multiple rows for each category if categorical)
# with two columns of the weighted summaries for the levels in the dependent variable

summary_factorlist_wt <- function(data, dependent, explanatory){
  # Create list to put in summaries into each element
  summary_tbl_list <- list()
  
  for(i in 1:length(explanatory)){
    
    # Extract variable
    n <- data %>%
      pull(!!sym(explanatory[i]))
    
    # If numeric then make weighted mean
    if(is.numeric(n)) {
      z_mean <- data %>%
        group_by(!!sym(dependent)) %>%
        summarise(mean = round(weighted.mean(!!sym(explanatory[i]), w = eave_weight),1),
                  sd = round(sqrt(spatstat.geom::weighted.var(!!sym(explanatory[i]), w = eave_weight)),1)) %>%
        mutate(mean.sd = paste0(mean, " (",sd,")")) %>%
        select(-mean, -sd) %>%
        mutate("characteristic" = explanatory[i]) %>%
        pivot_wider(names_from = !!sym(dependent), values_from = mean.sd) %>%
        relocate(characteristic) %>%
        mutate(levels = "mean.sd")
      
      
      z_median <- data %>%
        group_by(!!sym(dependent)) %>%
        summarise(median = spatstat.geom::weighted.median(!!sym(explanatory[i]), w = eave_weight),
                  q1 = spatstat.geom::weighted.quantile(!!sym(explanatory[i]), w = eave_weight, probs = 0.25),
                  q3 = spatstat.geom::weighted.quantile(!!sym(explanatory[i]), w = eave_weight, probs = 0.75)) %>%
        mutate("characteristic" = explanatory[i]) %>%
        mutate(iqr = q3 -q1) %>%
        mutate(median.iqr = paste0(median, " (",iqr,")")) %>%
        select(-q1, -q3, -median, -iqr) %>%
        pivot_wider(names_from = !!sym(dependent), values_from = median.iqr) %>%
        relocate(characteristic) %>%
        mutate(levels = "median.iqr")
      
      # Combine!!
      summary_tbl_list[[i]] <- full_join(z_mean, z_median)
      
      
      # Else get sum of weights of each level
    } else {
      
      summary_tbl_list[[i]] <- data %>%
        group_by(!!sym(explanatory[i]), !!sym(dependent)) %>%
        summarise(n = sum(eave_weight)) %>%
        ungroup() %>%
        group_by(!!sym(dependent)) %>%
        mutate(perc = sprintf("%.1f",round(n/sum(n)*100,1))) %>%
        mutate(n_perc = paste0(round(n,0), " (", perc,"%)")) %>%
        select(-n, -perc) %>%
        pivot_wider(names_from = !!sym(dependent), values_from = n_perc) %>%
        rename("levels"=explanatory[i]) %>%
        mutate("characteristic" = explanatory[i]) %>%
        relocate(characteristic)
      
    }
    
    
    
  }
  
  # Combine list together to make dataset
  summary_tbl_wt <- summary_tbl_list %>%
    reduce(full_join)
  
  summary_tbl_wt
}




##### Covariate balance #####
# Calculates the standardised mean differences (smd) between the uv and vacc for each of 
# the categorical explanatory variables for a vaccine type.

# Input:
# - data = the dataset - either the matched cohort (df_cc_desc) or entire cohort (z_chrt_desc)
# - explanatory = a string of characters of the explanatory variables with labels as their names
# - z_vacc_type = vaccine type (PB or AZ)

# Output:
# A table with each explanatory variable and their categories as a row with
# the weighted means and weighted variances. Differences between the two vaccination groups
# (vacc and uv) of the means and the pooled sds, to create the standardised mean differences.

# Table output to be used to plot comparisons between the matched and overall population 
# (before matching - crude)


cov_bal_vacc_fn <- function(data, explanatory, z_vacc_type){
  
  covariate_balance <- list()
  
  ## Find whether it is the matched or original dataset
  data_tbl <- table(data$vacc, data$vacc_type)
  
  # If summaries for uv by vacc type = 0 then this refers to the original dataset
  if(data_tbl[1,1]==0 &data_tbl[1,2]==0){
    
    # Replace the vacc_type for the uv as the z_vacc_type
    data <- data %>%
      filter(vacc_type == z_vacc_type | vacc == "uv") %>%
      mutate(vacc_type = z_vacc_type)
    
    
  } else {
    data <- data %>%
      filter(vacc_type == z_vacc_type)
  }
  
  
  # For each explanatory variable, get the standardised mean differences
  for(i in 1:length(explanatory)){
    
    # Identify whether binary or not
    n <- data %>%
      pull(!!sym(explanatory[i]))
    
    # If binary then choose a reference level and use this to compare between the other group
    if(length(unique(n)) == 2){
      
      level1 <- max(n)
      
      
      covariate_balance[[i]] <- data %>%
        mutate(var = ifelse(!!sym(as.character(explanatory)[i]) == level1, 1, 0)) %>%
        group_by(vacc) %>%
        summarise(weighted.mean = weighted.mean(var, w= eave_weight, na.rm = T),
                  var = spatstat.geom::weighted.var(var, w = eave_weight)) %>%
        mutate(mean.diff = weighted.mean[vacc=="uv"] -  weighted.mean[vacc=="vacc"],
               sd.pooled = sqrt((var[vacc=="uv"] + var[vacc=="vacc"])/2),
               smd = mean.diff/sd.pooled) %>%
        #select(smd) %>%
        mutate(characteristic = names(explanatory[i])) %>%
        mutate(levels = level1) %>%
        distinct()
      
      
    } else
      # If more than one group, take each level at a time and use as binary variable
      if(length(unique(n)) > 2){
        
        n_levels <- unique(n)
        
        
        covariate_balance_multi <- list()
        
        for(j in 1:length(n_levels)){
          level1 <- n_levels[j]
          
          if(!is.na(level1)) {
            covariate_balance_multi[[j]] <- data %>%
              mutate(var = ifelse(!!sym(explanatory[i]) == level1, 1, 0)) %>%
              group_by(vacc) %>%
              summarise(weighted.mean = weighted.mean(var, w= eave_weight, na.rm=T),
                        var = spatstat.geom::weighted.var(var, w = eave_weight, na.rm=T)) %>%
              mutate(mean.diff = weighted.mean[vacc=="uv"] -  weighted.mean[vacc=="vacc"],
                     sd.pooled = sqrt((var[vacc=="uv"] + var[vacc=="vacc"])/2),
                     smd = mean.diff/sd.pooled) %>%
              #select(smd) %>%
              mutate(characteristic = names(explanatory[i])) %>%
              mutate(levels = level1) %>%
              distinct() 
          } else {
            
            covariate_balance_multi[[j]] <- data %>%
              mutate(var = replace_na(!!sym(explanatory[i]), 1)) %>%
              mutate(var = ifelse(var == 1, 1, 0)) %>%
              group_by(vacc) %>%
              summarise(weighted.mean = weighted.mean(var, w= eave_weight, na.rm=T),
                        var = spatstat.geom::weighted.var(var, w = eave_weight, na.rm=T)) %>%
              mutate(mean.diff = weighted.mean[vacc=="uv"] -  weighted.mean[vacc=="vacc"],
                     sd.pooled = sqrt((var[vacc=="uv"] + var[vacc=="vacc"])/2),
                     smd = mean.diff/sd.pooled) %>%
              #select(smd) %>%
              mutate(characteristic = names(explanatory[i])) %>%
              mutate(levels = level1) %>%
              distinct() 
            
          }
          
          
          
        }
        
        covariate_balance [[i]] <- covariate_balance_multi %>%
          reduce(full_join)
        
        
      }
    
  }
  
  covariate_balance_all <- covariate_balance %>%
    reduce(full_join) %>%
    mutate(label = paste0(characteristic, ": ", levels)) %>%
    relocate(label) %>%
    arrange(desc(characteristic)) %>%
    mutate(vacc_type = z_vacc_type)
  
  covariate_balance_all
  
  
}


