##########################################################
## Title: 1st dose COVID-19 vaccine waning
## Code author(s): Rachel Mulholland <rachel.mulholland@ed.ac.uk> 
##                 Chris Robertson <chrisobertson@nhs.net>
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

