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
        mutate("Characteristic" = explanatory[i]) %>%
        pivot_wider(names_from = !!sym(dependent), values_from = mean.sd) %>%
        relocate(Characteristic) %>%
        mutate(Levels = "mean.sd")
      
      
      z_median <- data %>%
        group_by(!!sym(dependent)) %>%
        summarise(median = spatstat.geom::weighted.median(!!sym(explanatory[i]), w = eave_weight),
                  q1 = spatstat.geom::weighted.quantile(!!sym(explanatory[i]), w = eave_weight, probs = 0.25),
                  q3 = spatstat.geom::weighted.quantile(!!sym(explanatory[i]), w = eave_weight, probs = 0.75)) %>%
        mutate("Characteristic" = explanatory[i]) %>%
        mutate(iqr = q3 -q1) %>%
        mutate(median.iqr = paste0(median, " (",iqr,")")) %>%
        select(-q1, -q3, -median, -iqr) %>%
        pivot_wider(names_from = !!sym(dependent), values_from = median.iqr) %>%
        relocate(Characteristic) %>%
        mutate(Levels = "median.iqr")
      
      # Combine!!
      summary_tbl_list[[i]] <- full_join(z_mean, z_median)
      
      
      # Else get sum of weights of each level
    } else if (length(unique(data %>% pull(!!sym(dependent) ) ) ) ==1) {
      
      # This is for when there is only one level in the dependent variable
      summary_tbl_list[[i]]   <- data %>%
        group_by(!!sym(explanatory[i])) %>%
        summarise(n = sum(eave_weight)) %>%
        ungroup() %>%
        mutate(perc = sprintf("%.1f",round(n/sum(n)*100,1))) %>%
        mutate_if(is.numeric, ~formatC(round(.,0), format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
        mutate(n_perc := paste0(n, " (", perc,"%)")) %>%
        select(-n, -perc) %>%
        rename("Levels"=explanatory[i], !!dependent := n_perc) %>%
        mutate("Characteristic" = explanatory[i]) %>%
        relocate(Characteristic)
      
    } else {

      summary_tbl_list[[i]] <- data %>%
        group_by(!!sym(explanatory[i]), !!sym(dependent)) %>%
        summarise(n = sum(eave_weight)) %>%
        ungroup() %>%
        group_by(!!sym(dependent)) %>%
        mutate(perc = sprintf("%.1f",round(n/sum(n)*100,1))) %>%
        mutate_if(is.numeric, ~formatC(round(.,0), format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
        mutate(n_perc := paste0(n, " (", perc,"%)")) %>%
        select(-n, -perc) %>%
        pivot_wider(names_from = !!sym(dependent), values_from = n_perc) %>%
        rename("Levels"=explanatory[i]) %>%
        mutate("Characteristic" = explanatory[i]) %>%
        relocate(Characteristic)
      
    }
  }
  
  # Combine list together to make dataset
  summary_tbl_wt <- summary_tbl_list %>%
    reduce(full_join)
  
  summary_tbl_wt <- mutate(summary_tbl_wt, 
                  Characteristic = ifelse(duplicated(Characteristic), '', Characteristic))
  
  summary_tbl_wt
}

##### Table of events and event rates by vaccine category
# Calculates the standardised mean differences (smd) between the uv and vacc for each of 
# the categorical explanatory variables for a vaccine type.

# Input:
# - data = the cohort descriptive dataset - z_chrt_desc

# Output:
# A table with weighted millions of person years spent with in each vaccination category in
# the cohort, and weighted count of events and event rates per million person years for
# hospitaliation, death, and hospitalisation or death post vaccination, and more than 14
# days post vaccination

# Table output to be used to plot comparisons between the matched and overall population 
# (before matching - crude)

event_summary_wt <- function(data){
  
  summary_tbl_list <- list()
  
  # First row is person years spent with each vaccination status in cohort
  first_row <-  t(select(data, starts_with('days'))) %*% pull(data, eave_weight)/(365.21 * 1000) 
  
  dependent <- grep('vacc_at', names(data), value = TRUE)
  
  for (i in 1:length(dependent)){
    summary_tbl_list[[i]] <- data %>%
      group_by(!!sym(dependent[i]) ) %>%
      summarise(n = sum(eave_weight)) %>%
      na.omit() %>%
      mutate(rate = sprintf('%.2f',n/first_row))  %>% 
      # format numbers with commas every 3 digits,  
      mutate_if(is.numeric, ~formatC(., format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
      mutate(n = paste0( n, ' (', rate, ')') )  %>%
      select(-rate) %>%
      pivot_wider(names_from = !!sym(dependent[i]), values_from = n) %>%
      mutate(Event = dependent[i])  
  }
  
  # Combine list together to make dataset
  summary_tbl_wt <- summary_tbl_list %>% 
    reduce(full_join) %>%
    mutate(Event = c('Hospitalisation',
                     '14 days prior to hospitalisation',
                     'Death',
                     '14 days prior to death',
                     'Hospitliation or death',
                     '14 days prior to hospitalisation or death')) 
  
  first_row <- formatC(sprintf('%.2f',first_row), format = "f", big.mark = ",", drop0trailing = TRUE)
  
  names(first_row) <- names(summary_tbl_wt)[1:5]
  
  first_row <- data.frame(as.list(first_row), stringsAsFactors = FALSE) %>% 
    mutate(Event = 'Person years (thousands)') %>%
    relocate(Event)
  
  summary_tbl_wt <- bind_rows(first_row, summary_tbl_wt) %>% relocate(uv, .after = Event)
  
  names(summary_tbl_wt) <- c('Event', 'Unvaccinated', 'First dose ChAdOx1', 
                             'Second dose ChAdOx1', 'First dose BNT162b2', 
                             'Second dose BNT162b2')
  
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
# of the means and the pooled sds, to create the standardised mean differences.

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


##### Weekly GLM relative risk ratios by a variable ####
# Calculates weekly relative risk of event between uv and vacc matches by a subset

# Input:
# - variable = variable to subset
# - level = level of interest to subset variable to

GLM_rr_var <- function(variable, level){
  
  z.fmla <- as.formula(paste("Surv(time_to_event, event)",
                             " ~ vacc + z.yr + vacc_type +",
                             paste(variable, collapse= "+")))
  
  z.agg <- pyears(z.fmla,
                  data=df_cc_ps_matches , scale=1, data.frame=TRUE)
  
  df_res <- z.agg$data
  df_res <- df_res %>% 
    mutate(pyears =round(pyears/365.25,1)) %>%
    mutate(Rate = round(event/pyears*1000)) %>% 
    #mutate(RR = Rate/first(Rate)) %>%
    mutate(id = paste0("z.yr", z.yr, ":vacc",vacc, ":vacc_type", vacc_type)) %>%
    select(-n) %>%
    # Subset to function 
    filter(!!sym(variable) == level)
  df_res
  
  z_pois <- z.agg$data %>%
    # Subset to function 
    filter(!!sym(variable) == level)
  
  
  # PB
  z_glm_pb <- glm(event ~ offset(log(pyears)) + -1 + z.yr + z.yr:vacc, 
                  family=poisson, data=z_pois, subset= vacc_type == "PB")
  summary(z_glm_pb)
  
  z_glm_pb_est <- data.frame(round(exp(cbind(z_glm_pb$coefficients, confint.default(z_glm_pb) ) ), 3)) %>%
    rownames_to_column(var ="id") %>%
    mutate(id = paste0(id, ":vacc_typePB")) %>%
    rename(est=2, lwr = 3, upr = 4) %>%   
    mutate(RR_est = paste0(round(est,2), " (", round(lwr,2), ", ", round(upr, 2), ")")) %>%   
    select(-c(est, lwr, upr))
  
  
  
  # AZ
  z_glm_az <- glm(event ~ offset(log(pyears)) + -1 + z.yr + z.yr:vacc, 
                  family=poisson, data=z_pois, subset= vacc_type == "AZ")
  summary(z_glm_az)
  
  z_glm_az_est <- data.frame(round(exp(cbind(z_glm_az$coefficients, confint.default(z_glm_az) ) ), 3)) %>%
    rownames_to_column(var ="id") %>%
    mutate(id = paste0(id, ":vacc_typeAZ")) %>%
    rename(est=2, lwr = 3, upr = 4) %>%   
    mutate(RR_est = paste0(round(est,2), " (", round(lwr,2), ", ", round(upr, 2), ")")) %>%  
    select(-c(est, lwr, upr))
  
  
  # Combine
  z_glm_vacc_output <- df_res %>%
    left_join(full_join(z_glm_pb_est,z_glm_az_est)) %>%
    rename(group = !!sym(variable)) %>%
    select(-id)
  
  z_glm_vacc_output
  
  
}


##### Weekly GAM relative risk ratios by vacc type and a variable ####
# Calculates GAM of daily relative risk of event between uv and vacc matches by a subset split
# by vaccine type

# Input:
# - z_vacc_type = vaccine type (AZ or PB)
# - variable = variable to subset
# - level = level of interest to subset variable to

GAM_rr_var <- function(z_vacc_type, variable, level){
  
  # Vaccc tpye title
  if(z_vacc_type == "PB"){
    z_vacc_title = "BNT162b2"
  } else {
    z_vacc_title = "ChAdOx1"
  }
  
  # Day
  z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,seq(0,max(df_cc_ps_matches$time_to_event),by=1) ))
  
  
  # Formula
  z.fmla <- as.formula(paste("Surv(time_to_event,event)",
                             " ~ vacc + z.yr + vacc_type +",
                             paste(variable, collapse= "+")))
  
  
  z.agg <- pyears(z.fmla,
                  data=df_cc_ps_matches , scale=1, data.frame=TRUE)
  
  
  
  z_pois <- z.agg$data %>%
    mutate(vacc = recode(vacc, "v1" = "uv", "v2" = "vacc")) %>%
    mutate(day =  as.numeric(z.yr)-1) %>%
    filter(vacc_type == z_vacc_type) %>%
    # Subset to function 
    filter(!!sym(variable) == level)
  
  myknots <-  seq(14,70, by =7)
  
  z_gam_vacc <- gam(event ~ offset(log(pyears)) + vacc + s(day, by=vacc, k=length(myknots)), 
                    knots = list(day = myknots), 
                    family=poisson, data=z_pois, subset= vacc_type == z_vacc_type)
  
  
  # Create cut off so length of vacc and uv match
  z <- z_pois %>%
    group_by(z.yr) %>%
    summarise(n=n()) %>%
    filter(n == 2)
  
  
  
  ## Plot GAM RR
  #for using for predictions
  z_nd <- z_pois %>%
    filter(z.yr %in% z$z.yr)
  
  #this is where the modelling starts
  z_p <- mgcv::predict.gam(z_gam_vacc, newdata=z_nd, type="lpmatrix")
  
  
  #columns
  c_uv <- grepl(":vaccuv", colnames(z_p)) #smoothcolumns for unvacc
  c_vacc <- grepl(":vaccvacc", colnames(z_p))  #smooth columns for vacc
  
  #rows
  r_uv <- with(z_nd, vacc=="uv")
  r_vacc <- with(z_nd, vacc=="vacc")
  
  #difference
  X <- z_p[r_vacc,] - z_p[r_uv,]
  
  #zero other columns
  X[,!(c_uv | c_vacc)] <- 0
  
  #set vaccine column to 1 to get the effect
  X[, "vaccvacc"] <- 1
  
  #offset difference
  offset_diff <- log(z_nd[r_uv,"pyears"]) - log(z_nd[r_vacc,"pyears"])
  
  
  
  difference <- X %*% coef(z_gam_vacc) + offset_diff
  
  se_difference <- sqrt(rowSums((X %*% vcov(z_gam_vacc, unconditional=TRUE)) * X))
  z_lcl <- difference - 1.96*se_difference
  z_ucl <- difference + 1.96*se_difference
  
  
  z_rr <- tibble(rr = exp(difference),
                 rr_lwr = exp(z_lcl),
                 rr_upr = exp(z_ucl),
                 day = 1:length(difference),
                 rr_ci = paste0(round(rr,2), " (", round(rr_lwr,2), ", ", round(rr_upr,2), ")"))
  
  z_rr <- z_rr %>%
    left_join(filter(z_pois, vacc == "vacc") %>%
                select(event, day) %>%
                rename(event_vacc =event))
  z_rr <- z_rr %>%
    left_join(filter(z_pois, vacc == "uv") %>%
                select(event, day) %>%
                rename(event_uv =event))
  
  
  ## P-value for quadratic term post 14 days
  z_quad <- glm(event ~ offset(log(pyears)) +  day*vacc + I(day^2)*vacc, 
                family=poisson, data=z_pois, subset = vacc_type == z_vacc_type & day > 14)
  
  zz_p <-predict.glm(z_quad, newdata=z_nd, type="response")
  z_nd$fit <- zz_p
  z_nd$prate <- z_nd$fit/z_nd$pyears
  
  z_quad_sum <- summary(z_quad)$coefficients
  upr <- summary(z_quad)$coefficients[6,1] + 1.96*summary(z_quad)$coefficients[6,2]
  lwr <- summary(z_quad)$coefficients[6,1] - 1.96*summary(z_quad)$coefficients[6,2]
  
  
  
  # Plot
  p_c <- z_rr %>%
    mutate(rr_upr = if_else(rr_upr>3, 3, rr_upr)) %>%
    ggplot() +
    geom_line(aes(x=day, y=rr), col=eave_blue) +
    geom_ribbon(aes(x=day, ymin=rr_lwr, ymax=rr_upr), alpha = 0.25, fill=eave_blue) +
    labs(x="Days to event", y="RR",
         subtitle = paste0(level),
         caption = paste0("Vaccinated quadratic coefficient after 14 days: ",round(z_quad_sum[6,1],4),
                          #" (", round(upr,3), ", ",round(lwr,3),"), 
                          ", p.value = ", round(z_quad_sum[6,4],3))) +
    geom_vline(xintercept = 14, linetype=2) +
    annotate("text", x=14.5, y=2, label = "14 days", hjust=0, size=3.5)+
    theme_light() +
    lims(y=c(0,3)) +
    geom_hline(yintercept = 1, linetype = 1) +
#    geom_hline(yintercept = 0.5, linetype = 3) +
#    annotate("text", x=80, y=0.6, label = "Waning threshold", size=3) +
    scale_x_continuous(breaks = seq(14,84, by = 7), 
                       limits = c(0,84))
  
  
  # Weekly z_rr
  z_rr_wkly <- z_rr %>%
    #mutate(rr_upr = ifelse(rr_upr==3, Inf, rr_upr)) %>%
    filter(day %in% seq(14,84, by = 14)) %>%
    mutate(rr_ci = paste0(round(rr,2), " (", round(rr_lwr,2), ", ", round(rr_upr,2), ")")) %>%
    mutate(y=1)
  
  
  p_c <- p_c +
    geom_text(aes(x=day, y=0, label = rr_ci), check_overlap = T, size=2, data=z_rr_wkly, inherit.aes = F)
  
  p_c
  
}



