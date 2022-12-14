##########################################################
## Title: 1st dose COVID-19 vaccine waning
## Code author(s): Steven Kerr <steven.kerr@ed.ac.uk>
## Description: Carries out meta-analysis of vaccine waning results
##########################################################


# Libraries
library(tidyverse)
library(mgcv)
library(lubridate)
library(survival)
library(survminer)
library(mgcv)
library(cowplot)

setwd('/conf/EAVE/GPanalysis/analyses/Covid-vaccine-waning-pooled')

# Colours 
eave_green <- rgb(54, 176, 136, maxColorValue = 255)
eave_blue <- rgb(71,93,167, maxColorValue = 255)
eave_blue2 <- rgb(0,192,209, maxColorValue = 255)
eave_gold <- rgb(255,192,0, maxColorValue = 255)
eave_orange <- rgb(244,143,32, maxColorValue = 255)

numeric_cols = c('n_uv', 'event_uv', 'n_vacc', 'event_vacc')

##### 1 Load data

# Scottish results
df_dose_1_scot <- read.csv('./data/first_dose_1/meta-analysis/tbl_gam_age_death_hosp.csv') 
df_dose_2_scot <- read.csv('./data/second_dose_1/meta-analysis/tbl_gam_age_death_hosp.csv') %>%
  rename(n_uv = n_v1, n_vacc = n_v2, event_uv = event_v1, event_vacc = event_v2)


# Norther Ireland results
df_dose_1_ni <- read.csv('./data/NI/section_6_gam_results_for_transfer_to_scotland_TRE_by_age_first_dose.csv') 
df_dose_2_ni <- read.csv('./data/NI/section_6_gam_results_for_transfer_to_scotland_TRE_by_age_second_dose.csv')

# Wales results
df_wales <- read.csv('./data/Wales/d_pdays_hosp_death_vacc_dose_age.csv')

df_dose_1_wales <- df_wales %>% filter(dose_cat =='First') %>% select(-dose_cat)
df_dose_2_wales <- df_wales %>% filter(dose_cat =='Second') %>% select(-dose_cat)

# England results
df_dose_1_eng <- read.csv('./data/England/daily_gam_output_age_dose1.csv') 
df_dose_2_eng <- read.csv('./data/England/daily_gam_output_age_dose2.csv')%>%
  rename(n_uv = n_v1, n_vacc = n_v2, event_uv = event_v1, event_vacc = event_v2)


##### 2 Functions

# Process NI data
process_ni <- function(df){
  
 df <- filter(df, endpoint_name == 'hospitalisation or death with COVID-19') %>%
        select(-endpoint_name) %>%
        rename(vacc_type = vaccine_type,
               event_uv = event_unvaccinated,
               event_vacc = event_vaccinated,
               n_uv = n_unvaccinated,
               n_vacc = n_vaccinated) %>%
        mutate( vacc_type = str_replace(vacc_type, 'az', 'AZ'),
                vacc_type = str_replace(vacc_type, 'pfizer', 'PB'),
                age_group = as.character(age_group)) %>%
        mutate( age_group = as.factor( str_replace(age_group, '<65', '18-64' ))) %>%
        select(vacc_type, age_group, day, numeric_cols) %>%
        arrange(vacc_type, age_group, day) %>%
        mutate_at(numeric_cols, ~replace(., is.na(.), 0))
}

# Process Scottish data
process_scot <- function(df){
  
  df <- rename(df, age_group = Subgroup,
               day = z.yr) %>%
        arrange(vacc_type, age_group, day)%>%
        mutate_at(numeric_cols, ~replace(., is.na(.), 0))
}

# Process Welsh data
process_wales <- function(df){
  
  df <- rename(df, vacc_type = case_vacc_name,
               age_group = age_cat,
               day = x_day,
               n_uv = control_pdays,
               event_uv = control_event,
               n_vacc = case_pdays,
               event_vacc = case_event) %>%
    mutate(age_group = as.character(age_group)) %>%
    mutate(age_group = as.factor( ifelse(age_group == ('80-110'), '80+', age_group))) %>%
    arrange(vacc_type, age_group, day)%>%
    mutate_at(numeric_cols, ~replace(., is.na(.), 0))
}

# Process English data
process_eng <- function(df){
  
  df <- select(df, -X) %>%
    rename(vacc_type = vacc_type_1,
           age_group = age_grp2,
           day = z.yr
           ) %>%
    mutate( vacc_type = str_replace(vacc_type, 'AstraZeneca', 'AZ'),
                            vacc_type = str_replace(vacc_type, 'Pfizer', 'PB'),
                            age_group = as.character(age_group)) %>%
    arrange(vacc_type, age_group, day)%>%
    mutate_at(numeric_cols, ~replace(., is.na(.), 0))
}


#Sum over age groups
sum_over_age_groups <- function(df){
  df %>% group_by(vacc_type, day) %>% 
    summarise(n_uv = sum(n_uv), event_uv = sum(event_uv),
              n_vacc = sum(n_vacc), event_vacc = sum(event_vacc)) %>%
    ungroup()  
}

# Get the maximum number of days that anyone survived in df by grouping variables ...
get_max_days <- function(df, ...){
  df %>% group_by(...) %>%
    summarise(max_day = max(day)) %>% 
    as.data.frame()
}

# Get minimum over countries of the maximum days that anyone survived by
# We have to censor at this time because we don't know if someone survived longer 
# than their last date in any country's data.
get_min_days <- function(max_days_df_list){
  
  df <- max_days_df_list[[1]]
  
  for (i in 2:length(max_days_df_list)){
    df$max_day <- pmin(df$max_day, pull(max_days_df_list[[i]], max_day))
  } 
  df
}


# This creates a dataframe that is censored as above by vacc_type and age group
truncate_days <-function(df, min_days){
  
  df %>% left_join(min_days) %>% 
    filter(day <= max_day) %>%
    select(-max_day)
}

# Make data longer with, with target trial group as a column
separate_groups <- function(df){
 
  rbind(select(df, -n_vacc, -event_vacc) %>% mutate(vacc = 'uv') %>%
          rename(n = n_uv, event = event_uv),  
        select(df, -n_uv, -event_uv) %>% mutate(vacc = 'vacc') %>%
          rename(n = n_vacc, event = event_vacc)) %>%
    mutate(vacc = as.factor(vacc)) 
}

###########################################################

# Pre-precessing
df_dose_1_ni = process_ni(df_dose_1_ni)
df_dose_2_ni = process_ni(df_dose_2_ni)

df_dose_1_scot = process_scot(df_dose_1_scot)
df_dose_2_scot = process_scot(df_dose_2_scot)

df_dose_1_wales = process_wales(df_dose_1_wales)
df_dose_2_wales = process_wales(df_dose_2_wales)

df_dose_1_eng = process_eng(df_dose_1_eng)
df_dose_2_eng = process_eng(df_dose_2_eng)


# List of dataframes 
df_list_dose_1 <- list(df_dose_1_ni, df_dose_1_scot, df_dose_1_wales, df_dose_1_eng)
df_list_dose_2 <- list(df_dose_2_ni, df_dose_2_scot, df_dose_2_wales, df_dose_2_eng)

# List of dataframes summed over age groups
df_list_dose_1_all <- lapply(df_list_dose_1, sum_over_age_groups)
df_list_dose_2_all <- lapply(df_list_dose_2, sum_over_age_groups)


# Get lists of minimum days, where we truncate
min_days_dose_1  <- lapply(df_list_dose_1, function(x) {get_max_days(x, vacc_type, age_group)}) %>% 
                    get_min_days()
  
min_days_dose_2 <- lapply(df_list_dose_2, function(x) {get_max_days(x, vacc_type, age_group)})%>%
                           get_min_days()

min_days_dose_1_all <- lapply(df_list_dose_1_all, function(x) {get_max_days(x, vacc_type)}) %>%
                               get_min_days()
min_days_dose_2_all <- lapply(df_list_dose_2_all, function(x) {get_max_days(x, vacc_type)})%>%
                               get_min_days()



# Truncate by days
df_list_dose_1 <- lapply(df_list_dose_1, function(x){ truncate_days(x, min_days_dose_1)})
df_list_dose_2 <- lapply(df_list_dose_2, function(x){ truncate_days(x, min_days_dose_2)})

df_list_dose_1_all <- lapply(df_list_dose_1_all, function(x){ truncate_days(x, min_days_dose_1_all)})
df_list_dose_2_all <- lapply(df_list_dose_2_all, function(x){ truncate_days(x, min_days_dose_2_all)})



# Add together numeric columns
df_dose_1 = df_list_dose_1[[1]]
df_dose_1_all = df_list_dose_1_all[[1]]

df_dose_2 = df_list_dose_2[[1]]
df_dose_2_all = df_list_dose_2_all[[1]]

for (i in 2:length(df_list_dose_1)){
  
  df_dose_1[, numeric_cols] = df_dose_1[, numeric_cols] + df_list_dose_1[[i]][, numeric_cols]
  
  df_dose_1_all[, numeric_cols] = df_dose_1_all[, numeric_cols] + df_list_dose_1_all[[i]][, numeric_cols]

  
  df_dose_2[, numeric_cols] = df_dose_2[, numeric_cols] + df_list_dose_2[[i]][, numeric_cols]
  
  df_dose_2_all[, numeric_cols] = df_dose_2_all[, numeric_cols] + df_list_dose_2_all[[i]][, numeric_cols]    
  
}





# Separate out the two target trial arms
df_dose_1 <- separate_groups(df_dose_1)
df_dose_2 <- separate_groups(df_dose_2)

df_dose_1_all <- separate_groups(df_dose_1_all)
df_dose_2_all  <- separate_groups(df_dose_2_all)


#### 3 GAMs ####
# Final figures for publication by vacc type



# ### C: GAM
# # Daily z.yr
# z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,seq(0,max(df_cc_ps_matches$time_to_event),by=1) ))
# 
# # Person years
# z.agg <- pyears(Surv(time_to_event,event) ~ vacc + z.yr +vacc_type,
#                 data=df_cc_ps_matches , scale=1, data.frame=TRUE)
# 
# # Calculate the number of events post 70 days (cut off due to low numbers)
# z.agg$data %>% mutate(long = as.numeric(z.yr) >70) %>% 
#   group_by(vacc, vacc_type, long) %>% summarise(N=sum(event))
# 
# # Take data and subset to vaccine type
# z_pois <- z.agg$data %>%
#   # Recode v1 -> uv and v2 -> vacc for convenicene, so same functions can be used
#   mutate(vacc = recode(vacc, "v1" = "uv", "v2" = "vacc")) %>%
#   dplyr::mutate(day =  as.numeric(z.yr)-1) %>%
#   dplyr::filter(vacc_type == z_vacc_type) 
# 
# # There are some people whose match had an event earlier than them, and there is
# # no follow up time for the other group. Add entries for those that are missing.
# grid <- expand.grid( min(z_pois$day):max(z_pois$day), c('uv', 'vacc'))
# 
# names(grid) <- c('day', 'vacc')
# 
# z_pois <- full_join(z_pois, grid)
# 
# z_pois <- mutate_at(z_pois,  c('n', 'event', 'pyears'), ~replace(., is.na(.), 0))

plot_gam <- function(df, vacc, time_limit){

  # Assign vaccine title
  if(vacc == "PB"){
    z_vacc_title = "BNT162b2"
  } else {
    z_vacc_title = "ChAdOx1"
  }
  
  myknots <-  seq(14,time_limit-7, by =7)
  
  z_gam_vacc <- gam(event ~ offset(log(n)) + vacc + s(day, by=vacc, k=length(myknots)), 
                    knots = list(day = myknots), 
                    family=poisson, data=df)
  
  ## Calculate GAM RR
  # Uses methods from https://fromthebottomoftheheap.net/2017/10/10/difference-splines-i/
  
  z_nd <- df #for using for predictions
  
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
  
  
  difference <- X %*% coef(z_gam_vacc) 
  
  # Obtain 95% CIs
  se_difference <- sqrt(diag(X %*% vcov(z_gam_vacc, unconditional=TRUE) %*% t(X)))
  z_lcl <- difference - 1.96*se_difference
  z_ucl <- difference + 1.96*se_difference
  
  # Create dataset with RRs and 95% CIs for plotting
  z_rr <- tibble(rr = exp(difference),
                 rr_lwr = exp(z_lcl),
                 rr_upr = exp(z_ucl),
                 day = 1:length(difference),
                 rr_ci = paste0(round(rr,2), " (", round(rr_lwr,2), ", ", round(rr_upr,2), ")"))
  
  # Join on no. events for vacc and uv
  z_rr <- z_rr %>%
    left_join(filter(df, vacc == "vacc") %>%
                select(event, day) %>%
                rename(event_vacc =event))
  z_rr <- z_rr %>%
    left_join(filter(df, vacc == "uv") %>%
                select(event, day) %>%
                rename(event_uv =event))
  
  
  
  # Save table
  #write.csv(z_rr, paste0("./output/second_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/poisson_gam_", z_vacc_type, ".csv"))
  
  ## Waning p-value
  # p-value for quadratic term post 14 days
  z_quad <- glm(event ~ offset(log(n)) +  day*vacc + I(day^2)*vacc,
                family=poisson, data=df, subset =  day > 14)
  z_linear <- glm(event ~ offset(log(n)) +  day*vacc,
                  family=poisson, data=df, subset = day > 14)
  
  
  summary(z_quad)
  # Check
  zz_p <-predict.glm(z_quad, newdata=z_nd, type="response")
  z_nd$fit <- zz_p
  z_nd$prate <- z_nd$fit/z_nd$n
  z_nd %>% ggplot(aes(x=day, y=prate, colour=vacc)) + geom_line()
  
  # Get estimates and 95% CIs
  z_quad_sum <- summary(z_quad)$coefficients
  upr <- summary(z_quad)$coefficients[6,1] + 1.96*summary(z_quad)$coefficients[6,2]
  lwr <- summary(z_quad)$coefficients[6,1] - 1.96*summary(z_quad)$coefficients[6,2]
  
  z_title <- "COVID-19 hospitalisations or deaths"
  
  ## Plot RRs
  p_c <- z_rr %>%
    mutate(rr_upr = if_else(rr_upr>3, 3, rr_upr)) %>%
    ggplot() +
    geom_line(aes(x=day, y=rr), col=eave_blue) +
    geom_ribbon(aes(x=day, ymin=rr_lwr, ymax=rr_upr), alpha = 0.25, fill=eave_blue) +
    labs(x="Days to event", y="RR",
         subtitle = paste0("Rate ratios of ", z_title," for ", z_vacc_title),
         caption = paste0("Vaccinated quadratic coefficient after 14 days: ",round(z_quad_sum[6,1],4),
                          #" (", round(upr,3), ", ",round(lwr,3),"), 
                          ", p.value = ", round(z_quad_sum[6,4],3))) +
    geom_vline(xintercept = 14, linetype=2) +
    annotate("text", x=14.5, y=2, label = "14 days", hjust=0, size=3.5)+
    theme_light() +
    lims(y=c(0,3)) +
    geom_hline(yintercept = 1, linetype = 1) +
    #  geom_hline(yintercept = 0.5, linetype = 3) +
    #  annotate("text", x=80, y=0.6, label = "Waning threshold", size=3) +
    scale_x_continuous(breaks = seq(14,time_limit, by = 7), 
                       limits = c(0,time_limit+7))
  
  
  ## Weekly RRs from z_rr to overlap ontop of p_c
  z_rr_wkly <- z_rr %>%
    #mutate(rr_upr = ifelse(rr_upr==3, Inf, rr_upr)) %>%
    filter(day %in% seq(14,time_limit, by = 14)) %>%
    mutate(rr_ci = paste0(round(rr,2), " (", round(rr_lwr,2), ", ", round(rr_upr,2), ")")) %>%
    mutate(y=1)
  
  # Combine p_c and z_rr_wkly
  p_c <- p_c +
    geom_text(aes(x=day, y=0, label = rr_ci), check_overlap = T, size=3, data=z_rr_wkly, inherit.aes = F)

  p_c
}

time_limit = 98

save_plots <- function(df, dose){
  
  for (row in 1:nrow(min_days_dose_1)){
    
    vaccine = min_days_dose_1$vacc_type[row]
    age = min_days_dose_1$age_group[row]
    
    print(vaccine)
    print(age)
    
    plot <- plot_gam( filter(df, vacc_type == vaccine, age_group == age), vaccine, time_limit)
    
    ggsave(paste0('./output/gam_dose_', dose, '_', vaccine, '_', age, '.png'), plot)
    
  }
  
}


# By age group
save_plots(df_dose_1, 1)
save_plots(df_dose_2, 2)


# All ages
plot_gam(df_dose_1_all %>% filter(vacc_type == 'AZ'), 'AZ', time_limit)
ggsave(paste0('./output/gam_dose_1_AZ_all.png'))

plot_gam(df_dose_1_all %>% filter(vacc_type == 'PB'), 'PB', time_limit)
ggsave(paste0('./output/gam_dose_1_PB_all.png'))

plot_gam(df_dose_2_all %>% filter(vacc_type == 'AZ'),'AZ', time_limit)
ggsave(paste0('./output/gam_dose_2_AZ_all.png'))

plot_gam(df_dose_2_all %>% filter(vacc_type == 'PB'), 'PB', time_limit)
ggsave(paste0('./output/gam_dose_2_PB_all.png'))

