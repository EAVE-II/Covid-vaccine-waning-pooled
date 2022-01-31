##########################################################
## Title: 1st dose COVID-19 vaccine waning
## Code author(s): Rachel Mulholland <rachel.mulholland@ed.ac.uk> 
## Description: 07_meta_analysis - Creates outputs for meta-analysis
##########################################################

#### 0 - Set up ####
# Run 0 - Set up in 03d_matching_modelling.R to set up df_cc_ps_matches to outcome

# Repeat for each outcome
z_event_endpoint <- "death_hosp"
#z_event_endpoint <- "hosp_covid"
#z_event_endpoint <- "death_covid"

# This controls the save path
#dose <- 'first_dose' 
dose <- 'second_dose'

#### 1 - Table of daily follow-up #####

z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,seq(0,max(df_cc_ps_matches$time_to_event),by=1) ))

## All vaccines
z.agg <- pyears(Surv(time_to_event,event) ~ vacc + z.yr,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)

daily_gam_all <- z.agg$data %>%
  mutate(z.yr =  as.numeric(z.yr)-1)

daily_gam_all <- daily_gam_all %>%
  select(-pyears) %>%
  pivot_wider(names_from = vacc, values_from = c("n", "event")) %>%
  mutate(vacc_type = "All") %>%
  select(6, 1, 2, 4, 3, 5)

head(daily_gam_all)


## Separate vaccines
z.agg <- pyears(Surv(time_to_event,event) ~ vacc + z.yr +vacc_type,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)

daily_gam_vacc <- z.agg$data %>%
  mutate(z.yr =  as.numeric(z.yr)-1)

daily_gam_vacc <- daily_gam_vacc %>%
  select(-pyears) %>%
  pivot_wider(names_from = vacc, values_from = c("n", "event")) %>%
  select(2, 1, 3, 5, 4, 6)



## Output
daily_gam_output <- bind_rows(daily_gam_all, daily_gam_vacc) 

write.csv(daily_gam_output, paste0("./data/", dose, "_", multiplicity_limit, "/meta-analysis/tbl_gam_", z_event_endpoint, ".csv"),
          row.names = F)



##### 2 - Table of daily follow-up by Age and Sex #####

df_cohort <- readRDS('./data/df_cohort.rds')

## Link age and sex into data
# Group age into 2 groups to allow for sufficient sample sizes
df_cc_ps_matches <- df_cc_ps_matches %>%
  #mutate(age_grp2 = ifelse(age_grp == "18-64", "18-64", "65+")) %>%
  left_join(select(df_cohort, EAVE_LINKNO, Sex))

# Function (if helpful)
GAM_tbl_var <- function(variable){
  
  # Day
  z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,seq(0,max(df_cc_ps_matches$time_to_event),by=1) ))
  
  
  # Formula
  z.fmla <- as.formula(paste("Surv(time_to_event,event)",
                             " ~ vacc + z.yr + vacc_type +",
                             paste(variable, collapse= "+")))
  
  
  z.agg <- pyears(z.fmla,
                  data=df_cc_ps_matches , scale=1, data.frame=TRUE)
  

  daily_gam_all <- z.agg$data %>%
    mutate(z.yr =  as.numeric(z.yr)-1) %>% 
    select(-pyears) %>%
    pivot_wider(names_from = vacc, values_from = c("n", "event")) %>%
    select(2, 3, 1, 4, 6, 5, 7 )
  
  
  daily_gam_all
  
}

# Sex
#daily_gam_output_sex <- GAM_tbl_var("Sex") %>%
 # rename(Subgroup = Sex)

# Age
daily_gam_output_age <- GAM_tbl_var("age_grp")%>%
  rename(Subgroup = age_grp)

# Stack ontop
#daily_gam_output_agesex <- bind_rows(daily_gam_output_sex, daily_gam_output_age)

# Output
write.csv(daily_gam_output_age, paste0("./data/", dose, '_', multiplicity_limit, "/meta-analysis/tbl_gam_age_", z_event_endpoint, ".csv"),
          row.names = F)


