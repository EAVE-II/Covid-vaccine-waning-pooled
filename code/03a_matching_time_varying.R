##########################################################
## Title: 1st dose COVID-19 vaccine waning
## Code author(s): Rachel Mulholland <rachel.mulholland@ed.ac.uk> 
##                 Chris Robertson <chrisrobertson@nhs.net>
## Description: 03a_matching_time_varying - Performs time-varying
##              propensity score matching
##########################################################

# Colours 
eave_green <- rgb(54, 176, 136, maxColorValue = 255)
eave_blue <- rgb(71,93,167, maxColorValue = 255)
eave_blue2 <- rgb(0,192,209, maxColorValue = 255)
eave_gold <- rgb(255,192,0, maxColorValue = 255)
eave_orange <- rgb(244,143,32, maxColorValue = 255)

#### 0 - Set up ####

# Upper limit on the number of times anyone is used as a match
# If you don't want an upper limit, set equal to ''.
multiplicity_limits <- c(1, 5, 10)

# Load in df_cohort and df_vaccinations data 
df_cohort <- readRDS("./data/df_cohort.rds")
df_vaccinations <- readRDS("./data/df_vaccinations.rds")

## Options list (for conditions)
output_list <- list()
# Exclude previous positive testing?
output_list$prev_pos <- "keep"
# Exclude care homes?
output_list$carehomes <- "remove"

# What event end point file to look at?
z_event_endpoint <- "death_hosp"
#z_event_endpoint <- "death" 
#z_event_endpoint <- "hosp_covid"
#z_event_endpoint <- "positive_test"

# Assign z_event to outcome of interest and assign z_test_event to number of days to
# add to specimen data if event time is missing
if (z_event_endpoint =="hosp_covid") {z_event <- covid_hospitalisations
z_test_event <- 7
z_title <- "COVID-19 hospitalisations"}
if (z_event_endpoint =="icu_death") {z_event <- covid_icu_death
z_test_event <- 14}
if (z_event_endpoint =="death") {z_event <- covid_death
z_test_event <- 21
z_title <- "COVID-19 deaths"}
if (z_event_endpoint =="death_hosp") {z_event <- covid_hosp_death
z_test_event_hosp <- 7
z_test_event_death <- 21
z_title <- "COVID-19 hospitalisations and deaths"}
if (z_event_endpoint =="positive_test") {z_event <- positive_test
z_title <- "COVID-19 positive test"}
colnames(z_event)

# Find end date according to admission date
a_end <- as.Date("2021-06-30")

# Filter event data to end date
z_event <- z_event %>%
  filter(admission_date <= a_end)

## NOTE: admission_date is generic - means hospital admission date and/or death date


##### 1 - Prepare exposed cohort data #####

## Vaccination data
# Make anyone vaccinated after the maximum endpoint time unvaccinated
z_vaccinations <- filter(df_vaccinations, 
                         date_vacc_1 <= a_end) %>% 
  mutate(vacc_type_2 = if_else(date_vacc_2 > a_end, NA_character_ , vacc_type_2),
         date_vacc_2 = ifelse(date_vacc_2 > a_end, as.Date(NA), date_vacc_2) )



## Whole cohort
# Group ages over 80s in whole cohort too (for linking/matching)
# Replace EAVE II link with _uv suffix
z_chrt <- df_cohort %>%
  dplyr::rename(EAVE_LINKNO_uv = EAVE_LINKNO) %>%
  filter(!is.na(ageYear)) %>%
  mutate(ur6_2016_name = replace_na(ur6_2016_name, "Unknown")) 

# Add in vacc variable
z_chrt <- z_chrt %>%
  # Create vacc flag if EAVE ID = Vacc EAVE ID
  left_join(select(z_vaccinations, EAVE_LINKNO, vacc_type, vacc_type_2,
                   date_vacc_1, flag_incon),
            by=c("EAVE_LINKNO_uv" = "EAVE_LINKNO")) %>%
  # Only keep people with consistent vaccination records
  filter( flag_incon == 0)  %>%
  mutate(vacc = case_when( !is.na(date_vacc_1) ~ 1,
                            TRUE ~ 0)) %>%
  # Our cohort will be uv, AZ or PB only. Remove moderna
  filter(vacc_type %in% c('AZ', 'PB', 'UNK') | is.na(vacc_type) ) %>%
  filter(vacc_type_2 %in% c('AZ', 'PB', 'UNK') | is.na(vacc_type_2) ) 


# Add age as individual year (truncated)
z_chrt <- z_chrt %>% mutate(ageYear = if_else(ageYear >= 100,100, ageYear)) %>% 
  mutate(ageYear = if_else(ageYear >= 95 & ageYear <= 99,97, ageYear)) %>% 
  mutate(ageYear = if_else(ageYear >= 90 & ageYear <= 94,92, ageYear)) %>% 
  mutate(ageYear = if_else(ageYear >= 80 & ageYear <= 89, trunc(ageYear/2)*2, ageYear)) 

# Add in household information
Cohort_Household <- readRDS("/conf/EAVE/GPanalysis/outputs/temp/Cohort_Household.rds") %>%
  mutate(n_hh_gp = cut(n_hh, breaks=c(0,1,2,5,10,30,100,max(n_hh)),
                       labels=c("1", "2", "3-5", "6-10", "11-30", "31-100", "101+")))%>% 
  mutate(ave_hh_age=if_else(is.na(ave_hh_age), mean(ave_hh_age, na.rm=T), ave_hh_age) )

z_chrt <- z_chrt %>%
  left_join(select(Cohort_Household, EAVE_LINKNO,
                   n_hh_gp, ave_hh_age, care_home_elderly),
            by= c("EAVE_LINKNO_uv"= "EAVE_LINKNO"))

# Eliminate carehomes (if required)
output_list$carehomes == "remove"
if (output_list$carehomes == "remove") {
  z_chrt <- z_chrt %>% 
    filter(care_home_elderly == 0) %>% 
    dplyr::select(-care_home_elderly)
}

# Total no. vacc
length(which(z_chrt$vacc==1))



#### 2 - Perform time-varying propensity scores ####
# Performs monthly propensity score matching models and excludes those not eligible for matching

# !! WARNING !! - Run time for code takes a long time depending on size of cohort

# Process: Nested loop:
# - Outer loop (j) = Monthly time-periods to fit time-varying propensity score
# - Inner loop (i) = Matching process split into 5 batches (see matching process line 138).
#                   This helps overcome issues with memory in the environment

# Propensity score matching model: age_gp*Sex + simd2020_sc_quintile + n_tests_gp + n_risk_gps +
# Council +ur6_2016_name + ave_hh_age + n_hh_gp + EAVE_Smoke + bmi_cat + in_hosp_status + prev_positive_status
# With in_hosp_status + prev_positive_status as time-varying covariates

# Matching process: Exact matching of 1% bands, individual age and Council Area

# Ineligible matching conditions: if controls they were vaccinated before the exposed, 
# if they had a COVID-19 hospitalisation or any death before the vaccination date of the exposed.


### Set up - Monthly time period
# Create loop for month
a_months <- seq(a_begin,a_end,by='months')
a_months[length(a_months)+1] <- a_end # Replace end date with end date (as sequence may miss this)
a_months[2] <- as.Date("2021-01-03") # Replace 2nd months start date with the 1st day of AZ

# Create list to assign monthly datasets to
z_merge_month <- list()


### Create a loop for the jth month and perform matching
# !! WARNING !!: Long run time
# Prints of statuses within the loop as been made to keep track of progress

for(j in 1:(length(a_months)-1)){
  
  # Extract start and end date within jth monthly time period
  z_strt <- a_months[j]
  z_end <- a_months[j+1]
  
  # Print to give progress on what time-period loop is on
  print(paste0("Time-period: ", j, " (", z_strt, " to ", z_end, ") - Started (0%)"))
   
  ## Calculate in hospital status
  # In hospital status = 1 if within 4 weeks of start date
  any_hosp_month <- all_hospitalisations %>%
    left_join(z_vaccinations) %>%
    mutate(intervals = interval(start = admission_date, end = discharge_date)) %>%
    mutate(in_hosp_status = int_overlaps(intervals, interval(start = z_strt-28, end =z_strt))) %>%
    select(EAVE_LINKNO, in_hosp_status) %>%
    group_by(EAVE_LINKNO) %>%
    summarise(in_hosp_status = sum(in_hosp_status)) %>%
    mutate(in_hosp_status = if_else(in_hosp_status == 0 | is.na(in_hosp_status),0,1))
  
  # Link in hospital status
  z_chrt_j <- z_chrt %>%
    left_join(any_hosp_month, by= c("EAVE_LINKNO_uv"= "EAVE_LINKNO")) %>%
    mutate(in_hosp_status = replace_na(in_hosp_status, 0 )) %>%
    mutate(in_hosp_status = as.character(in_hosp_status))
  
  
  ## Previous positive test
  # Previously testing positive status = 1 if tested positive before period
  z_chrt_j <- z_chrt_j %>%
    #mutate(prev_positive_status = if_else(SpecimenDate < date_vacc_1, 1, 0)) %>%
    mutate(prev_positive_status = if_else(SpecimenDate < z_strt, 1, 0))%>%
    mutate(prev_positive_status = replace_na(prev_positive_status, 0))%>%
    mutate(prev_positive_status = as.character(prev_positive_status))
  
  
  ## Adjust cohort to time-period
  # Only includes those who have not already been vaccinated, 
  # or have been vaccinated after the start date
  z_chrt_j <- z_chrt_j %>%
    # Filter anyone who has already been vaccinated
    filter(date_vacc_1 > z_strt | is.na(date_vacc_1)) %>%
    # If they are vaccinated after the time period, replace by 0
    mutate(vacc = if_else(date_vacc_1 > z_end | is.na(date_vacc_1), 0, vacc))%>%
    # If they were vaccinated after the time period, replace the date of vaccination by NA
    mutate(date_vacc_1 = ifelse(date_vacc_1 > z_end | is.na(date_vacc_1), NA, date_vacc_1)) %>%
    mutate(date_vacc_1 = as.Date(date_vacc_1, origin=as.Date("1970-01-01")))
  

  
  ## Propensity score 
  # Model 2: Age, sex, SIMD, n tests, n risk groups, council area, UR, 
  #average household age, number of people in household, EAVE II smoking status,
  # BMI, in-hospital status and previous positive status
  m_j <- glm(vacc ~ age_gp*Sex + simd2020_sc_quintile + n_tests_gp + n_risk_gps +
               Council +ur6_2016_name + ave_hh_age + n_hh_gp + EAVE_Smoke + bmi_cat +
               in_hosp_status + prev_positive_status,
             data=z_chrt_j, family=binomial)
  
  # Print progress
  print(paste0("Time-period: ", j, " - Propensity score model complete (20%)"))
  
  # Predict onto dataset and group into 1% bands
  z_ps <- predict(m_j, newdata=z_chrt_j, type="response")
  z_chrt_j <- z_chrt_j %>% 
    mutate(prop_score = z_ps) %>%
    mutate(ps_grp = cut(prop_score, breaks=seq(0,1, by=0.01)))
  
  
  ### Extract vaccination population only
  z_df <- z_chrt_j %>%
    filter(!is.na(date_vacc_1)) %>%
    dplyr::rename(EAVE_LINKNO_vacc = EAVE_LINKNO_uv)
  

  ### Prepare vaccination population and entire population for matching
  ## Only keep the info needed (to quicken up process) - prop score (grouped), Council area and age
  # Vaccinated population
  z_df <- z_df %>%
    select(EAVE_LINKNO_vacc, prop_score, ps_grp, Council, date_vacc_1, ageYear) %>%
    # Link event data
    left_join(select(z_event, EAVE_LINKNO, admission_date), by = c("EAVE_LINKNO_vacc"= "EAVE_LINKNO")) %>%
    rename(date_vacc_1_vacc = date_vacc_1, prop_score_vacc = prop_score, admission_date_vacc = admission_date)
  
  head(z_df)
  
  # Cohort population - also link vaccine data, event data and death data 
  z_chrt_j <- z_chrt_j %>%
    select(EAVE_LINKNO_uv,prop_score, ps_grp, Council, date_vacc_1, ageYear) %>%
    rename(date_vacc_1_uv = date_vacc_1, prop_score_uv=prop_score) %>%
    # Link event data
    left_join(z_event %>%
                select(-SpecimenDate) %>%
                rename(admission_date_uv = admission_date), 
              by=c("EAVE_LINKNO_uv" = "EAVE_LINKNO")) %>%
    # Link death data
    left_join(any_death %>%
                select(-SpecimenDate) %>%
                rename(any_death_date_uv = admission_date), by=c("EAVE_LINKNO_uv" = "EAVE_LINKNO"))
  
  head(z_chrt_j)
  
  # Print progress
  print(paste0("Time-period: ", j, " - Matching process beginning (30%)"))
  
  
  ### Perform the matching
  ## Exact match through propensity scores (grouped), individual age in years and Council area
  
  ## Split cohort into batches of 5
  # This is to help with memory as sample sizes can be large and R session may crash
  # Split number of rows into 5 and use in a loop
  nrows <- 1:nrow(z_df)
  loop_breaks5 <- split(nrows, cut(seq_along(nrows),
                                   5,
                                   labels = FALSE))
  
  # Create list to assign data through loop
  z_merge_ps_list <- list()
  
  
  # Begin loop
  for(i in 1:5){
    
    # Splitting the data 
    z_df_i <- z_df[loop_breaks5[[i]], ]
    
    # Merge the vaccination cohort with the whole cohort by the matching variables
    z_merge_ps_i <- z_df_i %>% 
      left_join(z_chrt_j, by=c("ps_grp","Council", "ageYear"))
    
    # Since merging with characteristics, this will assign all study IDs matching on these 
    # characteristics per vaccinated person in z_df (i.e. multiple IDs per vaccinated person)
    
    ## Vaccinated control
    # Remove any matches where data of vaccination of control is before vaccination case
    z_merge_ps_i <- z_merge_ps_i %>% 
      filter(date_vacc_1_vacc < date_vacc_1_uv | is.na(date_vacc_1_uv))
    
    ## Event for control occurs before paired vaccination
    # If event for potential match was before the vaccinated case was vaccinated 
    # then these are ineligible for matching
    # and therefore should be excluded
    z_merge_ps_i <- z_merge_ps_i %>% 
      filter(is.na(admission_date_uv) | admission_date_uv > date_vacc_1_vacc)
    
    ## Control death occurs before paired vaccination
    # If potential unvaccinated match died before the vaccinated case was vaccinated
    #ineligible for matching
    # admission_date - date of death
    z_merge_ps_i <- z_merge_ps_i %>% 
      filter(is.na(any_death_date_uv) | any_death_date_uv > date_vacc_1_vacc)
    
    
    ## Randomly assign match out of leftovers
    # Assings a random ID to all rows and arranges from highest to lowest
    # Once reordered, take a unique row per vaccinated ID
    set.seed(123)
    z_merge_ps_i <- z_merge_ps_i %>% 
      mutate(random_id = runif(nrow(z_merge_ps_i))) %>% 
      arrange(random_id) %>% #randomly rearrange to make sure that the same controls are not always selected
      filter(!duplicated(EAVE_LINKNO_vacc)) %>%  #get one match per vaccinated individuals
      dplyr::select(-random_id)
    
    
    #### Assign to list for loop
    z_merge_ps_list[[i]] <- z_merge_ps_i
    
    # Print progress
    print(paste0("Time-period: ", j, " - Matching batch ", i, " complete (", 40+(i-1)*10, "%)"))

  }
  
  ## Merge batched data into the one dataframe
  
  # Make into dataframe
  z_merge_ps_data <- z_merge_ps_list %>%
    reduce(full_join)
  
  # Print progress
  print(paste0("Time-period: ", j, " - Matched data combined for time period (90%)"))
  
  # Assign to list z_merge_month
  z_merge_month[[j]] <- z_merge_ps_data
  
  # Print progress
  print(paste0("Time-period: ", j, " - Complete (100%)"))
  
}


### Full join data
df_matches <- z_merge_month %>%
  reduce(full_join)

# Randomly sample matches when they are used more than the multiplicity limit
limit_matches <- function(df_matches, limit){
  
  if (limit != 'inf'){
  df_matches_less <- df_matches %>%
    group_by(EAVE_LINKNO_uv) %>%
    mutate(match_multiplicity = n()) %>%
    ungroup()
  
  df_matches_more <- filter(df_matches_less, match_multiplicity > limit)
  
  df_matches_less <- filter(df_matches_less, match_multiplicity <= limit)
  
  df_matches_more <- df_matches_more %>%
                            group_by(EAVE_LINKNO_uv) %>%
                            sample_n(limit) %>%
                            ungroup()
  
  df_matches <- bind_rows(df_matches_less, df_matches_more)
  
  saveRDS(df_matches, paste0("./data/df_matches_first_dose_", limit, '_', z_event_endpoint,".rds"))
  }
}

##### 3 - Checks ####
# Unique?
nrow(df_matches) == length(unique(df_matches$EAVE_LINKNO_vacc))

# % of matches that have been vaccinated
nrow(df_matches)/length(which(z_chrt$vacc==1))



#### 4 - Output ####


for (limit in multiplicity_limits){
  limit_matches(df_matches, limit)
}

rm(z_merge_ps_list, z_merge_ps_i, z_merge_ps_data, z_merge_month, z_df, z_df_i, m_j, loop_breaks5)
