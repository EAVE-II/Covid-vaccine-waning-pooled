##########################################################
## Title: 1st dose COVID-19 vaccine waning
## Code author(s): Rachel Mulholland <rachel.mulholland@ed.ac.uk> 
##                 Chris Robertson <chrisobertson@nhs.net>
## Description: 03b_matching_cc - Formats the matched data into the
##              case-control (cc) format
##########################################################

#### 0 - Set up ####

# Choose which multiplicity limit we're using
multiplicity_limit <- 1

## Load in matched and vaccinated data

# Event
z_event_endpoint <- "death_hosp"
#z_event_endpoint <- "positive_test"
#z_event_endpoint <- "hosp_covid"
#z_event_endpoint <- "death_covid"

### Load in data based on endpoint
if (z_event_endpoint =="hosp_covid") {z_event <- covid_hospitalisations
z_merge <- readRDS(paste0("./data/df_matches_second_dose_", multiplicity_limit, '_', z_event_endpoint,".rds"))
}
if (z_event_endpoint =="death_hosp") {z_event <- covid_hosp_death
z_merge <- readRDS(paste0("./data/df_matches_second_dose_", multiplicity_limit, '_', z_event_endpoint,".rds"))
}
if (z_event_endpoint =="positive_test") {z_event <- positive_test
z_merge <- readRDS(paste0("./data/df_matches_second_dose_", multiplicity_limit, '_', z_event_endpoint,".rds"))

}
if (z_event_endpoint =="death_covid") {z_event <- covid_death
z_merge <- readRDS(paste0("./data/df_matches_second_dose_", multiplicity_limit, '_', z_event_endpoint,".rds"))
}


## Vaccination data
# Make anyone first dose vaccinated after the maximum endpoint time unvaccinated
z_vaccinations <- filter(df_vaccinations, date_vacc_1 <= a_end) %>% 
  mutate(vacc_type_2 = if_else(date_vacc_2 > a_end, NA_character_ , vacc_type_2),
         date_vacc_2 = as.Date(ifelse(date_vacc_2 > a_end, NA, date_vacc_2), origin=as.Date("1970-01-01")) )


##### 1 - Separate cases and controls ####

# Cases
z_cc_v2 <- z_merge %>% 
  dplyr::select(EAVE_LINKNO_v2, vacc_type, date_vacc_2_v1, date_vacc_2_v2, Council,
                 admission_date_v1, admission_date_v2) %>% 
  mutate(EAVE_LINKNO = EAVE_LINKNO_v2) %>% 
  relocate(EAVE_LINKNO) 

# Controls
z_cc_v1 <- z_merge %>% 
  dplyr::select(EAVE_LINKNO_v1, EAVE_LINKNO_v2, vacc_type, date_vacc_2_v1, date_vacc_2_v2, 
                Council, admission_date_v1, admission_date_v2) %>% 
  dplyr::rename(EAVE_LINKNO = EAVE_LINKNO_v1)

##### 2 - Combine cases and controls together ####

## Stack vaccinated and unvaccinated on top
z_cc <- bind_rows(z_cc_v2, z_cc_v1) %>% 
  arrange(EAVE_LINKNO_v2) %>% 
  # Create vacc flag if EAVE ID = Vacc EAVE ID
  mutate(vacc = factor(if_else(EAVE_LINKNO==EAVE_LINKNO_v2, 1,0), 
                       levels=c(0,1),labels=c("v1","v2"))) %>%
  relocate(EAVE_LINKNO)

nrow(z_cc) == nrow(z_cc_v2) + nrow(z_cc_v1)

#### 3 - Eliminate matches with events that are ineligible ####
# For data flow diagram (run first to get numbers of pairs where at least one had a covid hosp/death on the same day or before vacc exposed)
z_cc2 <- z_cc %>%
  filter(admission_date_v2 <= date_vacc_2_v2 | admission_date_v1 <= date_vacc_2_v2 ) %>%
  select(-admission_date_v2, -admission_date_v1) %>%
  left_join(select(df_cohort, EAVE_LINKNO, eave_weight))


## Eliminate those with event on the same day or before vaccination date for both pairs (NEW)
z_cc <- z_cc %>%
  filter(admission_date_v2 > date_vacc_2_v2 | is.na(admission_date_v2)) %>%
  filter(admission_date_v1 > date_vacc_2_v2 | is.na(admission_date_v1)) %>%
  select(-admission_date_v2, -admission_date_v1)


#sum(z_cc2$eave_weight)
nrow(z_cc2)/2



#### 4 - Create time to event #####

## Initalise
z_cc <- z_cc %>% 
  # Initialise event date with end date
  mutate(event_date = a_end) 

## Merge event data
# Merge in the event date to calculate the event and event date
z_cc <- z_cc %>% 
  # Merge event date
  left_join(z_event, by="EAVE_LINKNO") %>%
  # If event occurs (admission date isnt NA) then put 1
  mutate(event = if_else(!is.na(admission_date), 1,0)) %>% 
  # Adjust event date if event happens before event date
  mutate(event_date = if_else( (event==1) & (admission_date <= event_date), 
                               admission_date, event_date)) %>%
  dplyr::select(-(SpecimenDate))

## Link death data
#link in any death and modify the event date.  There should not be any events to change as event date <= death date.  
z_cc <- z_cc %>%
  #Link death data
  left_join(any_death %>%
              rename("any_death_date" = "admission_date"), by="EAVE_LINKNO") %>%
  #Adjust event date if person dies before end date
  mutate(event_date = if_else(!is.na(any_death_date) & (any_death_date < event_date),
                              any_death_date, event_date))%>%
  dplyr::select(-SpecimenDate)

# Create a non-vaccine censored time to event. This is for calculating person years
# lost due to censoring because of control receiving second dose.
z_cc <- mutate(z_cc, event_date_non_vacc_censored = event_date) %>% 
  mutate(event_non_vacc_censored = if_else(event==1 & event_date < admission_date, 0, event)) 

## Censor if control receives second dose
z_cc <- z_cc %>% 
  # relocate the follow up to the date the control was second-dose vaccinated
  mutate(event_date = if_else(!is.na(date_vacc_2_v1) & (date_vacc_2_v1 < event_date), 
                              date_vacc_2_v1, event_date))

## Adjust event flag
#change the event marker from 1 to 0 for those whose admission_date is greater than the current event_date
z_cc <- z_cc %>% 
  mutate(event = if_else(event==1 & event_date < admission_date, 0, event)) 



## Adjust if endpoint is positive test
if(z_event_endpoint == "positive_test"){
  ## If positive infection as outcome, then look at those tested positive before 8th Dec
  positive_test_all <- EAVE_cohort %>% 
    dplyr::select(EAVE_LINKNO, SpecimenDate, result) %>%
    filter(result == 1) %>%
    select(-result)
  
  # Filter to only those who did not test positive or tested postive after second dose vaccination
  z_cc <- z_cc %>% 
    left_join(positive_test_all)  %>%
    rename(SpecimenDate_v1 = SpecimenDate) %>%
    left_join(positive_test_all, by=c("EAVE_LINKNO_v2"="EAVE_LINKNO"))%>%
    rename(SpecimenDate_v2 = SpecimenDate) %>%
    filter(SpecimenDate_v2 > date_vacc_2_v2 | is.na(SpecimenDate_v2))
  
  
  z_id <- z_cc %>%
    filter(SpecimenDate_v1 <= date_vacc_2_v2) %>%
    pull(EAVE_LINKNO) %>%
    unique()
  
  z_cc <- filter(z_cc, 
                 !(EAVE_LINKNO_v2 %in% z_id))
  
}


### Calculate time until event ###
# Time until event is the event date - vaccination date
z_cc <- z_cc %>%  
  # Whole time to event
  mutate(time_to_event = as.numeric(event_date-date_vacc_2_v2)) %>%
  # Time to event starting at 14 days
  mutate(time_to_event14 = ifelse(time_to_event < 14, NA, time_to_event)) %>%
  # Time to non-vacine censored event
  mutate(time_to_event_non_vacc_censored = as.numeric(event_date_non_vacc_censored-date_vacc_2_v2)) %>%
  # Put into time periods (currently up to 12 weeks)
  mutate(period = cut(time_to_event, 
                      breaks= c(-1,13,20,27,34,41, 48, 55, 62, 69, 76, 83, max(time_to_event, na.rm=T)),
                      labels=c("0:13","14:20","21:27","28:34","35:41", "42:47", "49:55", "56:62", "63:69", 
                               "70:76", "77:83", "84+"))) 



##### 5 - Fix data errors #####
#time on study possible negative for data errors - omit both vacc and matched unvacc
# Find negative time to events (i.e. event happened before vaccination = error)
# Might be people vaccinated in hospital - Check
z_errors <- filter(z_cc, time_to_event <0) # n = 11 
nrow(z_errors)

z_errors_ids <- unique(z_errors$EAVE_LINKNO_v2)
length(z_errors_ids)

# Delete from dataset
z_cc <- filter(z_cc, 
               !(EAVE_LINKNO_v2 %in% z_errors_ids))



#### 6 - Add in characteristics #####
df_cc <- z_cc %>% left_join(select(df_cohort, 
                   EAVE_LINKNO, age_gp, simd2020_sc_quintile, ageYear,
                   n_tests_gp, n_risk_gps, ur6_2016_name), by=c("EAVE_LINKNO" = "EAVE_LINKNO"))%>%
  mutate(age_grp = case_when(ageYear < 65 ~"18-64", 
                             ageYear < 80 ~"65-79",
                             TRUE ~ "80+"))

## Add in information about 2nd dose vaccinated person only
df_cc <- df_cc %>%
  #Age
  left_join(select(df_cohort, EAVE_LINKNO, ageYear),
            by=c("EAVE_LINKNO_v2" = "EAVE_LINKNO"), 
            suffix = c("", "_vacc"))


#### 7 - Repeat for hospitalisations and deaths separately (only if using composite outcomes) ####
# Repeats the above process for hospitalisations and deaths separately as unique time to event
# variables for the two outcomes

## Add in events for hospitalisations and deaths
df_cc <- df_cc %>% 
  mutate(event_date2 = a_end) %>%
  # relocate the follow up to the date the first dose vaccinated control received 2nd dose
  mutate(event_date2 = if_else(!is.na(date_vacc_2_v1) & (date_vacc_2_v1 < event_date2), 
                               date_vacc_2_v1, event_date2))

## Link death data
#link in any death and modify the event date.  There should not be any events to change as event date <= death date.  
df_cc <- df_cc %>%
  #Adjust event date if person dies before end date
  mutate(event_date2 = if_else(!is.na(any_death_date) & (any_death_date < event_date2),
                               any_death_date, event_date2))

## Merge event data -for hospitalisations and deaths

# Merge in the event date to calculate the event and event date
df_cc <- df_cc %>%
  # If event occurs (admission date isnt NA) then put 1
  mutate(event_hosp = if_else(!is.na(hosp_admission_date), 1,0)) %>% 
  # Adjust event date if event happens before event date
  mutate(event_date_hosp = if_else( (event_hosp==1) & (hosp_admission_date <= event_date2), 
                                    hosp_admission_date, event_date2)) %>%
  # If event occurs (admission date isnt NA) then put 1
  mutate(event_death = if_else(!is.na(NRS.Date.Death), 1,0)) %>% 
  # Adjust event date if event happens before event date
  mutate(event_date_death = if_else( (event_death==1) & (NRS.Date.Death <= event_date2), 
                                     NRS.Date.Death, event_date2)) 

### Calculate time to events

df_cc <- df_cc %>%  
  # Hosp
  mutate(time_to_hosp = as.numeric(event_date_hosp-date_vacc_2_v2)) %>%
  mutate(time_to_event14_hosp = ifelse(time_to_hosp < 14, NA, time_to_hosp)) %>%
  mutate(period_hosp = cut(time_to_hosp, 
                           breaks= c(-1,13,20,27,34,41, 48, 55, 62, 69, 76, 83, max(time_to_hosp, na.rm=T)),
                           labels=c("0:13","14:20","21:27","28:34","35:41", "42:47", "49:55", "56:62", "63:69", 
                                    "70:76", "77:83", "84+"))) %>%
  # Death
  mutate(time_to_death = as.numeric(event_date_death-date_vacc_2_v2)) %>%
  mutate(time_to_event14_death = ifelse(time_to_death < 14, NA, time_to_death)) %>%
  mutate(period_death = cut(time_to_death, 
                            breaks= c(-1,13,20,27,34,41, 48, 55, 62, 69, 76, 83, max(time_to_death, na.rm=T)),
                            labels=c("0:13","14:20","21:27","28:34","35:41", "42:47", "49:55", "56:62", "63:69", 
                                     "70:76", "77:83", "84+")))

## Adjust event flag for censoring
## If censored before event, event is 0
df_cc <- df_cc %>% 
  mutate(event_hosp = if_else(event_hosp==1 & event_date_hosp < hosp_admission_date, 0, event_hosp)) %>%
  mutate(event_death = if_else(event_death==1 & event_date_death < NRS.Date.Death, 0, event_death))


# Person years and events lost due to vaccine censoring
# Censored
z.agg <- pyears(Surv(time_to_event,event) ~ total ,
            data=df_cc %>% mutate(total = 'Censored') , scale=1, data.frame=TRUE)

df_res1 <- z.agg$data
df_res1 <- df_res1 %>% 
  mutate(pyears =round(pyears/365.25,1)) %>%
  select(-n) 

# Uncensored
z.agg <- pyears(Surv(time_to_event_non_vacc_censored, event_non_vacc_censored) ~ total,
                data=df_cc %>% mutate(total = 'Uncensored') , scale=1, data.frame=TRUE)

df_res2 <- z.agg$data
df_res2 <- df_res2 %>% 
  mutate(pyears =round(pyears/365.25,1)) %>%
  select(-n) 

# Combine
df_res <- bind_rows(df_res1, df_res2) %>%
  mutate_if(is.numeric, ~formatC(., format = "f", big.mark = ",", drop0trailing = TRUE))

names(df_res) <- c('', 'Person years', 'Events')

write.csv(df_res, paste0("./output/second_dose_", multiplicity_limit, "/final/matching_summary/pyears_lost.csv"))


##### 8 - Save as rds ####
saveRDS(df_cc, paste0("./data/df_cc_second_dose_", multiplicity_limit, '_',
                      z_event_endpoint, ".rds"))

rm(z_cc, z_cc_v1, z_cc_v2, df_cc, z_merge)
