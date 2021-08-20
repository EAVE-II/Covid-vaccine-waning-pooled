##########################################################
## Title: 1st dose COVID-19 vaccine waning
## Code author(s): Rachel Mulholland <rachel.mulholland@ed.ac.uk> 
##                 Chris Robertson <chrisobertson@nhs.net>
## Description: 03b_matching_cc - Formats the matched data into the
##              case-control (cc) format
##########################################################

#### 0 - Set up ####

## Load in matched and vaccinated data

# Event
z_event_endpoint <- "death_hosp"
#z_event_endpoint <- "positive_test"
#z_event_endpoint <- "hosp_covid"
#z_event_endpoint <- "death_covid"

### Load in data based on endpoint
if (z_event_endpoint =="hosp_covid") {z_event <- covid_hospitalisations
z_merge <- readRDS("./output/df_matches_death_hosp.rds")
}
if (z_event_endpoint =="death_hosp") {z_event <- covid_hosp_death
z_merge <- readRDS("./output/df_matches_death_hosp.rds")
}
if (z_event_endpoint =="positive_test") {z_event <- positive_test
z_merge <- readRDS("./output/df_matches_death_hosp.rds")

}
if (z_event_endpoint =="death_covid") {z_event <- covid_death
z_merge <- readRDS("./output/df_matches_death_hosp.rds")
}


## Vaccination data
# Make anyone vaccinated after the maximum endpoint time unvaccinated
z_vaccinations <- filter(df_vaccinations, date_vacc_1 <= a_end) %>% 
  mutate(vacc_type_2 = if_else(date_vacc_2 > a_end, NA_character_ , vacc_type_2),
         date_vacc_2 = as.Date(ifelse(date_vacc_2 > a_end, NA, date_vacc_2), origin=as.Date("1970-01-01")) )


##### 1 - Separate cases and controls ####

# Cases
z_cc_vacc <- z_merge %>% 
  dplyr::select(EAVE_LINKNO_vacc, date_vacc_1_vacc, date_vacc_1_uv, prop_score_vacc, admission_date_vacc, admission_date_uv) %>% 
  mutate(EAVE_LINKNO = EAVE_LINKNO_vacc) %>% 
  relocate(EAVE_LINKNO) %>%
  rename(prop_score=prop_score_vacc)


# Controls
z_cc_uv <- z_merge %>% 
  dplyr::select(EAVE_LINKNO_uv, EAVE_LINKNO_vacc, date_vacc_1_vacc, 
                date_vacc_1_uv, prop_score_uv, admission_date_vacc, admission_date_uv) %>% 
  dplyr::rename(EAVE_LINKNO = EAVE_LINKNO_uv, prop_score=prop_score_uv) 



##### 2 - Combine cases and controls together ####

## Stack vaccinated and unvaccinated on top
z_cc <- bind_rows(z_cc_vacc, z_cc_uv) %>% 
  arrange(EAVE_LINKNO_vacc) %>% 
  # Create vacc flag if EAVE ID = Vacc EAVE ID
  mutate(vacc = factor(if_else(EAVE_LINKNO==EAVE_LINKNO_vacc, 1,0), 
                       levels=c(0,1),labels=c("uv","vacc"))) %>%
  relocate(EAVE_LINKNO)

nrow(z_cc) == nrow(z_cc_vacc) + nrow(z_cc_uv)



#### 3 - Eliminate matches with events that are ineligible ####
# For data flow diagram (run first to get numbers of pairs where at least one had a covid hosp/death on the same day or before vacc exposed)
z_cc2 <- z_cc %>%
  filter(admission_date_vacc <= date_vacc_1_vacc | admission_date_uv <= date_vacc_1_vacc) %>%
  select(-admission_date_vacc, -admission_date_uv) %>%
  left_join(select(df_cohort, EAVE_LINKNO, eave_weight))


## Eliminate those with event on the same day or before vaccination date for both pairs (NEW)
z_cc <- z_cc %>%
  filter(admission_date_vacc > date_vacc_1_vacc | is.na(admission_date_vacc)) %>%
  filter(admission_date_uv > date_vacc_1_vacc | is.na(admission_date_uv)) %>%
  select(-admission_date_vacc, -admission_date_uv)


#sum(z_cc2$eave_weight)
nrow(z_cc2)/2




#### 4 - Create time to event #####

## Initalise
z_cc <- z_cc %>% 
  # Initialise event date with end date
  mutate(event_date = a_end) 

## Censor if unvaccinated becomes vaccinated
z_cc <- z_cc %>% 
  # relocate the follow up to the date the unvaccinated control was vaccinated - Censoring?
  mutate(event_date = if_else(!is.na(date_vacc_1_uv) & (date_vacc_1_uv < event_date), 
                              date_vacc_1_uv, event_date))

## Censor if vaccinated receives 2nd dose
z_cc <- z_cc %>%
  # If vaccinated person gets 2nd dose - censor
  left_join(select(df_vaccinations,EAVE_LINKNO, date_vacc_2), 
            by=c("EAVE_LINKNO_vacc" = "EAVE_LINKNO"), 
            suffix=c("","_vacc")) %>%
  # If there is a date for vaccine 2 and it is before the event date and it is after the 1st dose date
  mutate(event_date = if_else(!is.na(date_vacc_2) & # If vaccinated person has been vaccinated and
                                (date_vacc_2 < event_date) &  # they are vaccinated before the event date
                                (date_vacc_2 > date_vacc_1_vacc),  # And the date of vaccination is after the 1st dose
                              date_vacc_2, event_date)) # Then replace with 2nd dose

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
  # Link death data
  left_join(any_death %>%
              rename("any_death_date" = "admission_date"), by="EAVE_LINKNO") %>%
  # Adjust event date if person dies before end date
  mutate(event_date = if_else(!is.na(any_death_date) & (any_death_date < event_date), 
                              any_death_date, event_date))%>%  
  dplyr::select(-SpecimenDate)


## Adjust event flag
#change the event marker from 1 to 0 for those whose admission_date is greater than the current event_date
#these are instances where the unvaccinated control is vaccinated before the admission_date
#and hence censored at the vaccination date
z_cc <- z_cc %>% 
  mutate(event = if_else(event==1 & event_date < admission_date, 0, event))


## Adjust if endpoint is positive test
if(z_event_endpoint == "positive_test"){
  ## If positive infection as outcome, then look at those tested positive before 8th Dec
  positive_test_all <- EAVE_cohort %>% 
    dplyr::select(EAVE_LINKNO, SpecimenDate, result) %>%
    filter(result == 1) %>%
    select(-result)
  
  # Filter to only those who did not test positive or tested postive after the vaccination date
  z_cc <- z_cc %>% 
    left_join(positive_test_all)  %>%
    rename(SpecimenDate_uv = SpecimenDate) %>%
    left_join(positive_test_all, by=c("EAVE_LINKNO_vacc"="EAVE_LINKNO"))%>%
    rename(SpecimenDate_vacc = SpecimenDate) %>%
    filter(SpecimenDate_vacc > date_vacc_1_vacc | is.na(SpecimenDate_vacc))
  
  
  z_id <- z_cc %>%
    filter(SpecimenDate_uv <= date_vacc_1_vacc) %>%
    pull(EAVE_LINKNO) %>%
    unique()
  
  z_cc <- filter(z_cc, 
                 !(EAVE_LINKNO_vacc %in% z_id))
  
}




### Calculate time until event ###
# Time until event is the event date - vaccination date
z_cc <- z_cc %>%  
  # Whole time to event
  mutate(time_to_hosp = as.numeric(event_date-date_vacc_1_vacc)) %>%
  # Time to event starting at 14 days
  mutate(time_to_event14 = ifelse(time_to_hosp < 14, NA, time_to_hosp)) %>%
  # Put into time periods (currently up to 12 weeks)
  mutate(period = cut(time_to_hosp, 
                      breaks= c(-1,13,20,27,34,41, 48, 55, 62, 69, 76, 83, max(time_to_hosp, na.rm=T)),
                      labels=c("0:13","14:20","21:27","28:34","35:41", "42:47", "49:55", "56:62", "63:69", 
                               "70:76", "77:83", "84+"))) 



##### 5 - Fix data errors #####
#time on study possible negative for data errors - omit both vacc and matched unvacc
# Find negative time to events (i.e. event happened before vaccination = error)
# Might be people vaccinated in hospital - Check
z_errors <- filter(z_cc, time_to_hosp <0) # n = 6 (n=5 for infections)
nrow(z_errors)

z_errors_ids <- unique(z_errors$EAVE_LINKNO_vacc)
length(z_errors_ids)

# Delete from dataset
z_cc <- filter(z_cc, 
               !(EAVE_LINKNO_vacc %in% z_errors_ids))


#### 6 - Add in characteristics #####
#link to vaccination type - link both vaccinated and unvaccinated control so that
#selection on vacc type can be easily made.
df_cc <- z_cc %>%
  left_join(dplyr::select(z_vaccinations,
                          EAVE_LINKNO, vacc_type), by=c("EAVE_LINKNO_vacc" = "EAVE_LINKNO"),
            suffix = c("", "_vacc")) %>%
  left_join(select(df_cohort, 
                   EAVE_LINKNO, age_gp, ageYear, simd2020_sc_quintile,
                   n_tests_gp, n_risk_gps, Council, ur6_2016_name), by=c("EAVE_LINKNO" = "EAVE_LINKNO"),
            suffix = c("", "_vacc"))%>%
  mutate(age_grp = case_when(ageYear < 65 ~"18-64", 
                             ageYear < 80 ~"65-79",
                             TRUE ~ "80+"))

## Add in information about vaccinated person only
df_cc <- df_cc %>%
  #Age
  left_join(select(df_cohort, EAVE_LINKNO, ageYear),
            by=c("EAVE_LINKNO_vacc" = "EAVE_LINKNO"), 
            suffix = c("", "_vacc"))


#### 7 - Repeat for hospitalisations and deaths separately (only if using composite outcomes) ####
# Repeats the above process for hospitalisations and deaths separately as unique time to event
# variables for the two outcomes

## Add in events for hospitalisations and deaths
df_cc <- df_cc %>% 
  mutate(event_date2 = a_end) %>%
  # relocate the follow up to the date the unvaccinated control was vaccinated - Censoring?
  mutate(event_date2 = if_else(!is.na(date_vacc_1_uv) & (date_vacc_1_uv < event_date2), 
                               date_vacc_1_uv, event_date2))

# If vaccinated person becomes vaccinated then censor for this
df_cc <- df_cc %>%
  # If there is a date for vaccine 2 and it is before the event date and it is after the 1st dose date
  mutate(event_date2 = if_else(!is.na(date_vacc_2) & # If vaccinated person has been vaccinated and
                                 (date_vacc_2 < event_date2) &  # they are vaccinated before the event date
                                 (date_vacc_2 > date_vacc_1_vacc),  # And the date of vaccination is after the 1st dose
                               date_vacc_2, event_date2)) # Then replace with 2nd dose

## Link death data
#link in any death and modify the event date.  There should not be any events to change as event date <= death date.  
df_cc <- df_cc %>% 
  # Adjust event date if person dies before end date
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
  mutate(time_to_hosp_hosp = as.numeric(event_date_hosp-date_vacc_1_vacc)) %>%
  mutate(time_to_event14_hosp = ifelse(time_to_hosp_hosp < 14, NA, time_to_hosp_hosp)) %>%
  mutate(period_hosp = cut(time_to_hosp_hosp, 
                           breaks= c(-1,13,20,27,34,41, 48, 55, 62, 69, 76, 83, max(time_to_hosp_hosp, na.rm=T)),
                           labels=c("0:13","14:20","21:27","28:34","35:41", "42:47", "49:55", "56:62", "63:69", 
                                    "70:76", "77:83", "84+"))) %>%
  # Death
  mutate(time_to_hosp_death = as.numeric(event_date_death-date_vacc_1_vacc)) %>%
  mutate(time_to_event14_death = ifelse(time_to_hosp_death < 14, NA, time_to_hosp_death)) %>%
  mutate(period_death = cut(time_to_hosp_death, 
                            breaks= c(-1,13,20,27,34,41, 48, 55, 62, 69, 76, 83, max(time_to_hosp_death, na.rm=T)),
                            labels=c("0:13","14:20","21:27","28:34","35:41", "42:47", "49:55", "56:62", "63:69", 
                                     "70:76", "77:83", "84+")))

## Adjust event flag for censoring
## If censored before event, event is 0
df_cc <- df_cc %>% 
  mutate(event_hosp = if_else(event_hosp==1 & event_date_hosp < hosp_admission_date, 0, event_hosp)) %>%
  mutate(event_death = if_else(event_death==1 & event_date_death < NRS.Date.Death, 0, event_death))



##### 8 - Save as rds ####
saveRDS(df_cc, paste0("./output/df_cc_",
                      z_event_endpoint, ".rds"))

rm(z_cc, z_cc_uv, z_cc_vacc, df_cc, z_merge)



