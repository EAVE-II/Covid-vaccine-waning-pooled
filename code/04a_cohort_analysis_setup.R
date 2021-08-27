###### Cohort analysis ####

# Setting up cohort analysis for long data

#### Set up ####
output_list <- list()
z_event_endpoint <- "death_hosp" #hosp, death or icu_death any_nc_em

if (z_event_endpoint =="hosp_covid") z_event <- covid_hospitalisations
if (z_event_endpoint =="icu_death") z_event <- covid_icu_death
if (z_event_endpoint =="death_covid") z_event <- covid_death
if (z_event_endpoint =="any_death") z_event <- any_death
if (z_event_endpoint =="any_nc_non_emerg_adm") z_event <- all_nc_non_emerg_hosp
if (z_event_endpoint =="any_nc_emerg_adm") z_event <- all_nc_emerg_hosp
if (z_event_endpoint =="death_hosp") z_event <- covid_hosp_death %>%
  select(EAVE_LINKNO, SpecimenDate, admission_date)

# Extract end date
a_end <- as.Date("2021-04-30")

# Filter event data to end date
z_event <- z_event %>%
  filter(admission_date <= a_end)

# Event
output_list$endpoint <- z_event_endpoint
# Number of events
output_list$number_events <- nrow(z_event)
# End date
output_list$last_event <- a_end
#Take 100 controls per 1 case
z_n_controls_per_event <- 100
output_list$controls_per_event <- z_n_controls_per_event

#make anyone vaccinated after the maximum endpoint time unvaccinated
z_vaccinations <- filter(df_vaccinations, date_vacc_1 <= a_end) %>% 
  mutate(vacc_type_2 = if_else(date_vacc_2 > a_end, NA_character_ , vacc_type_2),
         date_vacc_2 = as.Date(ifelse(date_vacc_2 > a_end, NA, date_vacc_2), origin=as.Date("1970-01-01")) )

#in z_event admission date is generic and meansdate death for deaths, date icu admission for icu

print(z_event_endpoint)
print(nrow(z_event))

## Vaccinations in event data

#investigation of numbers and timings - events among vaccinated
#join events to vaccinations
z1 <- left_join(z_event, z_vaccinations, by="EAVE_LINKNO") %>% 
  mutate(vaccinated = if_else(!is.na(date_vacc_1), 1,0)) %>% 
  mutate(vacc_at_event = if_else(!is.na(date_vacc_1) & date_vacc_1 > admission_date, 0, vaccinated) )
table(z1$vaccinated, z1$vacc_at_event, exclude=NULL)

# Event IDs
z_event_ids <- unique(z_event$EAVE_LINKNO)

####  DF cohort to exclude care homes ####
# Total cohort 
z_chrt <- df_cohort %>%
  #dplyr::rename(EAVE_LINKNO_uv = EAVE_LINKNO) %>%
  mutate(ur6_2016_name = replace_na(ur6_2016_name, "NA")) %>%
  # Add in age group
  mutate(age_grp = case_when(ageYear < 65 ~"18-64", 
                             ageYear < 80 ~"65-79",
                             TRUE ~ "80+"))

# Add in household information
Cohort_Household <- readRDS("/conf/EAVE/GPanalysis/outputs/temp/Cohort_Household.rds") %>%
  mutate(n_hh_gp = cut(n_hh, breaks=c(0,1,2,5,10,30,100,max(n_hh)),
                       labels=c("1", "2", "3-5", "6-10", "11-30", "31-100", "101+")))%>% 
  mutate(ave_hh_age=if_else(is.na(ave_hh_age), mean(ave_hh_age, na.rm=T), ave_hh_age) )

# Eliminate care homes
z_chrt <- z_chrt %>%
  left_join(select(Cohort_Household, EAVE_LINKNO,
                   n_hh_gp, ave_hh_age, care_home_elderly),
            by= c("EAVE_LINKNO"= "EAVE_LINKNO")) %>% 
  filter(care_home_elderly == 0)%>% 
  dplyr::select(-care_home_elderly)

# check
sum(z_chrt$eave_weight)

# Everyone else in the cohort who did not have an event
z_rest <- filter(z_chrt, !(EAVE_LINKNO %in% z_event_ids))

# Non-event IDs
z_ids_rest <- unique(z_rest$EAVE_LINKNO)

# Take sample of non-event IDs 100 x the number of events
z_ids_rest <- sample(z_ids_rest, size=nrow(z_event)*z_n_controls_per_event, replace=FALSE)

# Combine the event IDs and the 100 controls IDs
z_ids <- c(z_event_ids, z_ids_rest)

# Filter df cohort to z_ids
z_df <- filter(z_chrt, EAVE_LINKNO %in% z_ids )

###### Run '04b_cohort_analysis_long_hosp.R' script #####




source("./04b_cohort_analysis_long_hosp.R")

table(z_vaccination_long$vacc_status, exclude=NULL)

#result is z_vaccination_long
#z_covariates is the covariate file


###### Run '04c_cohort_analysis_add_endpts.R' script #####


#add in the event end point
source("./Eles_code/04c_cohort_analysis_add_endpts.R")


#df is the data frame for analysis
#check the numbers
print(table(df$event))
print(nrow(z_event))
print(table(df$vacc_status, df$event))
print(sum(df$vacc_status != "uv" & df$event==1))


#### Add in vaccination data ####
#add in the vacccination type
df <- df %>% left_join(dplyr::select(df_vaccinations, EAVE_LINKNO, vacc_type) , by="EAVE_LINKNO") %>% 
  mutate(vacc_type = if_else(is.na(vacc_type), "uv", vacc_type))

df <- df %>% mutate(weight = if_else(event==1,1,nrow(z_chrt)/(z_n_controls_per_event*nrow(z_event)) ) )  %>% 
  mutate(weight = if_else(event==1,1, eave_weight*weight) ) 



#### Output #####

Location <- "/conf/"  # Server
project_path <- paste0(Location,"EAVE/GPanalysis/progs/RM/Vaccine/Vaccine_waning/")

#saveRDS(df, paste0(project_path,"output/df_long_",z_event_endpoint,".RDS"))
#saveRDS(df, paste0(project_path,"output/df_long_",z_event_endpoint,"_12wks.RDS"))
saveRDS(df, paste0(project_path,"output/df_long_",z_event_endpoint,"_10wks.RDS"))
saveRDS(df, paste0(".output/df_long_",z_event_endpoint,"_8wks.RDS"))
saveRDS(df, paste0("./output/df_long_death_hosp_8wks.RDS"))

