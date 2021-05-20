##########################################################
## Title: 1st dose COVID-19 vaccine waning
## Code author(s): Rachel Mulholland <rachel.mulholland@ed.ac.uk> 
##                 Chris Robertson <chrisobertson@nhs.net>
## Description: 01_data_setup - Sets up the two baseline datasets
##              df_cohort and df_vaccinations
##########################################################

#### 0 - Set-up #####

# Libraries
library(tidyverse)
library(mgcv)
library(lubridate)

# Server location
Location <- "/conf/"  # Server

# Cohort start date - 8th Dec (when vaccines first rolled out in Scotland)
a_begin <- as.Date("2020-12-08")


##### 1 - Baseline characteristic data (df_cohort)  ####

## EAVE Cohort - Endpoints as of 20th April 2021
# Load in End points dataset
EAVE_cohort <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/Cohort_Demog_Endpoints_Times2021-04-20.rds")) %>%
  dplyr::filter(!duplicated(EAVE_LINKNO))%>%
  #remove all who have died before the beginning
  dplyr::filter(is.na(NRS.Date.Death) | (!is.na(NRS.Date.Death) & 
                                           NRS.Date.Death >= a_begin)) %>%
  #remove under 18s
  dplyr::filter(ageYear >= 18) %>%
  #adjust inconsistencies in the endpoints and times - all hosp have an admission date
  mutate(death_covid = case_when(death_covid==1 & 
                                   is.na(NRS.Date.Death) ~ 0,
                                 TRUE ~ death_covid),
         icu_death = case_when(icu_death==1 & 
                                 is.na(date_icu_death) ~ 0,
                               TRUE ~ icu_death ) )


## Weights
EAVE_Weights <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_Weights.rds"))
# Join onto EAVE cohort
EAVE_cohort  <- EAVE_cohort %>% left_join(EAVE_Weights, by="EAVE_LINKNO")
# If missing, impute with the mean weight
EAVE_cohort$eave_weight[is.na(EAVE_cohort$eave_weight)] <- mean(EAVE_cohort$eave_weight, na.rm=T)


##  Risk groups
# EAVE II BP and smoker status
eave_rg_bpsmoke <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_RG_EAVE_BP_Smoke.rds")) %>%
  filter(!duplicated(EAVE_LINKNO))%>% 
  dplyr::select(EAVE_LINKNO, EAVE_Smoking_Status_Worst, EAVE_BP) %>% 
  dplyr::rename(EAVE_Smoke = EAVE_Smoking_Status_Worst)

# QCOVID risk groups
qcovid_rg <- readRDS("/conf/EAVE/GPanalysis/progs/CR/Vaccine/output/temp/Qcovid.rds") %>% 
  dplyr::select(-(Sex:age_gp), -Q_BMI) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  # Categorise BMI into groups
  mutate(bmi_cat = cut(bmi_impute, breaks = c(-Inf, 18.5,24.9,29.9,Inf),
                       labels=c("Underweight","Normal weight","Overweight","Obese")))

## Number of previous tests
prev_tests  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/Tests.rds")) %>% 
  dplyr::select(EAVE_LINKNO, n_tests)

# Remove duplicates - take the row with the highest number of tests
prev_tests <- prev_tests %>%
  group_by(EAVE_LINKNO) %>%
  summarise(n_tests = max(n_tests)) %>% 
  mutate(n_tests = if_else(is.na(n_tests),0L,n_tests) )



## Previous positive tests
pos_tests <- readRDS(paste0(Location,"EAVE/GPanalysis/data/Positive_Tests.rds")) %>%
  group_by(EAVE_LINKNO) %>%
  # Take minimum specimen date if multiple
  summarise(specimen_date = min(specimen_date)) %>%
  # Calculate days from start date to specimen date
  mutate(days = as.numeric(specimen_date - a_begin)) %>%
  # Group days into categories
  mutate(test_before_dec8 = cut(days, breaks = c((min(days)-1), -28, -21, -14, -7, 0, max(days)),
                                labels=c("1m+", "4w","3w","2w","0-6d","after-8dec"))) %>%
  mutate(test_before_dec8 = as.character(test_before_dec8))


## Geography lookup
datazones <- read_csv(paste0(Location,"/EAVE/GPanalysis/data/restored/map_files/Datazone2011Lookup.csv")) %>% 
  dplyr::select(DataZone, InterZone, Council, HB)



## Link all datasets
df_cohort <- EAVE_cohort %>% 
  # Select EAVE_cohort variables
  dplyr::select(EAVE_LINKNO:ur6_2016_name, age_gp, -eave_weight) %>% 
  # Link to RG and BPsmoke
  left_join(eave_rg_bpsmoke %>%
              select(-EAVE_Smoke, -EAVE_BP), by="EAVE_LINKNO") %>% 
  # Link QCOVID risk groups
  left_join(qcovid_rg, by="EAVE_LINKNO") %>% 
  # Link Previous tests
  left_join(prev_tests, by="EAVE_LINKNO") %>%
  # If NA then 0 tests were taken
  mutate(n_tests = if_else(is.na(n_tests),0L,n_tests) ) %>%
  # Group the number of tests into categories
  mutate(n_tests_gp = cut(n_tests, breaks = c(-1,0,1,2,3,9,100), 
                          labels=c("0","1","2","3","4-9","10+"))) %>%
  # Link previous positive tests
  left_join(pos_tests %>%
              dplyr::select(EAVE_LINKNO, test_before_dec8), by="EAVE_LINKNO") %>%
  # If NA positive tests then no positive tests
  mutate(test_before_dec8 = if_else(is.na(test_before_dec8), "no pos test",test_before_dec8)) %>%
  # Like Data zones
  left_join(datazones, by="DataZone") %>% 
  mutate(HB = if_else(is.na(HB),"Unknown", HB),
         InterZone = if_else(is.na(InterZone),"Unknown", InterZone),
         Council = if_else(is.na(Council),"Unknown", Council)) %>%
  # Group age
  mutate(age_gp_5yr = cut(ageYear, breaks = c(seq(5, 80, by=5), 100)))


# Quick checks on column names
colnames(df_cohort)
colnames(EAVE_cohort)
colnames(eave_rg_bpsmoke)
colnames(qcovid_rg)


#### 2 - Vaccine data (df_vaccinations) ####
# Inital load
Vaccinations  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/cleaned_data/C19vaccine.rds")) %>% 
  mutate(Date = as.Date(occurrence_time)) %>% 
  mutate(vacc_type = case_when(grepl("COURAGEOUS", type) ~ "PB",
                               grepl("TALENT", type) ~ "AZ",
                               type == "39114911000001105" ~ "AZ",
                               type == "39115611000001103" ~ "PB",
                               #type == "39115611000001105" ~ "PB",
                               TRUE ~ "UNK"), 
         dose_number = if_else(stage %in% c(0,3), 1L, stage))

# 1st dose
Vaccinations1 <- filter(Vaccinations, dose_number==1) %>% 
  dplyr::select(EAVE_LINKNO, Date, vacc_type, dose_number) %>% 
  arrange(EAVE_LINKNO, Date) %>% 
  filter(!duplicated(EAVE_LINKNO))

# 2nd dose
Vaccinations2 <- filter(Vaccinations, dose_number==2) %>% 
  dplyr::select(EAVE_LINKNO, Date, vacc_type, dose_number) %>% 
  arrange(EAVE_LINKNO, Date) %>% 
  filter(!duplicated(EAVE_LINKNO))

# Join 1st and 2nd dose as separate columns
df_vaccinations <- left_join(Vaccinations1,Vaccinations2, by="EAVE_LINKNO") %>% 
  mutate(date_vacc_1 = as.Date(Date.x), 
         date_vacc_2 = as.Date(Date.y) ) %>% 
  dplyr::rename(vacc_type=vacc_type.x,
                vacc_type_2=vacc_type.y) %>% 
  dplyr::select(-dose_number.x, -dose_number.y, -Date.x, -Date.y) %>%
  #omit inconsistent records
  filter(vacc_type %in% c("AZ","PB")) %>% 
  filter(vacc_type_2 %in% c("AZ","PB") | is.na(vacc_type_2)) %>% 
  filter( !(!is.na(vacc_type_2) & (vacc_type_2 != vacc_type))) %>%
  #Get week of vaccination
  mutate(week_start_vacc1 = floor_date(date_vacc_1, unit = "week", 
                                       week_start = 1)) %>%
  mutate(week_start_vacc2 = floor_date(date_vacc_2, unit = "week", 
                                       week_start = 1))



##### 3 - Key outcomes data #####
# All dates of events will be named admission_date i.e.
# admission_date = Date of event (hospital admission date and/or date of death)
# Each dataset has the IDs, specimen date and event date of those with an outcome

## COVID-19 hospitalisations
covid_hospitalisations <- EAVE_cohort %>% 
  dplyr::select(EAVE_LINKNO, SpecimenDate, hosp_covid, date_hosp_covid, NRS.Date.Death) %>% 
  filter(hosp_covid==1) %>% 
  filter(date_hosp_covid >= a_begin) %>% 
  dplyr::rename(admission_date = date_hosp_covid) %>% 
  mutate(admission_date = if_else(is.na(NRS.Date.Death) | !is.na(NRS.Date.Death)&(admission_date <= NRS.Date.Death), admission_date, NRS.Date.Death)) %>% 
  dplyr::select(-hosp_covid, -NRS.Date.Death)


## COVID-19 ICU or deaths (severe cases)
covid_icu_death <- EAVE_cohort %>% 
  dplyr::select(EAVE_LINKNO, SpecimenDate, icu_death, date_icu_death, NRS.Date.Death) %>% 
  filter(icu_death==1) %>% 
  filter(date_icu_death >= a_begin) %>% 
  dplyr::rename(admission_date = date_icu_death) %>% 
  mutate(admission_date = if_else(is.na(NRS.Date.Death) | !is.na(NRS.Date.Death)&(admission_date <= NRS.Date.Death), admission_date, NRS.Date.Death)) %>% 
  dplyr::select(-icu_death, -NRS.Date.Death)


## COVID-19 deaths
covid_death <- EAVE_cohort %>% 
  dplyr::select(EAVE_LINKNO, SpecimenDate, death_covid, NRS.Date.Death) %>% 
  filter(death_covid==1) %>% 
  filter(NRS.Date.Death >= a_begin) %>% 
  dplyr::rename(admission_date = NRS.Date.Death) %>% 
  dplyr::select(-death_covid)


## Any deaths
any_death <- EAVE_cohort %>% 
  dplyr::select(EAVE_LINKNO, SpecimenDate, death_covid, NRS.Date.Death) %>% 
  filter(!is.na(NRS.Date.Death)) %>% 
  filter(NRS.Date.Death >= a_begin) %>% 
  dplyr::rename(admission_date = NRS.Date.Death) %>% 
  dplyr::select(-death_covid)

## COVID-19 hospitalisations and deaths (composite outcome)
covid_hosp_death <- covid_death %>%
  dplyr::rename(NRS.Date.Death= admission_date) %>% # 
  full_join(covid_hospitalisations) %>%
  dplyr::rename(hosp_admission_date = admission_date) %>%
  # If hospital admission happens before death then put death date
  # Fill in the remaining with the non-NA dates
  mutate(admission_date = if_else(hosp_admission_date < NRS.Date.Death, hosp_admission_date,
                                  NRS.Date.Death)) %>%
  mutate(admission_date = if_else(is.na(hosp_admission_date), NRS.Date.Death,
                                  hosp_admission_date)) %>%
  mutate(outcome_date = if_else(hosp_admission_date < NRS.Date.Death, "hosp",
                                "death")) %>%
  mutate(outcome_date = if_else(is.na(hosp_admission_date), "death",
                                "hosp"))



## COVID-19 positive test
positive_test <- EAVE_cohort %>% 
  dplyr::select(EAVE_LINKNO, SpecimenDate, result) %>%
  filter(result == 1) %>%
  #filter(SpecimenDate >= a_begin) %>% 
  mutate(admission_date = SpecimenDate) %>% # So consistent with other datasets
  select(-result)

df_cohort <- df_cohort %>%
  #select(-SpecimenDate) %>%
  left_join(positive_test %>%
              select(SpecimenDate, EAVE_LINKNO), by=c("EAVE_LINKNO" = "EAVE_LINKNO"))




##### 4 - Hospitalisation data ####
# All hospitalisations
all_hospitalisations  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/any_hospitalisation_post_01022020 (2).RDS")) %>%
  left_join(covid_hospitalisations, by="EAVE_LINKNO", suffix=c("","_covid")) %>% 
  filter(is.na(admission_date_covid) | !is.na(admission_date_covid)&(admission_date_covid > admission_date + 14) ) %>% 
  dplyr::select(-SpecimenDate, - admission_date_covid)

# All emergency hospitalisations
all_nc_emerg_hosp <- filter(all_hospitalisations, emergency==TRUE) %>% 
  dplyr::select(-emergency)

# All non-emergency hospitalisations
all_nc_non_emerg_hosp <- filter(all_hospitalisations, emergency==FALSE) %>% dplyr::select(-emergency)

# Creating intervals of hospitalisation data
df_intervals <- all_hospitalisations %>% 
  mutate(intervals = interval(start = admission_date, end = discharge_date)) %>% 
  select(EAVE_LINKNO, intervals) %>% 
  group_by(EAVE_LINKNO) %>% 
  mutate(interval_num = row_number()) %>% 
  pivot_wider(id_cols = EAVE_LINKNO, values_from = intervals, names_from = interval_num, names_prefix = "adm_")








##### 5- Update weights #####
# Update all the weights of the IDs that have had an interaction with the healthcare service to 1

# Filter out those without an EAVE Weight
df_cohort <- filter(df_cohort, !is.na(eave_weight)) #omit any who - need to fix -  do not match


# Find everyone with an interaction 
z_ids <- c(Vaccinations$EAVE_LINKNO, 
           any_death$EAVE_LINKNO, 
           covid_hospitalisations$EAVE_LINKNO, 
           filter(EAVE_cohort, tested==1)$EAVE_LINKNO) %>% 
  unique()

# Adjust the weights to 1 and downweight the rest
z_N <- round(sum(df_cohort$eave_weight) )
z_k <- sum(df_cohort$EAVE_LINKNO %in% z_ids)
z_m <- round(sum(filter(df_cohort, (EAVE_LINKNO %in% z_ids))$eave_weight))
z <- df_cohort %>% mutate(ew = if_else(EAVE_LINKNO %in% z_ids, 1, eave_weight*(z_N - z_k)/(z_N - z_m)) )

# Link new weights
df_cohort <- z %>% dplyr::select(-eave_weight) %>% dplyr::rename(eave_weight=ew)

# Check sum of the weights (more accurate population size)
sum(df_cohort$eave_weight)



##### 6 - Dataset checking ####
summary(df_cohort)
table(df_cohort$test_before_dec8)
table(df_cohort$EAVE_Smoke)
table(df_cohort$EAVE_BP)

# Check for duplicate rows
nrow(df_cohort) == length(unique(df_cohort$EAVE_LINKNO))



##### 7 - Save to output folder ####
# Save df_cohort
saveRDS(df_cohort, paste0("./output/df_cohort.rds"))

# Save df_vaccinations
saveRDS(df_vaccinations, paste0("./output/df_vaccinations.rds"))



##### 8- Remove files no longer needed ####
rm(EAVE_cohort, eave_rg_bpsmoke, EAVE_Weights, pos_tests, prev_tests, datazones,
   qcovid_rg, Vaccinations1, Vaccinations2,Vaccinations, z, z_ids, z_N, z_k, z_m)
