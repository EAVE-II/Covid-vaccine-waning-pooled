##### Getting entire cohort weights (before eliminating children) ####
# Load in End points dataset
EAVE_cohort_all <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/Cohort_Demog_Endpoints_Times2021-04-20.rds")) %>%
  dplyr::filter(!duplicated(EAVE_LINKNO))%>%
  #remove all who have died before the beginning
  dplyr::filter(is.na(NRS.Date.Death) | (!is.na(NRS.Date.Death) & 
                                           NRS.Date.Death >= a_begin)) %>%
  #remove under 18s
  #dplyr::filter(ageYear >= 18) %>%
  #adjust inconsistencies in the endpoints and times - all hosp have an admission date
  mutate(death_covid = case_when(death_covid==1 & 
                                   is.na(NRS.Date.Death) ~ 0,
                                 TRUE ~ death_covid),
         icu_death = case_when(icu_death==1 & 
                                 is.na(date_icu_death) ~ 0,
                               TRUE ~ icu_death ) )


## Weights
EAVE_Weights <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_Weights.rds"))
EAVE_cohort_all  <- EAVE_cohort_all %>% left_join(EAVE_Weights, by="EAVE_LINKNO")
EAVE_cohort_all$eave_weight[is.na(EAVE_cohort_all$eave_weight)] <- mean(EAVE_cohort_all$eave_weight, na.rm=T)

# Filter out those without an EAVE Weight
EAVE_cohort_all <- filter(EAVE_cohort_all, !is.na(eave_weight)) #omit any who - need to fix -  do not match


## Update weights
# EAVE hospitalisations
covid_hospitalisations <- EAVE_cohort_all %>% 
  dplyr::select(EAVE_LINKNO, SpecimenDate, hosp_covid, date_hosp_covid, NRS.Date.Death) %>% 
  filter(hosp_covid==1) %>% 
  filter(date_hosp_covid >= a_begin) %>% 
  dplyr::rename(admission_date = date_hosp_covid) %>% 
  mutate(admission_date = if_else(is.na(NRS.Date.Death) | !is.na(NRS.Date.Death)&(admission_date <= NRS.Date.Death), admission_date, NRS.Date.Death)) %>% 
  dplyr::select(-hosp_covid, -NRS.Date.Death)

# use the EAVE death cases
any_death <- EAVE_cohort_all %>% 
  dplyr::select(EAVE_LINKNO, SpecimenDate, death_covid, NRS.Date.Death) %>% 
  filter(!is.na(NRS.Date.Death)) %>% 
  filter(NRS.Date.Death >= a_begin) %>% 
  dplyr::rename(admission_date = NRS.Date.Death) %>% 
  dplyr::select(-death_covid)


#all known to exist - give a weight of 1 and downweight the rest
z_ids <- c(df_vaccinations$EAVE_LINKNO, 
           any_death$EAVE_LINKNO, 
           covid_hospitalisations$EAVE_LINKNO, 
           filter(EAVE_cohort_all, tested==1)$EAVE_LINKNO) %>% 
  unique()

z_N <- round(sum(EAVE_cohort_all$eave_weight) )
z_k <- sum(EAVE_cohort_all$EAVE_LINKNO %in% z_ids)
z_m <- round(sum(filter(EAVE_cohort_all, (EAVE_LINKNO %in% z_ids))$eave_weight))
z <- EAVE_cohort_all %>% mutate(ew = if_else(EAVE_LINKNO %in% z_ids, 1, eave_weight*(z_N - z_k)/(z_N - z_m)) )

EAVE_cohort_all <- z %>% dplyr::select(-eave_weight) %>% dplyr::rename(eave_weight=ew)


n1 <- sum(EAVE_cohort_all$eave_weight)



### Remove children ####
n2 <- sum(df_cohort$eave_weight)

(n1-n2)

(n1-n2)/n1*100
(n2)/n1*100


### Remove elderly care-home residents ####
# Add in household information
Cohort_Household <- readRDS("/conf/EAVE/GPanalysis/outputs/temp/Cohort_Household.rds") %>%
  mutate(n_hh_gp = cut(n_hh, breaks=c(0,1,2,5,10,30,100,max(n_hh)),
                       labels=c("1", "2", "3-5", "6-10", "11-30", "31-100", "101+")))%>% 
  mutate(ave_hh_age=if_else(is.na(ave_hh_age), mean(ave_hh_age, na.rm=T), ave_hh_age) )

# Add in cohort household information
df_cohort <- df_cohort %>%
  left_join(select(Cohort_Household, EAVE_LINKNO,
                   n_hh_gp, ave_hh_age, care_home_elderly))

n2 <- sum(df_cohort$eave_weight[which(df_cohort$care_home_elderly==1)])
n2/n1*100


#### Eligible cohort ####
n2 <- sum(df_cohort$eave_weight)
sum(df_cohort$eave_weight)/sum(EAVE_cohort_all$eave_weight)



#### Vaccination information ####

# Add in vaccination data
df_cohort <- df_cohort %>%
  left_join(df_vaccinations) %>%
  mutate(dose_no = case_when((vacc_type %in% c("PB", "AZ") & vacc_type_2 %in% c("PB", "AZ")) ~ 2,
                             vacc_type %in% c("PB", "AZ") ~ 1,
                             TRUE ~ 0)) %>%
  mutate(v1_v2_days = date_vacc_2 - date_vacc_1) %>%
  mutate(vacc1 = if_else(vacc_type %in% c("PB", "AZ"), 1,0))  %>%
  mutate(vacc2 = if_else(vacc_type_2 %in% c("PB", "AZ"), 1,0))


## Vaccinated only
n3 <- sum(df_cohort$eave_weight[which(df_cohort$vacc1==1)])
n3/n2


#### Vaccine waning ####
### Non-care home cohort
df_cohort_nch <- df_cohort %>%
  filter(care_home_elderly == 0)

n4 <- sum(df_cohort_nch$eave_weight)

n4/n1*100

n5 <- sum(df_cohort_nch$eave_weight[which(df_cohort_nch$date_vacc_1<=as.Date("2021-04-30"))])
n5/n4

## Matched cohort
n6 <- nrow(df_cc_ps_matches)/2

# Percentages
n6/n5

n6/(n4-n5)

###### Failure cohort ####
df_cohort_vacc <- df_cohort %>%
  left_join(covid_hosp_death) %>%
  filter(vacc1 ==1) %>%
  mutate(time_to_event = as.numeric(date_vacc_1-admission_date)) #%>%
  #filter(is.na(time_to_event) | time_to_event > 14)

nrow(df_cohort_vacc)
sum(!is.na(df_cohort_vacc$time_to_event))

length(which(df_cohort_vacc$time_to_event < 15))





