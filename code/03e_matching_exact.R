##### Set up ####
# Options
output_list <- list()
# Exclude previous positive testing?
output_list$prev_pos <- "keep"
# What event end point file to look at?
z_event_endpoint <- "death" 
z_event_endpoint <- "hosp_covid"
z_event_endpoint <- "death_hosp"

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

# Chris's match:
z_cc_exact <- readRDS("/conf/EAVE/GPanalysis/progs/RM/Vaccine/Vaccine_waning/output/matched_vacc_unvacc_cc_31012021.RDS")

## Vaccination data
# Make anyone vaccinated after the maximum endpoint time unvaccinated
z_vaccinations <- filter(df_vaccinations, date_vacc_1 <= a_end) %>% 
  mutate(vacc_type_2 = if_else(date_vacc_2 >= a_end, NA_character_ , vacc_type_2),
         date_vacc_2 = as.Date(ifelse(date_vacc_2 >= a_end, NA, date_vacc_2), origin=as.Date("1970-01-01")) )


# Create event time
# Using the time of admission to hosp/icu/death to set up the event time
z_cc_exact <- z_cc_exact %>% 
  mutate(event_date = a_end) %>%
  # relocate the follow up to the date the unvaccinated control was vaccinated - Censoring?
  mutate(event_date = if_else(!is.na(date_vacc_1_uv) & (date_vacc_1_uv < event_date), 
                              date_vacc_1_uv, event_date))

# If vaccinated person becomes vaccinated then censor for this
z_cc_exact <- z_cc_exact %>%
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
z_cc_exact <- z_cc_exact %>% 
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
z_cc_exact <- z_cc_exact %>% 
  # Link death data
  left_join(any_death %>%
              rename("any_death_date" = "admission_date"), by="EAVE_LINKNO") %>%
  # Adjust event date if person dies before end date
  mutate(event_date = if_else(!is.na(any_death_date) & (any_death_date < event_date), 
                              any_death_date, event_date))%>%  
  dplyr::select(-SpecimenDate)


#change the event marker from 1 to 0 for those whose admission_date is greater than the current event_date
#these are instances where the unvaccinated control is vaccinated before the admission_date
#and hence censored at the vaccination date



## Calculate time until event
# Time until event is the event date - vaccination date
z_cc_exact <- z_cc_exact %>%  
  mutate(time_to_hosp = as.numeric(event_date-date_vacc_1_vacc))

##### Fix data errors #####
#time on study possible negative for data errors - omit both vacc and matched unvacc
# Find negative time to events (i.e. event happened before vaccination = error)
# Might be people vaccinated in hospital - Check
z_errors <- filter(z_cc_exact, time_to_hosp <0) # n =7
nrow(z_errors)

z_errors_ids <- unique(z_errors$EAVE_LINKNO_vacc)
length(z_errors_ids)


# Delete from dataset
z_cc_exact <- filter(z_cc_exact, 
               !(EAVE_LINKNO_vacc %in% z_errors_ids))


#link to vaccination type - link both vaccinated and unvaccinated control so that
#selection on vacc type can be easily made.
df_cc_exact <- z_cc_exact %>%
  left_join(dplyr::select(z_vaccinations,
                          EAVE_LINKNO, vacc_type), by=c("EAVE_LINKNO_vacc" = "EAVE_LINKNO"),
            suffix = c("", "_vacc")) %>%
  left_join(select(df_cohort, 
                   EAVE_LINKNO, age_gp, ageYear, simd2020_sc_quintile,
                   n_tests_gp, n_risk_gps, Council, ur6_2016_name), by=c("EAVE_LINKNO_vacc" = "EAVE_LINKNO"),
            suffix = c("", "_vacc"))%>%
  mutate(age_grp = case_when(ageYear < 65 ~"18-64", 
                             ageYear < 80 ~"65-79",
                             TRUE ~ "80+"))

# Recalculate time to event from 14 days
df_cc_exact <- df_cc_exact %>%
  mutate(time_to_event14 = as.numeric(event_date-(date_vacc_1_vacc+14)))

df_cc_exact$time_to_event14[which(df_cc_exact$time_to_event14 < 0)] <- NA



#### Put into time-periods
df_cc_exact <- df_cc_exact %>% 
  mutate(period = cut(time_to_hosp, 
                      breaks= c(-1, 6, 13, 20, 27, 34, 41, max(time_to_hosp, na.rm=T)),
                      labels=c("0:6","7:13","14:20","21:27","28:34","35:41", "42+")) ) 


#### Take out care home matches  ###
Cohort_Household <- readRDS("/conf/EAVE/GPanalysis/outputs/temp/Cohort_Household.rds") %>%
  mutate(n_hh_gp = cut(n_hh, breaks=c(0,1,2,5,10,30,100,max(n_hh)),
                       labels=c("1", "2", "3-5", "6-10", "11-30", "31-100", "101+")))

df_cc_exact <- df_cc_exact %>%
  left_join(select(Cohort_Household, EAVE_LINKNO, care_home_elderly),
            by=c("EAVE_LINKNO_vacc" = "EAVE_LINKNO"),
            suffix=c("", "_vacc")) %>%
  left_join(select(Cohort_Household, EAVE_LINKNO, care_home_elderly),
            by=c("EAVE_LINKNO" = "EAVE_LINKNO"),
            suffix=c("", "_uv"))


df_cc_exact_ch <- df_cc_exact %>%
  filter(care_home_elderly ==0 & care_home_elderly_uv == 0)


##### Cumulative plots  #####


#### Overall
par(mfrow=c(1,2))
z <- survfit(Surv(time_to_hosp, event) ~ vacc, data=df_cc_exact_ch)
plot(z, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, xlab="days from vaccination",ylab="cumulative risk")
legend("topleft",legend=levels(df_cc_exact$vacc), lty=1, col=c(1,2), cex=0.8)

# From 14 days
z <- survfit(Surv(time_to_event14, event) ~ vacc, data=df_cc_exact_ch)
plot(z, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk")
legend("topleft",legend=levels(df_cc_exact$vacc), lty=1, col=c(1,2), cex=0.8)


#### By age group
par(mfrow=c(3,2))
# 80+
z_age1 <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_exact_ch, subset=ageYear_vacc >= 80 & ageYear_vacc <= 110)
plot(z_age1, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "Age 80yrs +")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

z_age1 <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_exact_ch, subset=ageYear_vacc >= 80 & ageYear_vacc <= 110)
plot(z_age1, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "Age 80yrs +")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

# 65 to 79
z_age2 <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_exact_ch, subset=ageYear_vacc >= 65 & ageYear_vacc <= 79)
plot(z_age2, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "Age 65-79yrs")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

z_age2 <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_exact_ch, subset=ageYear_vacc >= 65 & ageYear_vacc <= 79)
plot(z_age2, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "Age 65-79yrs")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

# < 64
z_age3 <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_exact_ch, subset= ageYear_vacc < 65)
plot(z_age3, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "18-64yrs")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

z_age3 <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_exact_ch, subset= ageYear_vacc < 65)
plot(z_age3, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "18-64yrs")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)



### By Vaccine
par(mfrow=c(2,2))

# PB
z_pb <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_exact_ch, subset = (vacc_type=="PB"))
plot(z_pb, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "PB")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

z_pb <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_exact_ch, subset = (vacc_type=="PB"))
plot(z_pb, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "PB")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

# AZ
z_az <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_exact_ch, subset = (vacc_type=="AZ"))
plot(z_az, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "AZ")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

z_az <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_exact_ch, subset = (vacc_type=="AZ"))
plot(z_az, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "AZ")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)


#### Vaccine and age
par(mfrow=c(3,2))
## PB
# 80+
z_pb <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_exact_ch, 
                subset = (vacc_type=="PB") & ageYear >= 80 & ageYear <= 110)
plot(z_pb, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "PB and Age 80yrs +")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

z_pb <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_exact_ch, 
                subset = (vacc_type=="PB") & ageYear >= 80 & ageYear <= 110)
plot(z_pb, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "PB and Age 80yrs +")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

# 65-79
z_pb <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_exact_ch, 
                subset = (vacc_type=="PB") & ageYear >= 65 & ageYear <= 79)
plot(z_pb, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "PB and Age 65-79yrs")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

z_pb <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_exact_ch, 
                subset = (vacc_type=="PB") & ageYear >= 65 & ageYear <= 79)
plot(z_pb, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "PB and Age 65-79yrs")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

# <65
z_pb <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_exact_ch, 
                subset = (vacc_type=="PB") & ageYear < 65)
plot(z_pb, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "PB and Age 18-64yrs")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

z_pb <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_exact_ch, 
                subset = (vacc_type=="PB") & ageYear < 65)
plot(z_pb, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "PB and Age 18-64yrs")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)





## AZ
# 80+
z_az <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_exact_ch, 
                subset = (vacc_type=="AZ") & ageYear >= 80 & ageYear <= 110)
plot(z_az, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "AZ and Age 80yrs +")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

z_az <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_exact_ch, 
                subset = (vacc_type=="AZ") & ageYear >= 80 & ageYear <= 110)
plot(z_az, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "AZ and Age 80yrs +")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

# 65-79
z_az <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_exact_ch, 
                subset = (vacc_type=="AZ") & ageYear >= 65 & ageYear <= 79)
plot(z_az, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "AZ and Age 65-79yrs")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

z_az <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_exact_ch, 
                subset = (vacc_type=="AZ") & ageYear >= 65 & ageYear <= 79)
plot(z_az, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "AZ and Age 65-79yrs")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

# <65
z_az <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_exact_ch, 
                subset = (vacc_type=="AZ") & ageYear < 65)
plot(z_az, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "AZ and Age 18-64yrs")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)

z_az <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_exact_ch, 
                subset = (vacc_type=="AZ") & ageYear < 65)
plot(z_az, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "AZ and Age 18-64yrs")
legend("topleft",legend=levels(df_cc_exact_ch$vacc), lty=1, col=c(1,2), cex=0.8)





