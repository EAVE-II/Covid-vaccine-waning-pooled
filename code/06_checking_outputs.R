##### Checking outputs ####

#### Matches ####
# Checking if match has died before vaccinated match is vaccinated:
which(df_matches$any_death_date_uv < df_matches$date_vacc_1_vacc)

# Checking if match has been vaccinated before vaccinated match is vaccinated:
which(df_matches$date_vacc_1_uv < df_matches$date_vacc_1_vacc)

# Checking if match has event before vaccinated match is vaccinated:
which(df_matches$admission_date_uv < df_matches$date_vacc_1_vacc)

# Plot of prop_scores
colnames(df_matches)

ggplot(df_matches) +
  geom_point(aes(x=prop_score_vacc, y= prop_score_uv))




##### Match stacked data and time to event #####

## Event date is correct:
# minimum date of 1) Unvaccinated becomeing vaccinated, 2) Vaccinated becoming 2nd dose vacc, 3) event or 4) death
# If NA, this should be end date

df_cc_ps_matches_check <- df_cc_ps_matches %>%
  select(EAVE_LINKNO, EAVE_LINKNO_vacc, date_vacc_1_vacc, date_vacc_1_uv, date_vacc_2,
         admission_date, any_death_date, event_date) %>%
  group_by(EAVE_LINKNO) %>%
  mutate(event_date_test = min(date_vacc_1_uv, date_vacc_2, admission_date, any_death_date, na.rm = T))
  
df_cc_ps_matches_check <- df_cc_ps_matches_check %>%
  mutate(event_date_test2 = replace_na(event_date_test, a_end))


#### Checking censoring worked ####
## Composite outcome
# How many people have their event date before the event?
z <- which(df_cc_ps_matches$event_date < df_cc_ps_matches$admission_date)
length(z)

sum(df_cc_ps_matches$event[z])

head(df_cc_ps_matches %>%
       filter(EAVE_LINKNO_vacc %in% df_cc_ps_matches$EAVE_LINKNO[z]))

## Hospitalisations
# How many people have their event date before the event?
z <- which(df_cc_ps_matches$event_date_hosp < df_cc_ps_matches$hosp_admission_date)
length(z)

sum(df_cc_ps_matches$event_hosp[z])


## Deaths
# How many people have their event date before the event?
z <- which(df_cc_ps_matches$event_date_death < df_cc_ps_matches$NRS.Date.Death)
length(z)

sum(df_cc_ps_matches$event_death[z])

#### Events ####
df_cc_ps_matches_events <- df_cc_ps_matches %>%
  filter(event == 1)
