##### Checking outputs ####

#### Matches ####
z_merge <- df_matches_ps_m2_nonch 

# Checking if match has died before vaccinated match is vaccinated:
which(z_merge$any_death_date_uv < z_merge$date_vacc_1_vacc)

# Checking if match has been vaccinated before vaccinated match is vaccinated:
which(z_merge$date_vacc_1_uv < z_merge$date_vacc_1_vacc)

# Checking if match has event before vaccinated match is vaccinated:
which(z_merge$admission_date_uv < z_merge$date_vacc_1_vacc)

# Plot of prop_scores
colnames(z_merge)

ggplot(z_merge) +
  geom_point(aes(x=prop_score_vacc, y= prop_score_uv))




##### Match stacked data and time to event #####

## Event date is correct:
# minimum date of 1) Unvaccinated becomeing vaccinated, 2) Vaccinated becoming 2nd dose vacc, 3) event or 4) death
# If NA, this should be end date

df_cc_check <- df_cc %>%
  select(EAVE_LINKNO, EAVE_LINKNO_vacc, date_vacc_1_vacc, date_vacc_1_uv, date_vacc_2,
         admission_date, any_death_date, event_date) %>%
  group_by(EAVE_LINKNO) %>%
  mutate(event_date_test = min(date_vacc_1_uv, date_vacc_2, admission_date, any_death_date, na.rm = T))
  
df_cc_check <- df_cc_check %>%
  mutate(event_date_test2 = replace_na(event_date_test, a_end))






