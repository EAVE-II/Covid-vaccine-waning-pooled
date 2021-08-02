##########################################################
## Title: 1st dose COVID-19 vaccine waning
## Code author(s): Rachel Mulholland <rachel.mulholland@ed.ac.uk> 
## Description: 07_meta_analysis - Creates outputs for meta-analysis
##########################################################

#### 0 - Set up ####
# Run 0 - Set up in 03d_matching_modelling.R to set up df_cc_ps_matches to outcome


#### 1 - Table of daily follow-up #####

## All vaccines
z.agg <- pyears(Surv(time_to_hosp,event) ~ vacc + z.yr,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)

daily_gam_all <- z.agg$data %>%
  mutate(z.yr =  as.numeric(z.yr)-1)

daily_gam_all <- daily_gam_all %>%
  select(-pyears) %>%
  pivot_wider(names_from = vacc, values_from = c("n", "event")) %>%
  mutate(vacc_type = "All") %>%
  select(vacc_type, z.yr, n_uv, event_uv, n_vacc, event_vacc) 

head(daily_gam_all)


## Separate vaccines
z.agg <- pyears(Surv(time_to_hosp,event) ~ vacc + z.yr +vacc_type,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)

daily_gam_vacc <- z.agg$data %>%
  mutate(z.yr =  as.numeric(z.yr)-1)

daily_gam_vacc <- daily_gam_vacc %>%
  select(-pyears) %>%
  pivot_wider(names_from = vacc, values_from = c("n", "event")) %>%
  select(vacc_type, z.yr, n_uv, event_uv, n_vacc, event_vacc)



## Output
daily_gam_output <- bind_rows(daily_gam_all, daily_gam_vacc) %>%
  filter(z.yr <= 90)

write.csv(daily_gam_output, paste0("./output/meta-analysis/tbl_gam_", z_event_endpoint, ".csv"),
          row.names = F)
