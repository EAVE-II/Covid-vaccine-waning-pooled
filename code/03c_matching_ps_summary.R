##########################################################
## Title: 1st dose COVID-19 vaccine waning
## Code author(s): Rachel Mulholland <rachel.mulholland@ed.ac.uk> 
##                 Chris Robertson <chrisrobertson@nhs.net>
## Description: 03c_matching_ps_summary - Calculates summary stats
##              and plots for the matched cc cohort
##########################################################

#### 0 - Set up ####

# Libraries
library("finalfit")
library("survival")
library("survminer")

# Colours 
eave_green <- rgb(54, 176, 136, maxColorValue = 255)
eave_blue <- rgb(71,93,167, maxColorValue = 255)
eave_blue2 <- rgb(0,192,209, maxColorValue = 255)
eave_gold <- rgb(255,192,0, maxColorValue = 255)
eave_orange <- rgb(244,143,32, maxColorValue = 255)

# Event
z_event_endpoint <- "death_hosp"
#z_event_endpoint <- "positive_test"
#z_event_endpoint <- "hosp_covid"
#z_event_endpoint <- "death_covid"

### Load in data based on endpoint
if (z_event_endpoint =="hosp_covid") {z_event <- covid_hospitalisations
df_matches <- readRDS(paste0("./data/df_matches_first_dose_", multiplicity_limit, "_death_hosp.rds"))

df_cc_ps_matches <- readRDS(paste0("./data/df_cc_first_dose_", multiplicity_limit, "_death_hosp", ".rds")) %>%
  select(-c(event, time_to_event, time_to_event14, period)) %>%
  rename(event = event_hosp, time_to_event = time_to_hosp, time_to_event14 =time_to_event14_hosp,
         period = period_hosp)

z_title <- "COVID-19 hospitalisations"}

if (z_event_endpoint =="death_hosp") {z_event <- covid_hosp_death
df_matches <- readRDS(paste0("./data/df_matches_first_dose_", multiplicity_limit, "_death_hosp.rds"))

df_cc_ps_matches <- readRDS(paste0("./data/df_cc_first_dose_", multiplicity_limit, "_death_hosp", ".rds"))
z_title <- "COVID-19 hospitalisations or deaths"
}

if (z_event_endpoint =="positive_test") {z_event <- positive_test
df_cc_ps_matches <- readRDS(paste0("./data/df_cc_first_dose_", multiplicity_limit, "_", "positive_test.rds"))
z_title <- "COVID-19 positive infections"
}

if (z_event_endpoint =="death_covid") {z_event <- covid_death
df_matches <- readRDS(paste0("./data/df_matches_first_dose_", multiplicity_limit, "_death_hosp.rds"))

df_cc_ps_matches <- readRDS(paste0("./data/df_cc_first_dose_", multiplicity_limit, "_death_hosp", ".rds"))%>%
  select(-c(event, time_to_event, time_to_event14, period)) %>%
  rename(event = event_death, time_to_event = time_to_death, time_to_event14 =time_to_event14_death,
         period = period_death)
z_title <- "COVID-19 deaths"
}

# Find end date according to admission date
a_end <- as.Date("2021-06-30")

# Filter event data to end date
z_event <- z_event %>%
  filter(admission_date <= a_end)



#### 1 - Characteristics of cc cohort ####
# Adding in characteristic information to cc data

## Add in household information
# Household information (from Sept 2020)
Cohort_Household <- readRDS("/conf/EAVE/GPanalysis/outputs/temp/Cohort_Household.rds") %>%
  mutate(n_hh_gp = cut(n_hh, breaks=c(0,1,2,5,10,30,100,max(n_hh)),
                       labels=c("1", "2", "3-5", "6-10", "11-30", "31-100", "101+")))%>% 
  mutate(ave_hh_age=if_else(is.na(ave_hh_age), mean(ave_hh_age, na.rm=T), ave_hh_age) )

# QCOVID diagnoses 
qcovid_diags <- colnames(df_cohort)[startsWith(colnames(df_cohort), "Q")]

# Link characteristic info to df_cc_ps_matches
df_cc_desc <- df_cc_ps_matches %>%
  #select(-eave_weight) %>%
  # Baseline characteristics for everyone
  left_join(select(df_cohort, EAVE_LINKNO, Sex, test_before_dec8, EAVE_BP, EAVE_Smoke, HB, eave_weight,
                   qcovid_diags, bmi_cat),
            by=c("EAVE_LINKNO" = "EAVE_LINKNO")) %>% 
  # Baseline characteristics for vacc only
  left_join(select(df_cohort, EAVE_LINKNO, ageYear),
            by=c("EAVE_LINKNO_vacc" = "EAVE_LINKNO"), 
            suffix = c("", "_vacc")) %>% 
  # Event data for everyone
  left_join(covid_hosp_death %>%
              mutate(hosp_death = if_else(!is.na(NRS.Date.Death) & !is.na(hosp_admission_date), "both", outcome_date)) %>%
              select(EAVE_LINKNO, hosp_death), 
            by=c("EAVE_LINKNO" = "EAVE_LINKNO")) %>%
  # Household information
  left_join(select(Cohort_Household, EAVE_LINKNO,
                   n_hh_gp, ave_hh_age, care_home_elderly),
            by= c("EAVE_LINKNO"= "EAVE_LINKNO")) %>%
  # UR
  mutate(ur6_2016_name = replace_na(ur6_2016_name, "Unknown")) %>%
  # Add total column to get overall numbers
  mutate(Total = "total")



#### 2 - Summary table (counts) #####
# Use counts in matched cohort to demonstrate 1:1 matching ratio
# Totals
table(df_cc_desc$Total, df_cc_desc$vacc)
table(df_cc_desc$Total, df_cc_desc$vacc, df_cc_desc$vacc_type) # By vacc type

# Explanatory variables
explanatory <- c("event","Sex", "ageYear", "age_grp", 
                 "simd2020_sc_quintile", "ur6_2016_name", "n_risk_gps",
                 "n_tests_gp","test_before_dec8", "ave_hh_age", "n_hh_gp", "bmi_cat", "EAVE_Smoke",
                 qcovid_diags)

# Vaccination status
dependent <- "vacc"
tbl4_tot <- df_cc_desc %>%
  summary_factorlist(dependent, explanatory, add_col_totals = TRUE, p = F)


write.csv(tbl4_tot, paste0("./output/first_dose_", multiplicity_limit, "/final/matching_summary/tbl4_tot.csv"))


## Vaccination type
# AZ
tbl4_az <- df_cc_desc %>%
  filter(vacc_type == "AZ") %>%
  summary_factorlist(dependent, explanatory, add_col_totals = TRUE, p = F) %>%
  rename(uv_az = uv, vacc_az = vacc) 
head(tbl4_az)

write.csv(tbl4_az, paste0("./output/first_dose_", multiplicity_limit, "/final/matching_summary/tbl4_az.csv"))


# PB
tbl4_pb <- df_cc_desc %>%
  filter(vacc_type == "PB") %>%
  summary_factorlist(dependent, explanatory, add_col_totals = TRUE, p = F) %>%
  rename(uv_pb = uv, vacc_pb = vacc)
head(tbl4_pb)

write.csv(tbl4_pb, paste0("./output/first_dose_", multiplicity_limit, "/final/matching_summary/tbl4_pb.csv"))

# Bind tables together
#tbl4 <- bind_cols(tbl4_tot, select(tbl4_az, uv_az, vacc_az))
#tbl4 <- bind_cols(tbl4, select(tbl4_pb, uv_pb, vacc_pb))

tbl4 <- bind_cols(tbl4_az, select(tbl4_pb, uv_pb, vacc_pb))

write.csv(tbl4, paste0("./output/first_dose_", multiplicity_limit, "/final/matching_summary/tbl4.csv"))

## Total
df_cc_desc %>%
  group_by(vacc) %>%
  summarise(N = n())


## Medians
# Age year
df_cc_desc %>%
  group_by(vacc) %>%
  summarise(median = median(ageYear),
            iqr = IQR(ageYear))
df_cc_desc %>%
  group_by(vacc, vacc_type) %>%
  summarise(median = median(ageYear),
            iqr = IQR(ageYear))

# Average number
df_cc_desc %>%
  group_by(vacc) %>%
  summarise(median = median(ave_hh_age),
            iqr = IQR(ave_hh_age))
df_cc_desc %>%
  group_by(vacc, vacc_type) %>%
  summarise(median = median(ave_hh_age),
            iqr = IQR(ave_hh_age))





#### 3 - Rates of events using person years ####
## Total population
# Overall
z.agg <- pyears(Surv(time_to_event,event) ~ vacc,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)

df_res <- z.agg$data
df_res <- df_res %>% 
  mutate(pyears =round(pyears/365.25,1)) %>%
  mutate(Rate = round(event/pyears*1000,1)) %>% 
  #mutate(RR = Rate/first(Rate)) %>%
  select(-n) %>%
  mutate(label = paste0(event, " (", Rate,")"))
df_res

# Check
df_cc_desc %>%
  group_by(vacc) %>%
  summarise(event = sum(event))



## For 14 days only
z.agg <- pyears(Surv(time_to_event14,event) ~ vacc,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)

df_res <- z.agg$data
df_res <- df_res %>% 
  mutate(pyears =round(pyears/365.25,1)) %>%
  mutate(Rate = round(event/pyears*1000,1)) %>% 
  #mutate(RR = Rate/first(Rate)) %>%
  select(-n) %>%
  mutate(label = paste0(event, " (", Rate,")"))
df_res


## Vaccine type
# Overall
z.agg <- pyears(Surv(time_to_event,event) ~ vacc + vacc_type,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)

df_res <- z.agg$data
df_res <- df_res %>% 
  mutate(pyears =round(pyears/365.25,1)) %>%
  mutate(Rate = round(event/pyears*1000,1)) %>% 
  #mutate(RR = Rate/first(Rate)) %>%
  select(-n) %>%
  mutate(label = paste0(event, " (", Rate,")"))
df_res

# Check
df_cc_desc %>%
  group_by(vacc, vacc_type) %>%
  summarise(event = sum(event))



# For 14 days only
z.agg <- pyears(Surv(time_to_event14,event) ~ vacc + vacc_type,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)

df_res <- z.agg$data
df_res <- df_res %>% 
  mutate(pyears =round(pyears/365.25,1)) %>%
  mutate(Rate = round(event/pyears*1000,1)) %>% 
  #mutate(RR = Rate/first(Rate)) %>%
  select(-n) %>%
  mutate(label = paste0(event, " (", Rate,")"))
df_res




##### 4 - Follow-up #####

# Overall time to follow up
ggplot(df_cc_ps_matches) +
  geom_density(aes(x=time_to_event), adjust=2) +
  geom_vline(xintercept = 14, linetype=2)

df_cc_ps_matches %>%
  summarise(median=median(time_to_event),
            IQR = IQR(time_to_event))

# Labels for vacc type
vacc_type_label <- c("BNT162b2", "ChAdOx1")
names(vacc_type_label) <- c("PB", "AZ")

# Split by vaccine
png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/matching_summary/followup_vacc.png"),
    width =600, height=400)
ggplot(df_cc_ps_matches) +
  #geom_density(aes(x=time_to_event, fill=vacc_type), alpha=0.5, adjust=3, stat="count")+
  geom_histogram(aes(x=time_to_event, fill=vacc_type), position = "dodge")+
  geom_vline(xintercept = 14, linetype=2) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green), labels = vacc_type_label) +
  labs(x="Follow-up (days)", fill="Vaccine type", 
       subtitle = paste0("Follow-up of ", z_title)) +
  annotate("text", x=14.5, y=400000, label = "14 days", hjust=0, size=3.5) +
  scale_x_continuous(breaks = seq(0,max(df_cc_ps_matches$time_to_event), by = 7))+
  scale_y_continuous(labels = function(x) format(x, scientific = F))

dev.off()

# Split by age group and vaccine 
png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/matching_summary/followup_vacc_age.png"),
    width =600, height=800)
df_cc_ps_matches %>%
  ggplot() +
  #geom_density(aes(x=time_to_event, fill=vacc_type), alpha=0.5, adjust=2) +
  geom_histogram(aes(x=time_to_event, fill=vacc_type), position = "dodge")+
  facet_wrap(~age_grp, ncol=1, scales="free")+
  geom_vline(xintercept = 14, linetype=2) +
  labs(x="Follow-up (days)", fill="Vaccine type",
       subtitle = paste0("Follow-up of ", z_title, " by age")) +
  #annotate("text", x=14.5, y=250000, label = "14 days", hjust=0, size=3.5)+
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green),labels = vacc_type_label)+
  scale_y_continuous(labels = function(x) format(x, scientific = F))

dev.off()


# Age, event and vaccine
png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/matching_summary/followup_vacc_age_event.png"),
    width =800, height=800)

ggplot(df_cc_ps_matches) +
  #geom_density(aes(x=time_to_event, fill=vacc_type), alpha=0.5, adjust=2)+
  geom_histogram(aes(x=time_to_event, fill=vacc_type), position = "dodge")+
  facet_wrap(age_grp~event, scales="free", ncol=2)+
  geom_vline(xintercept = 14, linetype=2) +
  labs(x="Follow-up (days)", fill="Vaccine type", 
       subtitle = paste0("Follow-up of ", z_title, " by age and event")) +
  #annotate("text", x=14.5, y=0.04, label = "14 days", hjust=0, size=3.5)+
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green))+
  scale_y_continuous(labels = function(x) format(x, scientific = F))
  

dev.off()





##### 5 - Number of matches being used #####
# Checks how many matches were used repeatedly

z_chrt_desc <- readRDS('./data/z_chrt_desc.rds')

# Number of matches being used multiple times
match_multiplicity <-filter( df_cc_ps_matches, vacc == 'uv') %>%
  group_by(EAVE_LINKNO) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n = n())



png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/matching_summary/match_multiplicity_histogram.png"))

ggplot(match_multiplicity, aes(x= rownames(match_multiplicity), y= n )) + 
  geom_bar(stat="identity", fill= eave_green) + 
  labs(x = "Match multiplicity") + 
  scale_y_continuous(labels = scales::label_comma())

dev.off()

saveRDS(match_multiplicity, paste0("./output/first_dose_", multiplicity_limit, "/final/matching_summary/match_multiplicity_histogram.rds"))


# Number of vaccinated that get matched

# Get those who have been 1st dose vaccinated in the cohort time period
z_v1 <- filter(z_chrt_desc, date_vacc_1 <= a_end) 

n_v1 = nrow(z_v1)

n_matches = filter(df_cc_ps_matches, vacc == 'vacc') %>%
            select(EAVE_LINKNO) %>%
            nrow()

percent_matched <- n_matches/n_v1 * 100

match_stats <- data.frame('First dose vaccinated' =   format(n_v1 , nsmall=1, big.mark=","),
                             'Matched' =  paste0( format(n_matches , nsmall=1, big.mark=",") , ' (', 
                                                  round(percent_matched, 1), '%)'), 
                             check.names = FALSE )

write.csv(match_stats, paste0("./output/first_dose_", multiplicity_limit, "/final/matching_summary/match_stats.csv"))





##### Covariate balance #####
# Uses cov_bal_vacc_fn in 00_functions.R
# Requires z_chrt_desc from 02_descriptive.R

# Labels for vaccines
vacc_type_label <- c("B) BNT162b2", "A) ChAdOx1")
names(vacc_type_label) <- c("PB", "AZ")

### Main variables

### Save dataset
# Cohort descriptive dataframe required for covariate balance
z_chrt_desc <- readRDS("./data/z_chrt_desc.rds")



# Explanatory variables with labels as names
explanatory <- c("Sex"="Sex", "Age (grouped)" = "age_grp", "Deprivation status" = "simd2020_sc_quintile", 
                 "Urban/Rural index" ="ur6_2016_name", "No. risk groups"="n_risk_gps",
                 "No. previous tests" = "n_tests_gp", "BMI" ="bmi_cat",
                 "Smoker status"="EAVE_Smoke", "BP status"="EAVE_BP")

## PB
# Overall population (crude)
cb_pb_crude <- cov_bal_vacc_fn(data = z_chrt_desc %>%
                                 mutate(vacc = recode(vacc, "0" = "uv", "1" = "vacc")),
                               explanatory = explanatory, z_vacc_type = "PB")
# Matched population
cb_pb <- cov_bal_vacc_fn(data = df_cc_desc,
                         explanatory = explanatory, z_vacc_type = "PB")

## AZ
# Overall population (crude)
cb_az_crude <- cov_bal_vacc_fn(data = z_chrt_desc %>%
                                 mutate(vacc = recode(vacc, "0" = "uv", "1" = "vacc")),
                               explanatory = explanatory, z_vacc_type = "AZ")

# Matched population
cb_az <- cov_bal_vacc_fn(data = df_cc_desc,
                         explanatory = explanatory, z_vacc_type = "AZ")

## Combine tables
cb_both <- full_join(cb_pb, cb_az) %>%
  mutate(data = "Matched") %>%
  full_join(full_join(cb_pb_crude, cb_az_crude) %>%
              mutate(data = "Crude"))

cb_both

# Save
write.csv(cb_both, paste0("./output/first_dose_", multiplicity_limit, "/final/matching_summary/covariate_balance_all.csv"))


# Plot
png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/matching_summary/covariate_balance.png"),
    width =800, height=600)

ggplot(cb_both, aes(x=smd, y= label, shape = data, colour = data), size=3) +
  geom_vline(xintercept = 0) + 
  geom_vline(xintercept = c(-0.1, 0.1), linetype=2) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype=3) +
  #geom_errorbar(aes(xmin=lwr, xmax=upr), width = 0.1) +
  geom_point(size=2) +
  scale_color_manual(values = c(eave_blue, eave_green)) +
  theme_light() +
  facet_grid(~vacc_type, labeller= labeller(vacc_type = vacc_type_label)) + 
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black", hjust = 0, size=10)) +
  labs(x="Standardised mean differences", y="", title = "Covariate balance",
       shape = "Cohort", color = "Cohort",
       caption = "Positive standard mean differences suggests characteristic is more common in the unvaccinated group")


dev.off()


### Risk group variables (for supplementary)

explanatory <- qcovid_diags
names(explanatory) <- c("Atrial fibrillation", "Asthma", "Blood cancer", "Heart failure",
                        "Cerebalpalsy", "Coronary heart disease", "Cirrhosis", "Congenital  heart disease",
                        "COPD", "Dementia", "Diabetes type 1", "Diabetes type 2",
                        "Epilepsy", "Fracture", "Neurological disorder", "Parkinson's",
                        "Pulmonary hypertension", "Pulmonary rare", "Peripheral vascular disease", "Rheumatoid arthritis or SLE",
                        "Respiratory cancer", "Severe mental illness", "Sickle cell disease", "Stroke/TIA",
                        "Thrombosis or pulmonary embolus", "Care housing category", "Learning disability or Down's", "Kidney disease")
## PB
cb_pb_crude <- cov_bal_vacc_fn(data = z_chrt_desc %>%
                                 mutate(vacc = recode(vacc, "0" = "uv", "1" = "vacc")),
                               explanatory = explanatory, z_vacc_type = "PB")
cb_pb <- cov_bal_vacc_fn(data = df_cc_desc,
                         explanatory = explanatory, z_vacc_type = "PB")

## AZ
cb_az_crude <- cov_bal_vacc_fn(data = z_chrt_desc %>%
                                 mutate(vacc = recode(vacc, "0" = "uv", "1" = "vacc")),
                               explanatory = explanatory, z_vacc_type = "AZ")
cb_az <- cov_bal_vacc_fn(data = df_cc_desc,
                         explanatory = explanatory, z_vacc_type = "AZ")

# Join
cb_both <- full_join(cb_pb, cb_az) %>%
  mutate(data = "Matched") %>%
  full_join(full_join(cb_pb_crude, cb_az_crude) %>%
              mutate(data = "Crude"))


# Save
write.csv(cb_both, paste0("./output/first_dose_", multiplicity_limit, "/final/matching_summary/covariate_balance_rg.csv"))

# Plot
png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/matching_summary/covariate_balance_rg.png"),
    width =800, height=600)

ggplot(cb_both, aes(x=smd, y= label, shape = data, colour = data), size=3) +
  geom_vline(xintercept = 0) + 
  geom_vline(xintercept = c(-0.1, 0.1), linetype=2) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype=3) +
  #geom_errorbar(aes(xmin=lwr, xmax=upr), width = 0.1) +
  geom_point(size=2) +
  scale_color_manual(values = c(eave_blue, eave_green)) +
  theme_light() +
  facet_grid(~vacc_type, labeller= labeller(vacc_type = vacc_type_label)) + 
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black", hjust = 0, size=10)) +
  labs(x="Standardised mean differences", y="", title = "Covariate balance",
       shape = "Cohort", color = "Cohort",
       caption = "Positive standard mean differences suggests characteristic is more common in the unvaccinated group
Risk groups defined using QCOVID codes
Learning disability: 1 = Learning disability, 2 = Down's Syndrome
       Kidney disease: 1 = CKD3, 2 = CKD4, 3 = CKD5 without dialysis/transplant, 4 = CKD5 with dialysis, 5 = CKD5 with transplant
       Care housing category: 1 = Carehome, 2 = Homeless") 


dev.off()




