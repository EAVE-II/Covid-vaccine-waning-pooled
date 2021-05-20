#### Summarising matching datasets #####

#### Set up ####
library("finalfit")

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
df_matches <- readRDS("/conf/EAVE/GPanalysis/progs/RM/Vaccine/Vaccine_waning/output/df_matches_death_hosp.rds")

df_cc_ps_matches <- readRDS("/conf/EAVE/GPanalysis/progs/RM/Vaccine/Vaccine_waning/output/df_cc_death_hosp.rds") %>%
  select(-c(event, time_to_hosp, time_to_event14, period)) %>%
  rename(event = event_hosp, time_to_hosp = time_to_hosp_hosp, time_to_event14 =time_to_event14_hosp,
         period = period_hosp)

z_title <- "COVID-19 hospitalisations"}

if (z_event_endpoint =="death_hosp") {z_event <- covid_hosp_death
df_matches <- readRDS("/conf/EAVE/GPanalysis/progs/RM/Vaccine/Vaccine_waning/output/df_matches_death_hosp.rds")
df_cc_ps_matches <- readRDS("/conf/EAVE/GPanalysis/progs/RM/Vaccine/Vaccine_waning/output/df_cc_death_hosp.rds")
z_title <- "COVID-19 hospitalisations or deaths"
}

if (z_event_endpoint =="positive_test") {z_event <- positive_test
#df_matches <- readRDS("/conf/EAVE/GPanalysis/progs/RM/Vaccine/Vaccine_waning/output/df_matches_positive_test.rds")
df_cc_ps_matches <- readRDS("/conf/EAVE/GPanalysis/progs/RM/Vaccine/Vaccine_waning/output/df_cc_positive_test.rds")
z_title <- "COVID-19 positive infections"
}

if (z_event_endpoint =="death_covid") {z_event <- covid_death
df_matches <- readRDS("/conf/EAVE/GPanalysis/progs/RM/Vaccine/Vaccine_waning/output/df_matches_death_hosp.rds")
df_cc_ps_matches <- readRDS("/conf/EAVE/GPanalysis/progs/RM/Vaccine/Vaccine_waning/output/df_cc_death_hosp.rds")%>%
  select(-c(event, time_to_hosp, time_to_event14, period)) %>%
  rename(event = event_death, time_to_hosp = time_to_hosp_death, time_to_event14 =time_to_event14_death,
         period = period_death)
z_title <- "COVID-19 deaths"
}


## Matches out out total eligible vaccinated ###
nrow(df_matches)/nrow(z_df)

nrow(df_matches)/nrow(z_chrt)

# No. matches which were vaccinated
df_matches %>%
  filter(!is.na(date_vacc_1_vacc)) %>%
  nrow()


# No. not matched
nrow(df_matches)-nrow(z_df)
1-nrow(df_matches)/nrow(z_df)
(nrow(df_matches)-nrow(z_df))/nrow(z_df)*100



#### Adding in characteristic information to cc data #####
qcovid_diags <- colnames(df_cohort)[startsWith(colnames(df_cohort), "Q")]

df_cc_desc <- df_cc_ps_matches %>%
  #left_join(any_hosp_nov_distinct, by=c("EAVE_LINKNO" = "EAVE_LINKNO")) %>%
  select(-eave_weight) %>%
  left_join(select(df_cohort, EAVE_LINKNO, Sex, test_before_dec8, EAVE_BP, EAVE_Smoke, HB, eave_weight,
                   qcovid_diags, bmi_cat),
            by=c("EAVE_LINKNO" = "EAVE_LINKNO")) %>% 
  left_join(select(df_cohort, EAVE_LINKNO, ageYear),
            by=c("EAVE_LINKNO_vacc" = "EAVE_LINKNO"), 
            suffix = c("", "_vacc")) %>% 
  left_join(covid_hosp_death %>%
              mutate(hosp_death = if_else(!is.na(NRS.Date.Death) & !is.na(hosp_admission_date), "both", outcome_date)) %>%
              select(EAVE_LINKNO, hosp_death), 
            by=c("EAVE_LINKNO" = "EAVE_LINKNO")) %>%
  mutate(hosp_death = replace_na(hosp_death, "no event")) %>%
  left_join(select(Cohort_Household, EAVE_LINKNO,
                   n_hh_gp, ave_hh_age, care_home_elderly),
            by= c("EAVE_LINKNO"= "EAVE_LINKNO")) %>%
  mutate(ur6_2016_name = replace_na(ur6_2016_name, "Unknown"))



#### Summary table - counts #####
explanatory <- c("event","event_hosp","event_death","Sex", "ageYear", "age_grp", "simd2020_sc_quintile", "ur6_2016_name", "n_risk_gps",
                 "n_tests_gp", "EAVE_Smoke", "bmi_cat", "ave_hh_age", "n_hh_gp")

# Vaccination status
dependent <- "vacc"
tbl4_tot <- df_cc_desc %>%
  summary_factorlist(dependent, explanatory, p = F)

head(tbl4_tot)


# Vaccination type
# AZ
tbl4_az <- df_cc_desc %>%
  filter(vacc_type == "AZ") %>%
  summary_factorlist(dependent, explanatory, p = F) %>%
  rename(uv_az = uv, vacc_az = vacc) 
head(tbl4_az)

# PB
tbl4_pb <- df_cc_desc %>%
  filter(vacc_type == "PB") %>%
  summary_factorlist(dependent, explanatory, p = F) %>%
  rename(uv_pb = uv, vacc_pb = vacc)
head(tbl4_pb)

# Merge
tbl4 <- tbl4_tot %>%
  left_join(left_join(tbl4_az, tbl4_pb))



write.csv(tbl4_tot, "./output/final/matching_summary/tbl4_tot.csv")
write.csv(tbl4_az, "./output/final/matching_summary/tbl4_az.csv")
write.csv(tbl4_pb, "./output/final/matching_summary/tbl4_pb.csv")

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





#### Rates using person years ####
# Overall
z.agg <- pyears(Surv(time_to_hosp,event) ~ vacc ,
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



#### Summary table - weights  #####
explanatory <- c("Sex", "ageYear", "age_grp", "simd2020_sc_quintile", "ur6_2016_name", "n_risk_gps",
                 "n_tests_gp", "test_before_dec8", "ave_hh_age", "n_hh_gp")

# Sum of weights for vacc and uv
sum(df_cc_desc$eave_weight[df_cc_desc$vacc == "vacc"])
length(which(df_cc_desc$vacc == "vacc"))/nrow(z_df)

sum(df_cc_desc$eave_weight[df_cc_desc$vacc == "uv"])

sum(df_cc_desc$eave_weight[df_cc_desc$vacc_type == "AZ"])
sum(df_cc_desc$eave_weight[df_cc_desc$vacc_type == "PB"])

sum(df_cc_desc$eave_weight)

sum(df_cc_desc$eave_weight[df_cc_desc$vacc == "vacc" & df_cc_desc$event_hosp==1])
sum(df_cc_desc$eave_weight[df_cc_desc$vacc == "uv" & df_cc_desc$event_hosp==1])

sum(df_cc_desc$eave_weight[df_cc_desc$vacc == "vacc" & df_cc_desc$event_death==1])
sum(df_cc_desc$eave_weight[df_cc_desc$vacc == "uv" & df_cc_desc$event_death==1])

sum(df_cc_desc$eave_weight[df_cc_desc$event_hosp==1 & df_cc_desc$vacc_type == "PB"])
sum(df_cc_desc$eave_weight[df_cc_desc$event_hosp==1 & df_cc_desc$vacc_type == "AZ"])

sum(df_cc_desc$eave_weight[df_cc_desc$event_death==1 & df_cc_desc$vacc_type == "PB"])
sum(df_cc_desc$eave_weight[df_cc_desc$event_death==1 & df_cc_desc$vacc_type == "AZ"])


# Individual characteristics
summary_tbl_wt <- summary_factorlist_wt(df_cc_desc,"vacc", explanatory = explanatory)
head(summary_tbl_wt)
write.csv(summary_tbl_wt, "./output/final/matching_summary/summary_table_weights.csv")

summary_tbl_wt <- summary_factorlist_wt(df_cc_desc,"vacc_type", explanatory = explanatory)
write.csv(summary_tbl_wt, "./output/final/matching_summary/summary_table_vacc_weights.csv")




##### Follow-up #####
df_cc_ps_matches <- df_cc_ps_matches %>%
  mutate(age_grp = case_when(ageYear < 65 ~"18-64", 
                             ageYear < 80 ~"65-79",
                             TRUE ~ "80+")) %>%
  left_join(select(df_cohort, EAVE_LINKNO, ageYear),
            by=c("EAVE_LINKNO_vacc" = "EAVE_LINKNO"), 
            suffix = c("", "_vacc")) %>%
  mutate(age_grp_vacc = case_when(ageYear_vacc < 65 ~"18-64", 
                                  ageYear_vacc < 80 ~"65-79",
                                  TRUE ~ "80+"))

# Overall time to follow up
ggplot(df_cc_ps_matches) +
  geom_density(aes(x=time_to_hosp), adjust=2) +
  geom_vline(xintercept = 14, linetype=2)

df_cc_ps_matches %>%
  summarise(median=median(time_to_hosp),
            IQR = IQR(time_to_hosp))

# Labels for vacc type
vacc_type_label <- c("BNT162b2", "ChAdOx1")
names(vacc_type_label) <- c("PB", "AZ")

# Split by vaccine
png(file=paste0("./output/final/matching_summary/followup_vacc.png"),
    width =600, height=400)
ggplot(df_cc_ps_matches) +
  #geom_density(aes(x=time_to_hosp, fill=vacc_type), alpha=0.5, adjust=3, stat="count")+
  geom_histogram(aes(x=time_to_hosp, fill=vacc_type), position = "dodge")+
  geom_vline(xintercept = 14, linetype=2) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green), labels = vacc_type_label) +
  labs(x="Follow-up (days)", fill="Vaccine type", 
       subtitle = paste0("Follow-up of ", z_title)) +
  annotate("text", x=14.5, y=400000, label = "14 days", hjust=0, size=3.5) +
  scale_x_continuous(breaks = seq(0,max(df_cc_ps_matches$time_to_hosp), by = 7))+
  scale_y_continuous(labels = function(x) format(x, scientific = F))

dev.off()

# Split by age group and vaccine 
png(file=paste0("./output/final/matching_summary/followup_vacc_age.png"),
    width =600, height=800)
df_cc_ps_matches %>%
  ggplot() +
  #geom_density(aes(x=time_to_hosp, fill=vacc_type), alpha=0.5, adjust=2) +
  geom_histogram(aes(x=time_to_hosp, fill=vacc_type), position = "dodge")+
  facet_wrap(~age_grp_vacc, ncol=1, scales="free")+
  geom_vline(xintercept = 14, linetype=2) +
  labs(x="Follow-up (days)", fill="Vaccine type",
       subtitle = paste0("Follow-up of ", z_title, " by age")) +
  #annotate("text", x=14.5, y=250000, label = "14 days", hjust=0, size=3.5)+
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green),labels = vacc_type_label)+
  scale_y_continuous(labels = function(x) format(x, scientific = F))

dev.off()


# Age, event and vaccine
png(file=paste0("./output/final/matching_summary/followup_vacc_age_event.png"),
    width =800, height=800)

ggplot(df_cc_ps_matches) +
  #geom_density(aes(x=time_to_hosp, fill=vacc_type), alpha=0.5, adjust=2)+
  geom_histogram(aes(x=time_to_hosp, fill=vacc_type), position = "dodge")+
  facet_grid(age_grp_vacc~event, scales="free")+
  geom_vline(xintercept = 14, linetype=2) +
  labs(x="Follow-up (days)", fill="Vaccine type", 
       subtitle = paste0("Follow-up of ", z_title, " by age and event")) +
  annotate("text", x=14.5, y=0.04, label = "14 days", hjust=0, size=3.5)+
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green))+
  scale_y_continuous(labels = function(x) format(x, scientific = F))
  

dev.off()



## Other plots

vaccine_status_label <- c("Control", "Exposed")
names(vaccine_status_label) <- c("uv", "vacc")

ggplot(df_cc_ps_matches) +
  geom_density(aes(x=time_to_hosp, fill=vacc_type), alpha=0.5, adjust=2)+
  geom_vline(xintercept = 14, linetype=2) +
  facet_grid(~vacc, labeller = labeller(vacc = vaccine_status_label)) +
  labs(x="Follow-up (days)", fill="Vaccine type") +
  annotate("text", x=14.5, y=0.04, label = "14 days", hjust=0, size=3.5)+
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green))

p1
summary(as.numeric(df_cc$time_to_hosp))

df_cc %>%
  group_by(vacc_type) %>%
  summarise(median=median(time_to_hosp),
            IQR = IQR(time_to_hosp))



ggplot(df_cc_ps_matches) +
  geom_density(aes(x=time_to_hosp, fill=vacc_type), alpha=0.5, adjust=2)+
  geom_vline(xintercept = 14, linetype=2) +
  facet_grid(age_2~vacc, labeller = labeller(vacc = vaccine_status_label)) +
  labs(x="Follow-up (days)", fill="Vaccine type", 
       subtitle = "Follow-up distributions by match group and age group") +
  annotate("text", x=14.5, y=0.06, label = "14 days", hjust=0, size=3.5)+
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green))



ggplot(df_cc_ps_matches) +
  geom_density(aes(x=time_to_hosp, fill=vacc), alpha=0.5, adjust=2)+
  geom_vline(xintercept = 14, linetype=2) +
  facet_grid(age_2~vacc_type) +
  labs(x="Follow-up (days)", fill="Vaccine type", 
       subtitle = "Follow-up distributions by vaccine and age group") +
  annotate("text", x=14.5, y=0.04, label = "14 days", hjust=0, size=3.5)+
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green))



df_cc_ps_matches %>%
  group_by(vacc_type, age_2) %>%
  summarise(median=median(time_to_hosp),
            IQR = IQR(time_to_hosp),
            n=n())

gridExtra::grid.arrange(p1, p2)

# By vaccine and event
ggplot(df_cc_ps_matches) +
  geom_density(aes(x=time_to_hosp, fill=vacc_type), alpha=0.5, adjust=2)+
  facet_wrap(~event, scales="free")+
  geom_vline(xintercept = 14, linetype=2)

### By vaccine, event and age


# All

ggplot(df_cc_ps_matches) +
  geom_density(aes(x=time_to_hosp, fill=vacc), alpha=0.5, adjust=2)+
  facet_grid(age_grp_vacc~event, scales="free")+
  geom_vline(xintercept = 14, linetype=2) +
  labs(x="Follow-up (days)", fill="Match group", 
       subtitle = "Follow-up distributions by vaccine and age group") +
  annotate("text", x=14.5, y=0.04, label = "14 days", hjust=0, size=3.5)+
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green))

# PB
p2 <- df_cc_ps_matches %>%
  filter(vacc_type == "PB") %>%
  ggplot() +
  geom_density(aes(x=time_to_hosp, fill=vacc_type), alpha=0.5, adjust=2)+
  facet_grid(age_grp_vacc~event, scales="free")+
  geom_vline(xintercept = 14, linetype=2) +
  labs(x="Follow-up (days)", fill="Match group", 
       subtitle = "PB: Follow-up distributions by vaccine and age group") +
  annotate("text", x=14.5, y=0.04, label = "14 days", hjust=0, size=3.5)+
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green), labels=c("Control", "Exposed"))

# AZ
p3 <- df_cc_ps_matches %>%
  filter(vacc_type == "AZ") %>%
  ggplot() +
  geom_density(aes(x=time_to_hosp, fill=vacc_type), alpha=0.5, adjust=2)+
  facet_grid(age_grp_vacc~event, scales="free")+
  geom_vline(xintercept = 14, linetype=2) +
  labs(x="Follow-up (days)", fill="Match group", 
       subtitle = "AZ: Follow-up distributions by vaccine and age group") +
  annotate("text", x=14.5, y=0.04, label = "14 days", hjust=0, size=3.5)+
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green), labels=c("Control", "Exposed"))

cowplot::plot_grid(p1, p2, p3, labels = "AUTO", ncol=1)


ggplot(df_cc) +
  geom_histogram(aes(x=time_to_hosp, fill=vacc_type), position = "dodge", bins=50) +
  facet_wrap(~age_gp_5yr)




###### Follow up for all outcomes ######
df_cc_ps_matches <- df_cc_ps_matches %>%
  mutate(age_grp = case_when(ageYear < 65 ~"18-64", 
                             ageYear < 80 ~"65-79",
                             TRUE ~ "80+")) %>%
  left_join(select(df_cohort, EAVE_LINKNO, ageYear),
            by=c("EAVE_LINKNO_vacc" = "EAVE_LINKNO"), 
            suffix = c("", "_vacc")) %>%
  mutate(age_grp_vacc = case_when(ageYear_vacc < 65 ~"18-64", 
                                  ageYear_vacc < 80 ~"65-79",
                                  TRUE ~ "80+"))

# Overall time to follow up
ggplot(df_cc_ps_matches) +
  geom_density(aes(x=time_to_hosp), adjust=2) +
  geom_vline(xintercept = 14, linetype=2)

df_cc_ps_matches %>%
  summarise(median=median(time_to_hosp),
            IQR = IQR(time_to_hosp))

# Split by vaccine
png(file=paste0("./output/final/matching_summary/followup_vacc.png"),
    width =600, height=400)
ggplot(df_cc_ps_matches) +
  geom_density(aes(x=time_to_hosp, fill=vacc_type), alpha=0.5, adjust=2)+
  geom_vline(xintercept = 14, linetype=2) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green)) +
  labs(x="Follow-up (days)", fill="Vaccine type", 
       subtitle = "Follow-up of exposed and control by vaccine type of exposed") +
  annotate("text", x=14.5, y=0.04, label = "14 days", hjust=0, size=3.5)

dev.off()

# Split by age group and vaccine 
png(file=paste0("./output/final/matching_summary/followup_vacc_age.png"),
    width =600, height=800)
df_cc_ps_matches %>%
  ggplot() +
  geom_density(aes(x=time_to_hosp, fill=vacc_type), alpha=0.5, adjust=2) +
  facet_wrap(~age_grp_vacc, ncol=1)+
  geom_vline(xintercept = 14, linetype=2) +
  labs(x="Follow-up (days)", fill="Vaccine type",
       subtitle = "Follow-up of exposed and control by vaccine type and age group of exposed") +
  annotate("text", x=14.5, y=0.04, label = "14 days", hjust=0, size=3.5)+
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green))

dev.off()


# Age, event and vaccine
png(file=paste0("./output/final/matching_summary/followup_vacc_age_event.png"),
    width =800, height=800)

ggplot(df_cc_ps_matches) +
  geom_density(aes(x=time_to_hosp, fill=vacc_type), alpha=0.5, adjust=2)+
  facet_grid(age_grp_vacc~event, scales="free")+
  geom_vline(xintercept = 14, linetype=2) +
  labs(x="Follow-up (days)", fill="Vaccine type", 
       subtitle = "Follow-up distributions by vaccine and age group") +
  annotate("text", x=14.5, y=0.04, label = "14 days", hjust=0, size=3.5)+
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green))

dev.off()








#### Person years summary #####
z.agg <- pyears(Surv(time_to_hosp,event) ~ vacc,
                data=df_cc_desc , scale=1, data.frame=TRUE, weights = eave_weight)

z.agg$data
table(df_cc_desc$vacc)

###### Propensity score matching summary ######

#df_cc_ps_matches <- readRDS("/conf/EAVE/GPanalysis/progs/RM/Vaccine/Vaccine_waning/output/df_cc_ps_matches.rds")




##### Overall distribution #####

ggplot(df_cc_ps_matches) +
  geom_density(aes(x=prop_score, fill=vacc), alpha=0.5, adjust=2)


### Check number of people used again
df_cc_ps_matches %>%
  group_by(EAVE_LINKNO) %>%
  summarise(n=n()) %>%
  arrange(desc(n)) %>%
  head()
# Max number of times someone was used was 7 times (model 1)
# For model 2 - most number of times used was 11 (2 people)
# Chris's update on 10th March - 2 people were used 15 times
# Final 8th April - max number of times someone was used was 24

df_cc_ps_matches %>%
  group_by(EAVE_LINKNO) %>%
  summarise(n=n()) %>%
  arrange(desc(n)) %>%
  ggplot()+
  geom_histogram(aes(x=n))



##### Covariate balance #####


cov_bal_vacc_fn <- function(data, explanatory, z_vacc_type){
  
  covariate_balance <- list()
  
  ## Find whether it is the matched or original dataset
  data_tbl <- table(data$vacc, data$vacc_type)
  
  # If summaries for uv by vacc type = 0 then this refers to the original dataset
  if(data_tbl[1,1]==0 &data_tbl[1,2]==0){
    
    # Replace the vacc_type for the uv as the z_vacc_type
    data <- data %>%
      filter(vacc_type == z_vacc_type | vacc == "uv") %>%
      mutate(vacc_type = z_vacc_type)
    
    
  } else {
    data <- data %>%
      filter(vacc_type == z_vacc_type)
  }
  
  
  # For each explanatory variable, get the standardised mean differences
  for(i in 1:length(explanatory)){
    
    # Identify whether binary or not
    n <- data %>%
      pull(!!sym(explanatory[i]))
    
    # If binary then choose a reference level and use this to compare between the other group
    if(length(unique(n)) == 2){
      
      level1 <- max(n)
      
      
      covariate_balance[[i]] <- data %>%
        mutate(var = ifelse(!!sym(as.character(explanatory)[i]) == level1, 1, 0)) %>%
        group_by(vacc) %>%
        summarise(weighted.mean = weighted.mean(var, w= eave_weight, na.rm = T),
                  var = spatstat.geom::weighted.var(var, w = eave_weight)) %>%
        mutate(mean.diff = weighted.mean[vacc=="uv"] -  weighted.mean[vacc=="vacc"],
               sd.pooled = sqrt((var[vacc=="uv"] + var[vacc=="vacc"])/2),
               smd = mean.diff/sd.pooled) %>%
        #select(smd) %>%
        mutate(characteristic = names(explanatory[i])) %>%
        mutate(levels = level1) %>%
        distinct()
      
      
    } else
      # If more than one group, take each level at a time and use as binary variable
      if(length(unique(n)) > 2){
        
        n_levels <- unique(n)
        
        
        covariate_balance_multi <- list()
        
        for(j in 1:length(n_levels)){
          level1 <- n_levels[j]
          
          if(!is.na(level1)) {
            covariate_balance_multi[[j]] <- data %>%
              mutate(var = ifelse(!!sym(explanatory[i]) == level1, 1, 0)) %>%
              group_by(vacc) %>%
              summarise(weighted.mean = weighted.mean(var, w= eave_weight, na.rm=T),
                        var = spatstat.geom::weighted.var(var, w = eave_weight, na.rm=T)) %>%
              mutate(mean.diff = weighted.mean[vacc=="uv"] -  weighted.mean[vacc=="vacc"],
                     sd.pooled = sqrt((var[vacc=="uv"] + var[vacc=="vacc"])/2),
                     smd = mean.diff/sd.pooled) %>%
              #select(smd) %>%
              mutate(characteristic = names(explanatory[i])) %>%
              mutate(levels = level1) %>%
              distinct() 
          } else {
            
            covariate_balance_multi[[j]] <- data %>%
              mutate(var = replace_na(!!sym(explanatory[i]), 1)) %>%
              mutate(var = ifelse(var == 1, 1, 0)) %>%
              group_by(vacc) %>%
              summarise(weighted.mean = weighted.mean(var, w= eave_weight, na.rm=T),
                        var = spatstat.geom::weighted.var(var, w = eave_weight, na.rm=T)) %>%
              mutate(mean.diff = weighted.mean[vacc=="uv"] -  weighted.mean[vacc=="vacc"],
                     sd.pooled = sqrt((var[vacc=="uv"] + var[vacc=="vacc"])/2),
                     smd = mean.diff/sd.pooled) %>%
              #select(smd) %>%
              mutate(characteristic = names(explanatory[i])) %>%
              mutate(levels = level1) %>%
              distinct() 
            
          }
          
          
          
        }
        
        covariate_balance [[i]] <- covariate_balance_multi %>%
          reduce(full_join)
        
        
      }
    
  }
  
  covariate_balance_all <- covariate_balance %>%
    reduce(full_join) %>%
    mutate(label = paste0(characteristic, ": ", levels)) %>%
    relocate(label) %>%
    arrange(desc(characteristic)) %>%
    mutate(vacc_type = z_vacc_type)
  
  covariate_balance_all
  
  
}

vacc_type_label <- c("B) BNT162b2", "A) ChAdOx1")
names(vacc_type_label) <- c("PB", "AZ")

explanatory <- c("Sex"="Sex", "Age (grouped)" = "age_grp", "Deprivation status" = "simd2020_sc_quintile", 
                 "Urban/Rural index" ="ur6_2016_name", "No. risk groups"="n_risk_gps",
                 "No. previous tests" = "n_tests_gp", "BMI" ="bmi_cat",
                 "Smoker status"="EAVE_Smoke", "BP status"="EAVE_BP")


cb_pb_crude <- cov_bal_vacc_fn(data = z_chrt_desc %>%
                                 mutate(vacc = recode(vacc1, "0" = "uv", "1" = "vacc")),
                               explanatory = explanatory, z_vacc_type = "PB")
cb_pb <- cov_bal_vacc_fn(data = df_cc_desc,
                         explanatory = explanatory, z_vacc_type = "PB")
cb_az_crude <- cov_bal_vacc_fn(data = z_chrt_desc %>%
                                 mutate(vacc = recode(vacc1, "0" = "uv", "1" = "vacc")),
                               explanatory = explanatory, z_vacc_type = "AZ")
cb_az <- cov_bal_vacc_fn(data = df_cc_desc,
                         explanatory = explanatory, z_vacc_type = "AZ")


cb_both <- full_join(cb_pb, cb_az) %>%
  mutate(data = "Matched") %>%
  full_join(full_join(cb_pb_crude, cb_az_crude) %>%
              mutate(data = "Crude"))

cb_both




write.csv(cb_both, "./output/final/matching_summary/covariate_balance_all.csv")



png(file=paste0("./output/final/matching_summary/covariate_balance.png"),
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


#### Risk factors 
explanatory <- qcovid_diags
names(explanatory) <- c("Atrial fibrillation", "Asthma", "Blood cancer", "Heart failure",
                        "Cerebalpalsy", "Coronary heart disease", "Cirrhosis", "Congenital  heart disease",
                        "COPD", "Dementia", "Diabetes type 1", "Diabetes type 2",
                        "Epilepsy", "Fracture", "Neurological disorder", "Parkinson's",
                        "Pulmonary hypertension", "Pulmonary rare", "Peripheral vascular disease", "Rheumatoid arthritis or SLE",
                        "Respiratory cancer", "Severe mental illness", "Sickle cell disease", "Stroke/TIA",
                        "Thrombosis or pulmonary embolus", "Care housing category", "Learning disability or Down's", "Kidney disease")

cb_pb_crude <- cov_bal_vacc_fn(data = z_chrt_desc %>%
                                 mutate(vacc = recode(vacc1, "0" = "uv", "1" = "vacc")),
                               explanatory = explanatory, z_vacc_type = "PB")
cb_pb <- cov_bal_vacc_fn(data = df_cc_desc,
                         explanatory = explanatory, z_vacc_type = "PB")
cb_az_crude <- cov_bal_vacc_fn(data = z_chrt_desc %>%
                                 mutate(vacc = recode(vacc1, "0" = "uv", "1" = "vacc")),
                               explanatory = explanatory, z_vacc_type = "AZ")
cb_az <- cov_bal_vacc_fn(data = df_cc_desc,
                         explanatory = explanatory, z_vacc_type = "AZ")


cb_both <- full_join(cb_pb, cb_az) %>%
  mutate(data = "Matched") %>%
  full_join(full_join(cb_pb_crude, cb_az_crude) %>%
              mutate(data = "Crude"))



write.csv(cb_both, "./output/final/matching_summary/covariate_balance_rg.csv")



png(file=paste0("./output/final/matching_summary/covariate_balance_rg.png"),
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



###### Covariate balance by Vaccinations - Propensity scores not %s ####

explanatory <- c("Sex"="Sex", "Age (grouped)" = "age_grp", "Deprivation status" = "simd2020_sc_quintile", 
                 "Urban/Rural index" ="ur6_2016_name", "No. risk groups"="n_risk_gps",
                 "No. previous tests" = "n_tests_gp", "Most recent test before 8th Dec" = "test_before_dec8",
                 "Smoker status"="EAVE_Smoke", "BP status"="EAVE_BP")


cov_bal_vacc_fn <- function(data, explanatory, z_vacc_type){
  covariate_balance <- list()
  
  for(i in 1:length(explanatory)){
    covariate_balance[[i]] <- data %>%
      filter(vacc_type == z_vacc_type) %>%
      select(vacc, !!sym(explanatory[i]), prop_score) %>%
      group_by(!!sym(explanatory[i])) %>%
      mutate(characteristic = (explanatory[i])) %>%
      mutate(characteristic_label = names(explanatory[i])) %>%
      mutate(mean_diff = mean(prop_score[vacc == "uv"]) - mean(prop_score[vacc == "vacc"])) %>%
      mutate(upr = t.test(prop_score[vacc == "uv"],
                          prop_score[vacc == "vacc"])$conf.int[1]) %>%
      mutate(lwr = t.test(prop_score[vacc == "uv"],
                          prop_score[vacc == "vacc"])$conf.int[2]) %>%
      mutate(p.value = t.test(prop_score[vacc == "uv"],
                              prop_score[vacc == "vacc"])$p.value) %>%
      select(-vacc, -prop_score) %>%
      distinct() %>%
      rename(levels = !!sym(explanatory[i])) %>%
      mutate(levels = as.character(levels)) %>%
      arrange(desc(levels)) %>%
      relocate(characteristic)
    
  }
  
  covariate_balance_all <- covariate_balance %>%
    reduce(full_join) %>%
    mutate(label = paste0(characteristic_label, ": ", levels)) %>%
    relocate(label) %>%
    arrange(desc(characteristic))
  
  covariate_balance_all
}

data <- z_chrt_desc %>%
  mutate(vacc = recode(vacc1, "0" = "uv", "1" = "vacc"))
z_vacc_type <- "PB"

data <- df_cc_desc

cb_pb <- cov_bal_vacc_fn(data = df_cc_desc, explanatory, "PB")
cb_pb_crude <- cov_bal_vacc_fn(data = z_chrt_desc %>%
                                 mutate(vacc = recode(vacc1, "0" = "uv", "1" = "vacc")), explanatory, "PB")

cb_az <- cov_bal_vacc_fn(data = df_cc_desc, explanatory, "AZ")





p1 <- ggplot(cb_pb, aes(x=mean_diff, y= label)) +
  geom_vline(xintercept = 0) + 
  geom_vline(xintercept = c(-0.1, 0.1), linetype=2) +
  geom_errorbar(aes(xmin=lwr, xmax=upr), width = 0.1) +
  geom_point(col=eave_blue, size=2) +
  theme_light() +
  labs(x="Mean difference", y="", subtitle = "BNT162b2")



p2 <- ggplot(cb_az, aes(x=mean_diff, y= label)) +
  geom_vline(xintercept = 0) + 
  geom_vline(xintercept = c(-0.1, 0.1), linetype=2) +
  geom_errorbar(aes(xmin=lwr, xmax=upr), width = 0.1) +
  geom_point(col=eave_blue, size=2) +
  theme_light() +
  labs(x="Mean difference", y="", subtitle = "ChAdOx1")


png(file=paste0("./output/final/matching_summary/covariate_balance_vacc.png"), width =800, height=600)

cowplot::plot_grid(p1, p2, ncol=2, labels = "AUTO")

dev.off()



check <- data %>%
  #filter(vacc_type == z_vacc_type) %>%
  select(vacc, !!sym(explanatory[i]), prop_score) %>%
  group_by(!!sym(explanatory[i])) %>%
  mutate(characteristic = (explanatory[i])) %>%
  mutate(characteristic_label = names(explanatory[i])) %>%
  mutate(mean_diff = mean(prop_score[vacc == "uv"]) - mean(prop_score[vacc == "vacc"]))


mean(check$prop_score[check$vacc == "uv" & check$Sex =="F"])
mean(check$prop_score[check$vacc == "vacc"& check$Sex =="F"])

##### Visualise ps of characteristics #####

ps_fn <- function(dta, variable, var_name) {
  
  
  # Subset to variable
  dta <- dta %>%
    select(variable, "vacc", "prop_score")
  
  # List
  list_outputs <- list()
  
  # Plot
  list_outputs[[1]] <- ggplot(dta, aes(x = prop_score, fill = vacc)) +
    geom_boxplot()+
    theme_light() +
    facet_wrap(variable) +
    scale_fill_manual(values = c(eave_blue, eave_green), labels = c("Control", "Exposed")) +
    labs(x="Propensity score", fill="Matched group", ylab="",
         subtitle = var_name) 
    
  # Plot
  list_outputs[[3]] <- ggplot(dta, aes(x = prop_score, fill = vacc)) +
    geom_density(adjust=2, alpha=0.5)+
    theme_light() +
    facet_wrap(variable) +
    scale_fill_manual(values = c(eave_blue, eave_green), labels = c("Control", "Exposed")) +
    labs(x="Propensity score", fill="Matched group", ylab="",
         subtitle = var_name) 
  
  # Table
  summary_tbl <- dta %>%
    group_by(!!sym(variable), vacc) %>%
    summarise(mean_ps = mean(prop_score),
              sd_ps = sd(prop_score))
  
  list_outputs[[2]] <- summary_tbl
  
  list_outputs
}


ps_fn(df_cc_ps_matches, "age_gp", "Age")
ps_fn(df_cc_ps_matches, "Sex", "Sex")
ps_fn(df_cc_ps_matches, "simd2020_sc_quintile", "Deprivation quintile")
ps_fn(df_cc_ps_matches, "n_tests_gp", "Number of previous COVID-19 tests")
ps_fn(df_cc_ps_matches, "n_risk_gps", "Number of risk groups")
ps_fn(df_cc_ps_matches, "ur6_2016_name", "Urban/rural index")
ps_fn(df_cc_ps_matches, "HB", "NHS Health Boards")


# Individual risk groups
qcovid_diags <- colnames(df_cohort)[startsWith(colnames(df_cohort), "Q")]
eave_diags <- c("EAVE_Smoke", "EAVE_BP")

df_cc_ps_matches <- df_cc_ps_matches %>%
  left_join(select(df_cohort, EAVE_LINKNO, qcovid_diags, eave_diags), by = c("EAVE_LINKNO" = "EAVE_LINKNO"))


ps_fn(df_cc_ps_matches, "Q_DIAG_ASTHMA")
ps_fn(df_cc_ps_matches, "Q_DIAG_BLOOD_CANCER")
ps_fn(df_cc_ps_matches, "Q_DIAG_CCF")
ps_fn(df_cc_ps_matches, "Q_DIAG_DEMENTIA")
ps_fn(df_cc_ps_matches, "Q_DIAG_PARKINSONS")
ps_fn(df_cc_ps_matches, "EAVE_Smoke")
ps_fn(df_cc_ps_matches, "EAVE_BP")



# All variables summary
variables <- c("age_2", "ageYear", "Sex", "simd2020_sc_quintile", "n_tests_gp", "n_risk_gps", "Council", "ur6_2016_name", qcovid_diags)
summary_tbl <- CreateTableOne(vars= variables, strata = "vacc", data=df_cc_ps_matches, smd=T)

CreateTableOne(vars= "ur6_2016_name", strata = "vacc", data=df_cc_ps_matches, smd=T)

write.csv(summary_tbl, "./outputs/tables/ps_summary_tbl.csv")



write.csv(summary_tbl2, "./outputs/tables/ps_summary_tbl.csv")



##### 65 yrs + matches ####

# Get cohort
df_cc_65 <- df_cc_desc %>%
  filter(ageYear_vacc > 64 & ageYear_vacc < 80)

colnames(df_cc_65)


df_cc_65 %>%
  summary_factorlist(dep, exp)



# PB:
df_cc_65PB <- df_cc_65 %>%
  filter(vacc_type == "PB")

df_cc_65PB %>%
  summary_factorlist(dep, exp)


# AZ:
df_cc_65AZ <- df_cc_65 %>%
  filter(vacc_type == "AZ")

colnames(df_cc_65)


df_cc_65AZ %>%
  summary_factorlist(dep, exp)


##### 80 yrs + matches ####
# Get cohort
df_cc_80 <- df_cc_desc %>%
  filter(ageYear_vacc >= 80)


df_cc_80 %>%
  summary_factorlist(dep, exp)






