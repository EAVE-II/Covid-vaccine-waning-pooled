##########################################################
## Title: 1st dose COVID-19 vaccine waning
## Code author(s): Rachel Mulholland <rachel.mulholland@ed.ac.uk> 
##                 Chris Robertson <chrisobertson@nhs.net>
## Description: 02_descriptive - Descriptive analyses on the
##              baseline cohort
##########################################################

#### 0 - Set up ####

# Libraries
library("finalfit")
library("lemon")

# Colours 
eave_green <- rgb(54, 176, 136, maxColorValue = 255)
eave_blue <- rgb(71,93,167, maxColorValue = 255)
eave_blue2 <- rgb(0,192,209, maxColorValue = 255)
eave_gold <- rgb(255,192,0, maxColorValue = 255)
eave_orange <- rgb(244,143,32, maxColorValue = 255)

# Outcome 
z_event_endpoint <- "death_hosp" # = death/hospitalisation

# Secondary outcomes:
# z_event_endpoint <- "covid_hospitalisations"
# z_event_endpoint <- "covid_death"

# Load in df_cohort and df_vaccinations data 
df_cohort <- readRDS("./output/df_cohort.rds")
df_vaccinations <- readRDS("./output/df_vaccinations.rds")

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
if (z_event_endpoint =="positive_test") {z_event <- positive_test}

# Find end date according to admission date
a_end <- max(z_event$admission_date)



##### 1 - Data set-up ####
# Link in vaccination and event data into baseline characteristics (df_cohort)
# to get descriptive dataset

## Vaccine data
# Subset vaccine data from 8th Dec until end date
z_vaccinations <- filter(df_vaccinations, date_vacc_1 <= a_end) %>% 
  mutate(vacc_type_2 = if_else(date_vacc_2 >= a_end, NA_character_ , vacc_type_2),
         date_vacc_2 = as.Date(ifelse(date_vacc_2 >= a_end, NA, date_vacc_2), origin=as.Date("1970-01-01")) )


## Baseline data
z_chrt_desc <- df_cohort %>%
  dplyr::rename(EAVE_LINKNO_uv = EAVE_LINKNO) %>%
  # Add NA UR as own group
  mutate(ur6_2016_name = replace_na(ur6_2016_name, "Unknown")) %>%
  # Add in age group
  mutate(age_grp = case_when(ageYear < 65 ~"18-64", 
                             ageYear < 80 ~"65-79",
                             TRUE ~ "80+"))

## Link in data to baseline data
# Add in vaccination data
z_chrt_desc <- z_chrt_desc %>%
  left_join(z_vaccinations, by=c("EAVE_LINKNO_uv" = "EAVE_LINKNO")) %>%
  mutate(dose_no = case_when((vacc_type %in% c("PB", "AZ") & vacc_type_2 %in% c("PB", "AZ")) ~ 2,
                             vacc_type %in% c("PB", "AZ") ~ 1,
                             TRUE ~ 0)) %>%
  mutate(v1_v2_days = date_vacc_2 - date_vacc_1) %>%
  mutate(vacc1 = if_else(vacc_type %in% c("PB", "AZ"), 1,0))  %>%
  mutate(vacc2 = if_else(vacc_type_2 %in% c("PB", "AZ"), 1,0))

# Link in event data
z_chrt_desc <- z_chrt_desc %>%
  left_join(select(z_event, -SpecimenDate), by=c("EAVE_LINKNO_uv" = "EAVE_LINKNO"))



# Make event indicators post vaccination and overall for uv
z_chrt_desc <- z_chrt_desc %>%
  # For composite outcome
  mutate(event = ifelse(admission_date <= date_vacc_1, "0", "1")) %>%
  mutate(event = ifelse(is.na(date_vacc_1) & !is.na(admission_date), "1", event)) %>%
  mutate(event = replace_na(event, "0")) %>%
  # For composite outcome, 14 days after
  mutate(event14 = ifelse(admission_date <= date_vacc_1+13, "0", "1")) %>%
  mutate(event14 = replace_na(event14, "0")) %>%
  # For hospitalisations
  mutate(event_hosp = ifelse(hosp_admission_date <= date_vacc_1,"0", "1")) %>%
  mutate(event_hosp = ifelse(is.na(date_vacc_1) & !is.na(hosp_admission_date), "1", event_hosp)) %>%
  mutate(event_hosp = replace_na(event_hosp, "0")) %>%
  # For deaths
  mutate(event_death = ifelse(NRS.Date.Death <= date_vacc_1,"0", "1")) %>%
  mutate(event_death = ifelse(is.na(date_vacc_1) & !is.na(NRS.Date.Death), "1", event_death)) %>%
  mutate(event_death = replace_na(event_death, "0"))



## Add in household information
# Household information (from Sept 2020)
Cohort_Household <- readRDS("/conf/EAVE/GPanalysis/outputs/temp/Cohort_Household.rds") %>%
  mutate(n_hh_gp = cut(n_hh, breaks=c(0,1,2,5,10,30,100,max(n_hh)),
                       labels=c("1", "2", "3-5", "6-10", "11-30", "31-100", "101+")))%>% 
  mutate(ave_hh_age=if_else(is.na(ave_hh_age), mean(ave_hh_age, na.rm=T), ave_hh_age) )

# Add in cohort household information
z_chrt_desc <- z_chrt_desc %>%
  left_join(select(Cohort_Household, EAVE_LINKNO,
                   n_hh_gp, ave_hh_age, care_home_elderly),
            by= c("EAVE_LINKNO_uv"= "EAVE_LINKNO"))


## In-hospital status pre-vaccinations
# If vacc then hospital status = 1 if admission date was at most 4 weeks before date of vacc
# If unvaccinated then flag if in-hospital throughout study (not used)
all_hospitalisations2 <- all_hospitalisations %>%
  left_join(z_vaccinations) %>%
  mutate(intervals = interval(start = admission_date, end = discharge_date)) %>%
  mutate(in_hosp_status = if_else(int_end(intervals) - 28 < date_vacc_1, 1, 0)) %>%
  mutate(in_hosp_status = if_else(is.na(date_vacc_1) & int_end(intervals) - 28 < a_end, 1, in_hosp_status)) %>%
  select(EAVE_LINKNO, in_hosp_status) %>%
  group_by(EAVE_LINKNO) %>%
  summarise(in_hosp_status = sum(in_hosp_status)) %>%
  mutate(in_hosp_status = if_else(in_hosp_status == 0 | is.na(in_hosp_status),0,1))

# Link into linked baseline data
z_chrt_desc <- z_chrt_desc %>%
  left_join(all_hospitalisations2, by= c("EAVE_LINKNO_uv"= "EAVE_LINKNO")) %>%
  mutate(in_hosp_status = replace_na(in_hosp_status, 0 )) %>%
  mutate(in_hosp_status = as.character(in_hosp_status))
  


## Previous positive test
# Indicator for whether person previously tested positive for COVID-19
# If vacc - tested positve before vacc
# If uv - if tested positive at all (not used)
z_chrt_desc <- z_chrt_desc %>%
  mutate(prev_positive_status = if_else(SpecimenDate < date_vacc_1, 1, 0)) %>%
  mutate(prev_positive_status = if_else(is.na(date_vacc_1) & SpecimenDate < a_end, 1, prev_positive_status))%>%
  mutate(prev_positive_status = replace_na(prev_positive_status, 0))%>%
  mutate(prev_positive_status = as.character(prev_positive_status))

# Quick summaries
prop.table(table(z_chrt_desc$in_hosp_status, z_chrt_desc$vacc1), 2)
prop.table(table(z_chrt_desc$prev_positive_status, z_chrt_desc$vacc1), 2)


## Total variable for summaries
z_chrt_desc <- z_chrt_desc %>%
  mutate(Total = "Total")

### Output dataset
saveRDS(z_chrt_desc, paste0("./output/z_chrt_desc.rds"))




##### 2 - Cohort summaries ####
## Using updated weights to get number of people

# Number of participants
n_tot <- sum(z_chrt_desc$eave_weight)
n_tot

# No. vacc
n <- sum(z_chrt_desc$eave_weight[z_chrt_desc$vacc1==1])
n
n/n_tot


## Split by vaccine type
n1 <- sum(z_chrt_desc$eave_weight[which(z_chrt_desc$vacc_type=="PB")])
n2 <- sum(z_chrt_desc$eave_weight[which(z_chrt_desc$vacc_type=="AZ")])
n1
n2
n1+n2 == n

n1/n

n2/n





##### 3 - Summary table (weights) ####
# Uses function summary_factorlist_wt from 00_functions.R

## Explanatory variables
# QCOVID risk groups
qcovid_diags <- colnames(df_cohort)[startsWith(colnames(df_cohort), "Q")]
# Other explanatory variables + qcovid_diag
explanatory <- c("Total","Sex", "ageYear", "age_grp", "simd2020_sc_quintile", "ur6_2016_name", "n_risk_gps",
                 "n_tests_gp", "test_before_dec8", "ave_hh_age", "n_hh_gp", "care_home_elderly",
                 "bmi_cat", "EAVE_Smoke", 
                 "event", qcovid_diags)

## Total population summary tables

# Dependent = vaccination status 
summary_tbl_wt1 <- summary_factorlist_wt(z_chrt_desc %>%
                                           mutate_at(vars(qcovid_diags), function(x) as.character(x)) %>% 
                                           mutate(care_home_elderly = as.character(care_home_elderly)),
                                         "vacc1", explanatory = explanatory)
# Dependent = vaccination type
summary_tbl_wt2 <- summary_factorlist_wt(z_chrt_desc %>%
                                           mutate_at(vars(qcovid_diags), function(x) as.character(x)) %>% 
                                           mutate(care_home_elderly = as.character(care_home_elderly)),
                                         "vacc_type", explanatory = explanatory)
# Combine and save as csv
summary_tbl_wt <- left_join(summary_tbl_wt1, summary_tbl_wt2)
write.csv(summary_tbl_wt, "./output/final/descriptives/summary_table_weights.csv", row.names = F)



## Non-elderly care home population summary tables
# Dependent = vaccination status 
summary_tbl_wt1 <- summary_factorlist_wt(z_chrt_desc %>% 
                                           mutate_at(vars(qcovid_diags), function(x) as.character(x)) %>% 
                                           filter(care_home_elderly == 0),
                                         "vacc1", explanatory = explanatory)
# Dependent = vaccination type
summary_tbl_wt2 <- summary_factorlist_wt(z_chrt_desc %>% 
                                           mutate_at(vars(qcovid_diags), function(x) as.character(x)) %>% 
                                           filter(care_home_elderly == 0),
                                         "vacc_type", explanatory = explanatory)
# Combine and save
summary_tbl_wt <- left_join(summary_tbl_wt1, summary_tbl_wt2)
write.csv(summary_tbl_wt, "./output/final/descriptives/summary_table_non_carehome_weights.csv", row.names=F)








#### 4 - Descriptive plots of total population #####
# Plots of vaccine uptake and events over time


### Vaccine uptake plots
# Histograms split by each vaccine type

# Labels for vacc type
vacc_type_label <- c("BNT162b2", "ChAdOx1")
names(vacc_type_label) <- c("PB", "AZ")


# Overall vaccine uptake
png(file=paste0("./output/final/descriptives/vacc_uptake.png"),
    width = 800, height=400)
z_chrt_desc %>%
ggplot() +
  geom_histogram(aes(x=date_vacc_1, fill=vacc_type), position="dodge", binwidth=2) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green), labels = vacc_type_label) +
  scale_color_manual(values = c(eave_blue, eave_green), labels = vacc_type_label) +
  labs(x="Date of 1st vaccine (2020/21)", fill="Vaccine type",
       subtitle = paste0("All 1st dose vaccine dates"),
       caption = paste0(a_begin, " to ", a_end)) +
  xlim(a_begin, a_end)
  
dev.off()


# Vaccine uptake by age group
png(file=paste0("./output/final/descriptives/vacc_uptake_age.png"),
    width = 800, height=800)
ggplot(z_chrt_desc) +
  geom_histogram(aes(x=date_vacc_1, fill=vacc_type), position="dodge", binwidth=2) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green), labels = vacc_type_label) +
  scale_color_manual(values = c(eave_blue, eave_green), labels = vacc_type_label) +
  labs(x="Date of 1st vaccine (2020/21)", fill="Vaccine type",
       subtitle = paste0("All 1st dose vaccine dates"),
       caption = paste0(a_begin, " to ", a_end)) +
  facet_wrap(~ age_grp, ncol=1) +
  xlim(a_begin, a_end)

dev.off()


## 2nd dose uptake
ggplot(z_chrt_desc) +
  geom_histogram(aes(x=date_vacc_2, fill=vacc_type), position="dodge", binwidth=2) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green)) +
  scale_color_manual(values = c(eave_blue, eave_green)) +
  labs(x="Date of 1st vaccine (2020/21)", fill="Vaccine type",
       subtitle = paste0("All 1st dose vaccine dates"),
       caption = paste0(a_begin, " to ", a_end)) +
  xlim(a_begin, a_end)




## Events

## COVID hospitalisations
p1 <- ggplot(z_chrt_desc) +
 # geom_density(aes(x=hosp_admission_date, fill=age_grp), alpha = 0.5) +
  geom_histogram(aes(x=hosp_admission_date, fill=age_grp), position="dodge", binwidth=5) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green, eave_orange)) +
  scale_color_manual(values = c(eave_blue, eave_green, eave_orange)) +
  labs(x="Date of COVID-19 hospital admission", fill="Age group",
       subtitle = "COVID-19 hospitalisation")



## COVID mortality
p2 <- ggplot(z_chrt_desc) +
  # geom_density(aes(x=hosp_admission_date, fill=age_grp), alpha = 0.5) +
  geom_histogram(aes(x=NRS.Date.Death, fill=age_grp), position="dodge", binwidth=5) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green, eave_orange)) +
  scale_color_manual(values = c(eave_blue, eave_green, eave_orange)) +
  labs(x="Date of COVID-19 death", fill="Age group",
       subtitle = "COVID-19 Death")


## COVID hospitalisations and deaths (Composite outcome)
p3 <- ggplot(z_chrt_desc) +
  geom_histogram(aes(x=admission_date, fill=age_grp), position="dodge", binwidth=5) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green, eave_orange)) +
  scale_color_manual(values = c(eave_blue, eave_green, eave_orange)) +
  labs(x="Date of COVID-19 hospitalisation or death", fill="Age group",
       subtitle = "Composite outcome: COVID-19 hospitalisation or death")

# COVID positive tests (not using)
p4 <- ggplot(z_chrt_desc) +
  geom_histogram(aes(x=SpecimenDate, fill=age_grp), position="dodge", binwidth=5) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green, eave_orange)) +
  scale_color_manual(values = c(eave_blue, eave_green, eave_orange)) +
  labs(x="Date of COVID-19 hospitalisation or death", fill="Age group",
       subtitle = "Positive COVID-19 test") +
  xlim(a_begin, a_end)


png(file=paste0("./output/final/descriptives/events_age.png"),
    width = 1000, height=300)
lemon::grid_arrange_shared_legend(p1, p2, p3, ncol = 3)
dev.off()







#### 5 - Descriptive plots of non-elderly care home population #####
# Same descriptives as 4 but for non-elderly care home residents

png(file=paste0("./output/final/descriptives/vacc_uptake_nonch.png"),
    width = 800, height=400)
z_chrt_desc %>%
  filter(care_home_elderly ==0) %>%
  ggplot() +
  geom_histogram(aes(x=date_vacc_1, fill=vacc_type), position="dodge", binwidth=2) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green), labels = vacc_type_label) +
  scale_color_manual(values = c(eave_blue, eave_green), labels = vacc_type_label) +
  labs(x="Date of 1st vaccine (2020/21)", fill="Vaccine type",
       subtitle = paste0("All 1st dose vaccine dates"),
       caption = paste0(a_begin, " to ", a_end, "
                        Excluding elderly care home residents"))

dev.off()


# Vaccinations over time by age
png(file=paste0("./output/final/descriptives/vacc_uptake_age_nonch.png"),
    width = 800, height=800)
z_chrt_desc %>%
  filter(care_home_elderly ==0) %>%
  ggplot() +
  geom_histogram(aes(x=date_vacc_1, fill=vacc_type), position="dodge", binwidth=2) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green), labels = vacc_type_label) +
  scale_color_manual(values = c(eave_blue, eave_green), labels = vacc_type_label) +
  labs(x="Date of 1st vaccine (2020/21)", fill="Vaccine type",
       subtitle = paste0("All 1st dose vaccine dates"),
       caption = paste0(a_begin, " to ", a_end, "
                        Excluding elderly care home residents")) +
  facet_wrap(~ age_grp, ncol=1)

dev.off()


## Event summary 

## COVID hospitalisations
p1 <- z_chrt_desc %>%
  filter(care_home_elderly ==0) %>%
  ggplot() +
  # geom_density(aes(x=hosp_admission_date, fill=age_grp), alpha = 0.5) +
  geom_histogram(aes(x=hosp_admission_date, fill=age_grp), position="dodge", binwidth=5) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green, eave_orange)) +
  scale_color_manual(values = c(eave_blue, eave_green, eave_orange)) +
  labs(x="Date of COVID-19 hospital admission", fill="Age group",
       subtitle = "COVID-19 hospitalisation")



## COVID mortality
p2 <- z_chrt_desc %>%
  filter(care_home_elderly ==0) %>%
  ggplot() +
  # geom_density(aes(x=hosp_admission_date, fill=age_grp), alpha = 0.5) +
  geom_histogram(aes(x=NRS.Date.Death, fill=age_grp), position="dodge", binwidth=5) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green, eave_orange)) +
  scale_color_manual(values = c(eave_blue, eave_green, eave_orange)) +
  labs(x="Date of COVID-19 death", fill="Age group",
       subtitle = "COVID-19 Death")


## Composite outcome
p3 <- z_chrt_desc %>%
  filter(care_home_elderly ==0) %>%
  ggplot() +
  geom_histogram(aes(x=admission_date, fill=age_grp), position="dodge", binwidth=5) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green, eave_orange)) +
  scale_color_manual(values = c(eave_blue, eave_green, eave_orange)) +
  labs(x="Date of COVID-19 hospitalisation or death", fill="Age group",
       subtitle = "Composite outcome: COVID-19 hospitalisation or death",
       caption = "Excluding elderly care home residents")

# save
png(file=paste0("./output/final/descriptives/events_age_nonch.png"),
    width = 1000, height=300)
lemon::grid_arrange_shared_legend(p1, p2, p3, ncol = 3)
dev.off()






##### 6 - Summary table (counts) - not used #####
# Provides same idea as 3 - Summary table (weights) but uses counts from 'finalfit' library

# Dependent = vaccination status 
z_chrt_desc %>%
  mutate(vacc1 = as.character(vacc1)) %>%
  summary_factorlist("vacc1", explanatory)

# Dependent = vaccination type 
z_chrt_desc %>%
  summary_factorlist(vacc_type, explanatory)

