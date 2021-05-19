##### Cohort analysis visualisations ####
cohort_analysis_tbl_12wks <- read_csv("output/final/cohort_analysis/cohort_analysis_tbl_12wks.csv") 
cohort_analysis_tbl_8wks <- read_csv("output/final/cohort_analysis/cohort_analysis_tbl_8wks.csv") 


head(cohort_analysis_tbl_8wks)
colnames(cohort_analysis_tbl_8wks)

# Replace high numbers by 3
cohort_analysis_tbl_8wks <- cohort_analysis_tbl_8wks %>%
  mutate_at(vars(matches("OR")), function(x) if_else(x > 3, 3, x))%>%
  filter(vacc_type != "Any") %>%
  filter(Vacc.Status != "uv")
  
lookup <- tibble(Vacc.Status = sort(unique(cohort_analysis_tbl_8wks$Vacc.Status)),
                 days = c(0, seq(14, 70, by=7)),
                 Vacc.Status2 = str_remove(Vacc.Status, "v1_"))



z_title <- "COVID-19 hospitalisations or deaths"

cohort_analysis_tbl_8wks <- cohort_analysis_tbl_8wks %>%
  left_join(lookup)

png(file=paste0("./output/final/cohort_analysis/cohort_analysis_RRs.png"),
    width =600, height=300)

cohort_analysis_tbl_8wks %>%
  ggplot() +
  # RR adjusted by age
  geom_point(aes(x=Vacc.Status2, y= OR_FULLADJ), col = eave_blue, size =2) +
  geom_errorbar(aes(ymin=OR_LCL_FULLADJ, ymax = OR_UCL_FULLADJ, x= Vacc.Status2), col = eave_blue) +
  facet_wrap(~vacc_type) +
  #geom_smooth(aes(x=days, y= OR_ADJ), method="gam", formula = y ~ s(x, bs = "tp"), col = eave_green, fill = eave_green) +
  #stat_smooth(aes(x=days, y= OR_ADJ), method ="lm", formula = y ~ x + I(x^2), col = eave_orange, fill = eave_orange) +
  labs(subtitle = paste0("Rate ratios of ", z_title, " based on cohort analysis"),
       y="Adjusted RR* (95% CI)", x = "Follow-up period",
       caption = "* Adjusted for time-period, age (groups of 5 years), sex, SIMD, individual risk groups, no. previous tests and Council Area.") + 
  theme_light() +
  ylim(0,1) +
  geom_hline(yintercept = 1, linetype=1) + 
  geom_hline(yintercept = 0.5, linetype = 3) +
  annotate("text", x=9, y=0.55, label = "Waning threshold", size=3)
  
dev.off()



cohort_analysis_tbl_8wks %>%
  ggplot() +
  # RR adjusted by age
  geom_point(aes(x=Vacc.Status2, y= OR_AGE), col = eave_blue, size =2) +
  geom_errorbar(aes(ymin=OR_LCL_AGE, ymax = OR_UCL_AGE, x= Vacc.Status2), col = eave_blue) +
  facet_wrap(~vacc_type) +
  #geom_smooth(aes(x=days, y= OR_ADJ), method="gam", formula = y ~ s(x, bs = "tp"), col = eave_green, fill = eave_green) +
  #stat_smooth(aes(x=days, y= OR_ADJ), method ="lm", formula = y ~ x + I(x^2), col = eave_orange, fill = eave_orange) +
  labs(subtitle = paste0("Rate ratios of ", z_title, " based on cohort analysis"),
       y="Adjusted RR* (95% CI)", x = "Follow-up period",
       caption = "* Adjusted for time-period, age (groups of 5 years), sex, SIMD, individual risk groups, no. previous tests and Council Area.") + 
  theme_light() +
  ylim(0,1) +
  geom_hline(yintercept = 1, linetype=1) + 
  geom_hline(yintercept = 0.5, linetype = 3) +
  annotate("text", x=9, y=0.55, label = "Waning threshold", size=3)

##### Checking summaries ####
z_out_results %>%
  filter(vacc_type == "PB" & Vacc.Status != "uv") %>%
  summarise(n=sum(Events))


z_pois %>%
  filter(vacc_status != "uv")%>%
  summarise(n=sum(event))

df_long_death_hosp_8wks %>%
  filter(vacc_status != "uv") %>%
  group_by(vacc_type) %>%
  summarise(n = count(EAVE_LINKNO))


#### 8 WEEKS #####
df_long_death_hosp_8wks <- df_long_death_hosp_8wks %>%
  left_join(select(Cohort_Household, EAVE_LINKNO,
                   n_hh_gp, ave_hh_age, care_home_elderly)) %>%
  filter(care_home_elderly == 0)

# No. people
length(unique(df_long_death_hosp_8wks$EAVE_LINKNO[which(df_long_death_hosp_8wks$vacc_status != "uv")]))


length(unique(df_long_death_hosp_8wks$EAVE_LINKNO[which(df_long_death_hosp_8wks$vacc_status != "uv" &
                                             df_long_death_hosp_8wks$vacc_type == "PB")]))


length(unique(df_long_death_hosp_8wks$EAVE_LINKNO[which(df_long_death_hosp_8wks$vacc_status != "uv")]))


# No. events
length(unique(df_long_death_hosp_8wks$EAVE_LINKNO[which(df_long_death_hosp_8wks$vacc_status != "uv" &
                                                          df_long_death_hosp_8wks$event == 1)]))
length(unique(df_long_death_hosp_8wks$EAVE_LINKNO[which(df_long_death_hosp_8wks$vacc_status != "uv" &
                                                          df_long_death_hosp_8wks$event == 1 &
                                                          df_long_death_hosp_8wks$vacc_type == "PB")]))
length(unique(df_long_death_hosp_8wks$EAVE_LINKNO[which(df_long_death_hosp_8wks$vacc_status != "uv" &
                                                          df_long_death_hosp_8wks$event == 1 &
                                                          df_long_death_hosp_8wks$vacc_type == "AZ")]))


