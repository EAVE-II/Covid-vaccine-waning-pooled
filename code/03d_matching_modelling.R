##########################################################
## Title: 1st dose COVID-19 vaccine waning
## Code author(s): Rachel Mulholland <rachel.mulholland@ed.ac.uk> 
##                 Chris Robertson <chrisrobertson@nhs.net>
## Description: 03d_matching_modelling - Carries out statistical
##              modelling (e.g. cumulative risks, GLMs, GAMs)
##########################################################

#### 0 - Set up ####

# Note: Paper uses sections 3, 6 and 7

# Libraries
library("survival")
library("survminer")
library("mgcv")
library("cowplot")

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

##### 1 - Cumulative risk incidence #####


#### Overall

png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/cumulative_risks/cumulative_overall.png"),
    width = 900, height=400)
par(mfrow=c(1,2))

z <- survfit(Surv(time_to_event, event) ~ vacc, data=df_cc_ps_matches)
plot(z, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, xlab="days from vaccination",ylab="cumulative risk")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

# From 14 days
z <- survfit(Surv(time_to_event14, event) ~ vacc, data=df_cc_ps_matches)
plot(z, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

dev.off()





### By Vaccine
png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/cumulative_risks/cumulative_vacc.png"),
    width = 900, height=600)


par(mfrow=c(2,2))

# PB
z_pb <- survfit(Surv(time_to_event, event) ~ vacc  , data=df_cc_ps_matches, subset = (vacc_type=="PB"))
plot(z_pb, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "PB")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

z_pb <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_ps_matches, subset = (vacc_type=="PB"))
plot(z_pb, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "PB")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

# AZ
z_az <- survfit(Surv(time_to_event, event) ~ vacc  , data=df_cc_ps_matches, subset = (vacc_type=="AZ"))
plot(z_az, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "AZ")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

z_az <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_ps_matches, subset = (vacc_type=="AZ"))
plot(z_az, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "AZ")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

dev.off()

#### Vaccine and age
png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/cumulative_risks/cumulative_pb_age.png"),
    width = 900, height=900)

par(mfrow=c(3,2))
## PB
# 80+
z_pb <- survfit(Surv(time_to_event, event) ~ vacc  , data=df_cc_ps_matches, 
                subset = (vacc_type=="PB") & ageYear >= 80 & ageYear <= 110)
plot(z_pb, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "PB and Age 80yrs +")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

z_pb <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_ps_matches, 
                subset = (vacc_type=="PB") & ageYear >= 80 & ageYear <= 110)
plot(z_pb, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "PB and Age 80yrs +")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

# 65-79
z_pb <- survfit(Surv(time_to_event, event) ~ vacc  , data=df_cc_ps_matches, 
                subset = (vacc_type=="PB") & ageYear >= 65 & ageYear <= 79)
plot(z_pb, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "PB and Age 65-79yrs")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

z_pb <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_ps_matches, 
                subset = (vacc_type=="PB") & ageYear >= 65 & ageYear <= 79)
plot(z_pb, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "PB and Age 65-79yrs")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

# <65
z_pb <- survfit(Surv(time_to_event, event) ~ vacc  , data=df_cc_ps_matches, 
                subset = (vacc_type=="PB") & ageYear < 65)
plot(z_pb, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "PB and Age 18-64yrs")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

z_pb <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_ps_matches, 
                subset = (vacc_type=="PB") & ageYear < 65)
plot(z_pb, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "PB and Age 18-64yrs")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

dev.off()



## AZ
png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/cumulative_risks/cumulative_az_age.png"),
    width = 900, height=900)
par(mfrow=c(3,2))
# 80+
z_az <- survfit(Surv(time_to_event, event) ~ vacc  , data=df_cc_ps_matches, 
                subset = (vacc_type=="AZ") & ageYear >= 80 & ageYear <= 110)
plot(z_az, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "AZ and Age 80yrs +")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

z_az <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_ps_matches, 
                subset = (vacc_type=="AZ") & ageYear >= 80 & ageYear <= 110)
plot(z_az, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "AZ and Age 80yrs +")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

# 65-79
z_az <- survfit(Surv(time_to_event, event) ~ vacc  , data=df_cc_ps_matches, 
                subset = (vacc_type=="AZ") & ageYear >= 65 & ageYear <= 79)
plot(z_az, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "AZ and Age 65-79yrs")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

z_az <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_ps_matches, 
                subset = (vacc_type=="AZ") & ageYear >= 65 & ageYear <= 79)
plot(z_az, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "AZ and Age 65-79yrs")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

# <65
z_az <- survfit(Surv(time_to_event, event) ~ vacc  , data=df_cc_ps_matches, 
                subset = (vacc_type=="AZ") & ageYear < 65)
plot(z_az, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "AZ and Age 18-64yrs")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

z_az <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_ps_matches, 
                subset = (vacc_type=="AZ") & ageYear < 65)
plot(z_az, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "AZ and Age 18-64yrs")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

dev.off()






##### 2 - Weekly GLM - Overall ####
# Performs poisson regression using person years of the vacc vs uv at weekly time periods

# Create weekly time periods of time to event
#z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,14,21,28,35,42, max(df_cc_ps_matches$time_to_event) ) )
z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,13,20,27,34,41, 48, 55, 62, 69, 76, 83, max(df_cc_ps_matches$time_to_event) ) )

# Calculate aggregated person years 
z.agg <- pyears(Surv(time_to_event,event) ~ vacc + z.yr,
                data=df_cc_ps_matches,scale=1,  data.frame=TRUE)

# Extract data
df_res <- z.agg$data

# Calculate rates and prepare for linkage with results
df_res <- df_res %>% 
  mutate(pyears =sprintf("%.2f",round(pyears/365.25,1))) %>%
  mutate(Rate = round( as.numeric(event)/ as.numeric(pyears)*1000,0)) %>% 
  mutate(vacc_z.yr = paste0("z.yr", z.yr, ":vacc",vacc)) %>%
  select(-n)

df_res

# Extract data for GLM
z_pois <- z.agg$data

# GLM model
z_glm <- glm(event ~ offset(log(pyears)) + -1 + z.yr + z.yr:vacc, 
             family=poisson, data=z_pois)
summary(z_glm)

# Obtain summary statistics and 95% CIs
z_glm_est <- as.data.frame(round(exp(cbind(z_glm$coefficients, confint.default(z_glm) ) ), 3)) %>%
  rownames_to_column(var ="vacc_z.yr") %>%
  rename(est=2, lwr = 3, upr = 4) %>%
  mutate(est = sprintf("%.2f",round(est,2)),
         upr = sprintf("%.2f",round(upr,2)),
         lwr = sprintf("%.2f",round(lwr,2))) %>%
  mutate(RR_est = paste0(est, " (", lwr, ", ", upr, ")")) %>%  
  select(-c(est, lwr, upr))

z_glm_est

# Combine GLM outputs with person years summary
z_glm_output <- df_res %>%
  left_join(z_glm_est, by = c('z.yr' = 'vacc_z.yr')) %>%
  select(-vacc_z.yr)

z_glm_output

# Save
write.csv(z_glm_output, paste0("./output/first_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/poisson_overall.csv"))



#### 3 - Weekly GLM - Vaccine Type ####
# Performs poisson regression using person years of the vacc vs uv at weekly time periods
# split by vaccine type

# Create weekly time periods of time to event
#z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,13,20,27,34,41, 48, 55, 62, 69, 76, 83, max(df_cc_ps_matches$time_to_event) ) )
z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,13,20,27,34,41, 48, 55, 62, 69, max(df_cc_ps_matches$time_to_event) ) )


# Calculate aggregated person years 
z.agg <- pyears(Surv(time_to_event,event) ~ vacc + z.yr +vacc_type,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)


# Extract data
df_res <- z.agg$data

# Calculate rates and prepare for linkage with results
df_res <- df_res %>% 
  mutate(pyears =round(pyears/365.25,1)) %>%
  mutate(Rate = round(event/pyears*1000)) %>% 
  #mutate(RR = Rate/first(Rate)) %>%
  mutate(id = paste0("z.yr", z.yr, ":vacc",vacc, ":vacc_type", vacc_type)) %>%
  select(-n)
df_res

# Check sum of events for vaccinated
sum(df_res$event[which(df_res$vacc == "vacc")])
sum(df_res$event[which(df_res$vacc == "uv")])

# Extract for GLM
z_pois <- z.agg$data

## Carry out GLM for each vaccine
## PB
# GLM
z_glm_pb <- glm(event ~ offset(log(pyears)) + -1 + z.yr + z.yr:vacc, 
             family=poisson, data=z_pois, subset= vacc_type == "PB")
summary(z_glm_pb)

# GLM outputs
z_glm_pb_est <- data.frame(round(exp(cbind(z_glm_pb$coefficients, confint.default(z_glm_pb) ) ), 3)) %>%
  rownames_to_column(var ="id") %>%
  mutate(id = paste0(id, ":vacc_typePB")) %>%
  rename(est=2, lwr = 3, upr = 4) %>%
  mutate(est = sprintf("%.2f",round(est,2)),
         upr = sprintf("%.2f",round(upr,2)),
         lwr = sprintf("%.2f",round(lwr,2))) %>%
  mutate(RR_est = paste0(est, " (", lwr, ", ", upr, ")")) %>%  
  select(-c(est, lwr, upr))


## AZ
# GLM
z_glm_az <- glm(event ~ offset(log(pyears)) + -1 + z.yr + z.yr:vacc, 
                family=poisson, data=z_pois, subset= vacc_type == "AZ")
summary(z_glm_az)

# GLM outputs
z_glm_az_est <- data.frame(round(exp(cbind(z_glm_az$coefficients, confint.default(z_glm_az) ) ), 3)) %>%
  rownames_to_column(var ="id") %>%
  mutate(id = paste0(id, ":vacc_typeAZ")) %>%
  rename(est=2, lwr = 3, upr = 4) %>%
  mutate(est = sprintf("%.2f",round(est,2)),
         upr = sprintf("%.2f",round(upr,2)),
         lwr = sprintf("%.2f",round(lwr,2))) %>%
  mutate(RR_est = paste0(est, " (", lwr, ", ", upr, ")")) %>%  
  select(-c(est, lwr, upr))


## Combine vaccine outputs
z_glm_vacc_output <- df_res %>%
  left_join(full_join(z_glm_pb_est,z_glm_az_est)) %>%
  select(-id)

# Save
write.csv(z_glm_vacc_output, 
          paste0("./output/first_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/poisson_vacc.csv"),
          row.names = F)




#### 4 - Weekly GLM - Vaccine Type by Age and Sex ####
# Performs poisson regression using person years of the vacc vs uv at weekly time periods
# split by vaccine type, age (18-64 and 65+) and Sex (female and male)

# Create weekly time periods of time to event
z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,13,20,27,34,41, 48, 55, 62, 69, max(df_cc_ps_matches$time_to_event) ) )

# Group age into 2 groups to allow for sufficient sample sizes
df_cc_ps_matches <- df_cc_ps_matches %>%
  mutate(age_grp2 = ifelse(age_grp == "18-64", "18-64", "65+")) %>%
  left_join(select(df_cohort, EAVE_LINKNO, Sex))



## Calculate weekly RR for age group and sex
# Use GLM_rr_var function in 00_functions.R

# Get unique levels
age_list <- unique(df_cc_ps_matches$age_grp2)
sex_list <- unique(df_cc_ps_matches$Sex)

# Assign RRs to list
z_glm_list <- list()   
    
z_glm_list[[1]] <- GLM_rr_var("Sex", sex_list[1])
z_glm_list[[2]] <- GLM_rr_var("Sex",sex_list[2])
z_glm_list[[3]] <- GLM_rr_var("age_grp2",age_list[1])
z_glm_list[[4]] <- GLM_rr_var("age_grp2",age_list[2])
    
# Create into data frame 
z_glm_outputs <- z_glm_list %>%
  reduce(full_join)
  
z_glm_outputs

# Save
write.csv(z_glm_outputs, paste0("./output/first_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/poisson_vacc_age_sex.csv"))


##### 5 - GAM - Overall ####

# time_to_event =z.yr as days rather than weeks
z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1:max(df_cc_ps_matches$time_to_event) ) )

#aggregate for overall
z.agg <- pyears(Surv(time_to_event,event) ~ vacc + z.yr,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)


z_pois <- z.agg$data


# GAM
z_gam <- gam(event ~ offset(log(pyears)) + -1 + vacc + s(as.numeric(z.yr), by=vacc), 
             family=poisson, data=z_pois)


# Predict risk of event
z_pois_pred <- z_pois
z_pred <- predict(z_gam, newdata=z_pois_pred, type="response", se =T)
z_pois_pred$est <- z_pred$fit
z_pois_pred$mu <- (z_pred$fit)^2
z_pois_pred$var.fit <- (z_pred$se.fit)^2
z_pois_pred$upr <- z_pred$fit + 1.96*z_pred$se.fit
z_pois_pred$lwr <- z_pred$fit - 1.96*z_pred$se.fit


## Calculate ratios
# Uses approximations for mean and variable of ratio (https://www.stat.cmu.edu/~hseltman/files/ratio.pdf)
# Did not use so can ignore

z_pois_pred_rr <- z_pois_pred %>%
  group_by(z.yr) %>%
  mutate(RR = est[vacc == "vacc"]/est[vacc == "uv"]) %>%
  mutate(var1 = ((est[vacc == "vacc"])^2)/((est[vacc == "uv"])^2)) %>%
  mutate(var2 = (var.fit[vacc == "vacc"])/((est[vacc == "vacc"])^2) + 
           (var.fit[vacc == "uv"])/((est[vacc == "uv"])^2)) %>%
  mutate(var = var1*var2) %>%
  select(-c(vacc, pyears, n, event, est, upr, lwr, var.fit, var1, var2)) %>%
  distinct() %>%
  mutate(upr = RR + 1.96*sqrt(var)) %>%
  mutate(lwr=RR + 1.96*sqrt(var))

z_pois_pred_rr

# Plot GAMs themselves
png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/gam_overall.png"),
    width = 700, height=400)

ggplot(z_pois_pred) +
  geom_line(aes(x=as.numeric(z.yr), y=est, col=vacc)) +
  geom_ribbon(aes(x=as.numeric(z.yr), ymin=lwr, ymax=upr, fill = vacc), alpha = 0.5) +
  scale_color_manual(values = c(eave_blue, eave_green), labels=c("Control", "Exposed"))+
  scale_fill_manual(values = c(eave_blue, eave_green), labels=c("Control", "Exposed")) +
  labs(x="Days to event", fill="Vaccination group",col="Vaccination group", 
       subtitle = "") +
  geom_vline(xintercept = 14, linetype=2) +
  #annotate("text", x=14.5, y=150, label = "14 days", hjust=0, size=3.5)+
  theme_light() 
  

dev.off()


# Plot GAM RRs
png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/gam_overall_rr.png"),
    width = 700, height=400)
ggplot(z_pois_pred_rr) +
  geom_line(aes(x=as.numeric(z.yr), y=RR), col=eave_blue) +
  geom_ribbon(aes(x=as.numeric(z.yr), ymin=lwr, ymax=upr), alpha = 0.25, fill=eave_blue) +
  #scale_color_manual(values = c(eave_blue, eave_green), labels=c("Control", "Exposed"))+
  #scale_fill_manual(values = c(eave_blue, eave_green), labels=c("Control", "Exposed")) +
  labs(x="Days to event", fill="Vaccination group", 
       subtitle = "") +
  geom_vline(xintercept = 14, linetype=2) +
  #annotate("text", x=14.5, y=150, label = "14 days", hjust=0, size=3.5)+
  theme_light() 
dev.off()




#### 6 - GAM and Cumulative plots by Vaccine type (Final plots used in paper) ####
# Final figures for publication by vacc type
# Was a function but didn't work with saving files automatically

# Select vaccine type
# PB has to be done first!
#z_vacc_type <- "PB"
z_vacc_type <- "AZ"

# Assign vaccine title
if(z_vacc_type == "PB"){
  z_vacc_title = "BNT162b2"
} else {
  z_vacc_title = "ChAdOx1"
}


### A: Overall cumulative plot
p_b_list <- list()
  
z <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_ps_matches, subset = (vacc_type==z_vacc_type))
p_b_list[[1]] <- ggsurvplot(z, data = df_cc_ps_matches, fun = "event", conf.int = T, 
           palette = c(eave_blue, eave_green),
           axes.offset = F,
           legend.labs = c("Control", "Recipient"),
           subtitle = "All ages",
           title = paste0("Cumulative risks of ",z_title," for ",z_vacc_title))


### B: Cumulative plot by age 

z <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_ps_matches, 
             subset = (vacc_type==z_vacc_type & ageYear_vacc >= 80 & ageYear_vacc <= 110))
p_b_list[[2]] <- ggsurvplot(z, data = df_cc_ps_matches, fun = "event", conf.int = T, 
                  palette = c(eave_blue, eave_green),
                  axes.offset = F,
                  legend.labs = c("Control", "Recipient"),
                  subtitle = "80yrs+")

z <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_ps_matches, 
             subset = (vacc_type==z_vacc_type & ageYear_vacc >= 65 & ageYear_vacc <= 79))
p_b_list[[3]] <- ggsurvplot(z, data = df_cc_ps_matches, fun = "event", conf.int = T, 
                   palette = c(eave_blue, eave_green),
                   axes.offset = F,
                   legend.labs = c("Control", "Recipient"),
                   subtitle = "65-79yrs")

z <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_ps_matches, 
             subset = (vacc_type==z_vacc_type & ageYear_vacc < 65))
p_b_list[[4]] <- ggsurvplot(z, data = df_cc_ps_matches, fun = "event", conf.int = T, 
                   palette = c(eave_blue, eave_green),
                   axes.offset = F,
                   legend.labs = c("Control", "Recipient"),
                   subtitle = "18-64yrs")

p_b <- arrange_ggsurvplots(p_b_list, ncol=1, nrow=4, print=T)


# Save
png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/plots/cumulative_risks_",z_vacc_type,".png"),
    width = 600, height=900)
p_b
dev.off()



### C: GAM
# Daily z.yr
z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,seq(0,max(df_cc_ps_matches$time_to_event),by=1) ))

# Person years
z.agg <- pyears(Surv(time_to_event,event) ~ vacc + z.yr +vacc_type,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)

# Calculate the number of events post 70 days (cut off due to low numbers)
z.agg$data %>% mutate(long = as.numeric(z.yr) >70) %>% 
  group_by(vacc, vacc_type, long) %>% summarise(N=sum(event))

# Take data and subset to vaccine type
z_pois <- z.agg$data %>%
  dplyr::mutate(day =  as.numeric(z.yr)-1) %>%
  dplyr::filter(vacc_type == z_vacc_type) 

# GAM
myknots <-  seq(14,70, by =7)

z_gam_vacc <- gam(event ~ offset(log(pyears)) + vacc + s(day, by=vacc, k=length(myknots)), 
                  knots = list(day = myknots), 
                  family=poisson, data=z_pois, subset= vacc_type == z_vacc_type)



## Calculate GAM RR
# Uses methods from https://fromthebottomoftheheap.net/2017/10/10/difference-splines-i/

z_nd <- z_pois  #for using for predictions

#this is where the modelling starts
z_p <- mgcv::predict.gam(z_gam_vacc, newdata=z_nd, type="lpmatrix")


#columns
c_uv <- grepl(":vaccuv", colnames(z_p)) #smoothcolumns for unvacc
c_vacc <- grepl(":vaccvacc", colnames(z_p))  #smooth columns for vacc

#rows
r_uv <- with(z_nd, vacc=="uv")
r_vacc <- with(z_nd, vacc=="vacc")

#difference
X <- z_p[r_vacc,] - z_p[r_uv,]

#zero other columns
X[,!(c_uv | c_vacc)] <- 0

#set vaccine column to 1 to get the effect
X[, "vaccvacc"] <- 1


difference <- X %*% coef(z_gam_vacc) 

# Obtain 95% CIs
se_difference <- sqrt(diag(X %*% vcov(z_gam_vacc, unconditional=TRUE) %*% t(X)))
z_lcl <- difference - 1.96*se_difference
z_ucl <- difference + 1.96*se_difference

# Create dataset with RRs and 95% CIs for plotting
z_rr <- tibble(rr = exp(difference),
               rr_lwr = exp(z_lcl),
               rr_upr = exp(z_ucl),
               day = 1:length(difference),
               rr_ci = paste0(round(rr,2), " (", round(rr_lwr,2), ", ", round(rr_upr,2), ")"))

# Join on no. events for vacc and uv
z_rr <- z_rr %>%
  left_join(filter(z_pois, vacc == "vacc") %>%
              select(event, day) %>%
              rename(event_vacc =event))
z_rr <- z_rr %>%
  left_join(filter(z_pois, vacc == "uv") %>%
              select(event, day) %>%
              rename(event_uv =event))



# Save table
write.csv(z_rr, paste0("./output/first_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/poisson_gam_", z_vacc_type, ".csv"))

## Waning p-value
# p-value for quadratic term post 14 days
z_quad <- glm(event ~ offset(log(pyears)) +  day*vacc + I(day^2)*vacc,
              family=poisson, data=z_pois, subset = vacc_type == z_vacc_type & day > 14)
z_linear <- glm(event ~ offset(log(pyears)) +  day*vacc,
              family=poisson, data=z_pois, subset = vacc_type == z_vacc_type & day > 14)


summary(z_quad)
# Check
zz_p <-predict.glm(z_quad, newdata=z_nd, type="response")
z_nd$fit <- zz_p
z_nd$prate <- z_nd$fit/z_nd$pyears
z_nd %>% ggplot(aes(x=day, y=prate, colour=vacc)) + geom_line()

# Get estimates and 95% CIs
z_quad_sum <- summary(z_quad)$coefficients
upr <- summary(z_quad)$coefficients[6,1] + 1.96*summary(z_quad)$coefficients[6,2]
lwr <- summary(z_quad)$coefficients[6,1] - 1.96*summary(z_quad)$coefficients[6,2]



## Plot RRs
p_c <- z_rr %>%
  mutate(rr_upr = if_else(rr_upr>3, 3, rr_upr)) %>%
  ggplot() +
  geom_line(aes(x=day, y=rr), col=eave_blue) +
  geom_ribbon(aes(x=day, ymin=rr_lwr, ymax=rr_upr), alpha = 0.25, fill=eave_blue) +
  labs(x="Days to event", y="RR",
       subtitle = paste0("Rate ratios of ", z_title," for ", z_vacc_title),
       caption = paste0("Vaccinated quadratic coefficient after 14 days: ",round(z_quad_sum[6,1],4),
                        #" (", round(upr,3), ", ",round(lwr,3),"), 
                        ", p.value = ", round(z_quad_sum[6,4],3))) +
  geom_vline(xintercept = 14, linetype=2) +
  annotate("text", x=14.5, y=2, label = "14 days", hjust=0, size=3.5)+
  theme_light() +
  lims(y=c(0,3)) +
  geom_hline(yintercept = 1, linetype = 1) +
#  geom_hline(yintercept = 0.5, linetype = 3) +
#  annotate("text", x=80, y=0.6, label = "Waning threshold", size=3) +
  scale_x_continuous(breaks = seq(14,84, by = 7), 
                     limits = c(0,84))
  

## Weekly RRs from z_rr to overlap ontop of p_c
z_rr_wkly <- z_rr %>%
  #mutate(rr_upr = ifelse(rr_upr==3, Inf, rr_upr)) %>%
  filter(day %in% seq(14,84, by = 14)) %>%
  mutate(rr_ci = paste0(round(rr,2), " (", round(rr_lwr,2), ", ", round(rr_upr,2), ")")) %>%
  mutate(y=1)

# Combine p_c and z_rr_wkly
p_c <- p_c +
  geom_text(aes(x=day, y=0, label = rr_ci), check_overlap = T, size=3, data=z_rr_wkly, inherit.aes = F)

p_c

# If PB then make copy so it can be plotted with AZ
if(z_vacc_type == "PB"){
  p_c1 = p_c
}

# Plot individual plot
png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/plots/gam_RR_",z_vacc_type,".png"),
    width = 800, height=400)

p_c
dev.off() 



## Plot GAMS together (p_c = AZ and p_c1 = PB)
png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/plots/gam_RR_both.png"),
    width = 800, height=600)

cowplot::plot_grid(p_c1,p_c, labels = "AUTO", ncol=1)

dev.off() 



##### 7 - GAM - Age and Sex ######
# Use function GAM_rr_var from 00_functions.R

df_cc_ps_matches <- df_cc_ps_matches %>%
  mutate(age_grp2 = ifelse(age_grp == "18-64", "18-64", "65+")) %>%
  left_join(select(df_cohort, EAVE_LINKNO, Sex))


## AZ
p_title <- ggdraw() +
  draw_label(paste0("Rate ratios of ", z_title," for ChAdOx1 by age and sex"),
             size=10)
p1 <- GAM_rr_var("AZ", "age_grp2","18-64")
p2 <- GAM_rr_var("AZ", "age_grp2","65+")
p3 <- GAM_rr_var("AZ", "Sex","F")
p4 <- GAM_rr_var("AZ", "Sex","M")


png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/plots/gam_RR_agesex_AZ.png"),
    width = 1000, height=600)

plot_grid(p_title, plot_grid(p1, p2, p3, p4, labels = "AUTO", ncol=2), ncol=1,
          rel_heights = c(0.5,5))

dev.off() 


## PB
p_title <- ggdraw() +
  draw_label(paste0("Rate ratios of ", z_title," for BNT162b2 by age and sex"),
             size=10)
p1 <- GAM_rr_var("PB", "age_grp2","18-64")
p2 <- GAM_rr_var("PB","age_grp2", "65+")
p3 <- GAM_rr_var("PB","Sex", "F")
p4 <- GAM_rr_var("PB","Sex", "M")


png(file=paste0("./output/first_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/plots/gam_RR_agesex_PB.png"),
    width = 1000, height=600)

plot_grid(p_title, plot_grid(p1, p2, p3, p4, labels = "AUTO", ncol=2), ncol=1,
          rel_heights = c(0.5,6))

dev.off() 


##### 8 - GAMS with different knots (Experimental) #####
z_vacc_type <- "PB"
#z_vacc_type <- "AZ"

if(z_vacc_type == "PB"){
  z_vacc_title = "BNT162b2"
} else {
  z_vacc_title = "ChAdOx1"
}

z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,seq(0,max(df_cc_ps_matches$time_to_event),by=1) ))


#aggregate for overall

z.agg <- pyears(Surv(time_to_event,event) ~ vacc + z.yr +vacc_type,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)

z.agg$data %>% mutate(long = as.numeric(z.yr) >70) %>% group_by(vacc, vacc_type, long) %>% summarise(N=sum(event))


z_pois <- z.agg$data %>%
  mutate(day =  as.numeric(z.yr)-1) %>%
  filter(vacc_type == z_vacc_type)

z_pois %>%
  filter(vacc == "vacc") %>%
  group_by(vacc_type) %>%
  summarise(events = sum(event))

myknots = seq(14,70, by =7)
z_gam_vacc <- gam(event ~ offset(log(pyears)) + vacc + s(day, by=vacc, k = length(myknots)), 
                  knots = list(day = myknots), 
                  family=poisson, data=z_pois, subset= vacc_type == z_vacc_type)

myknots <-  seq(14,70, by =7)

z_gam_vacc <- gam(event ~ offset(log(pyears)) + vacc + s(day, by=vacc, k=length(myknots)), 
                  knots = list(day = myknots), 
                  family=poisson, data=z_pois, subset= vacc_type == z_vacc_type)

my_knots <- z_gam_vacc$smooth[[1]]$xp
plot(x, y, col= "grey", main = "my knots");
lines(x, myfit$linear.predictors, col = 2, lwd = 2)
abline(v = my_knots, lty = 2)

par(mfrow=c(1,2))
plot(z_gam_vacc)

## Plot GAM RR

z_nd <- z_pois  #for using for predictions

#this is where the modelling starts
z_p <- mgcv::predict.gam(z_gam_vacc, newdata=z_nd, type="lpmatrix")


#columns
c_uv <- grepl(":vaccuv", colnames(z_p)) #smoothcolumns for unvacc
c_vacc <- grepl(":vaccvacc", colnames(z_p))  #smooth columns for vacc

#rows
r_uv <- with(z_nd, vacc=="uv")
r_vacc <- with(z_nd, vacc=="vacc")

#difference
X <- z_p[r_vacc,] - z_p[r_uv,]

#zero other columns
X[,!(c_uv | c_vacc)] <- 0

#set vaccine column to 1 to get the effect
X[, "vaccvacc"] <- 1

#offset difference
offset_diff <- log(z_nd[r_uv,"pyears"]) - log(z_nd[r_vacc,"pyears"])



difference <- X %*% coef(z_gam_vacc) + offset_diff

se_difference <- sqrt(rowSums((X %*% vcov(z_gam_vacc, unconditional=TRUE)) * X))
z_lcl <- difference - 1.96*se_difference
z_ucl <- difference + 1.96*se_difference


z_rr <- tibble(rr = exp(difference),
               rr_lwr = exp(z_lcl),
               rr_upr = exp(z_ucl),
               #day = 1:length(difference),
               day = min(z_pois$day):max(z_pois$day),
               rr_ci = paste0(round(rr,2), " (", round(rr_lwr,2), ", ", round(rr_upr,2), ")"))

z_rr <- z_rr %>%
  left_join(filter(z_pois, vacc == "vacc") %>%
              select(event, day) %>%
              rename(event_vacc =event))
z_rr <- z_rr %>%
  left_join(filter(z_pois, vacc == "uv") %>%
              select(event, day) %>%
              rename(event_uv =event))

#truncate for plot
#z_rr$rr_upr[z_rr$rr_upr>=3] <- 3

# Table
#write.csv(z_rr, paste0("./output/first_dose_", multiplicity_limit, "/final/modelling/", z_event_endpoint, "/poisson_gam_", z_vacc_type, ".csv"))

## P-value for quadratic term post 14 days
z_quad <- glm(event ~ offset(log(pyears)) +  day*vacc + I(day^2)*vacc, 
              family=poisson, data=z_pois, subset = vacc_type == z_vacc_type & day > 14)

zz_p <-predict.glm(z_quad, newdata=z_nd, type="response")
z_nd$fit <- zz_p
z_nd$prate <- z_nd$fit/z_nd$pyears
#z_nd %>% ggplot(aes(x=day, y=prate, colour=vacc)) + geom_line()

z_quad_sum <- summary(z_quad)$coefficients
upr <- summary(z_quad)$coefficients[6,1] + 1.96*summary(z_quad)$coefficients[6,2]
lwr <- summary(z_quad)$coefficients[6,1] - 1.96*summary(z_quad)$coefficients[6,2]



# Plot
p_c <- z_rr %>%
  mutate(rr_upr = if_else(rr_upr>3, 3, rr_upr)) %>%
  ggplot() +
  geom_line(aes(x=day, y=rr), col=eave_blue) +
  geom_ribbon(aes(x=day, ymin=rr_lwr, ymax=rr_upr), alpha = 0.25, fill=eave_blue) +
  labs(x="Days to event", y="RR",
       subtitle = paste0("Rate ratios of ", z_title," for ", z_vacc_title),
       caption = paste0("Vaccinated quadratic coefficient after 14 days: ",round(z_quad_sum[6,1],4),
                        #" (", round(upr,3), ", ",round(lwr,3),"), 
                        ", p.value = ", round(z_quad_sum[6,4],3))) +
  geom_vline(xintercept = 14, linetype=2) +
  annotate("text", x=14.5, y=2, label = "14 days", hjust=0, size=3.5)+
  theme_light() +
  lims(y=c(0,3)) +
  geom_hline(yintercept = 1, linetype = 1) +
#  geom_hline(yintercept = 0.5, linetype = 3) +
#  annotate("text", x=80, y=0.6, label = "Waning threshold", size=3) +
  scale_x_continuous(breaks = seq(14,84, by = 7), 
                     limits = c(0,84))


# Weekly z_rr
z_rr_wkly <- z_rr %>%
  #mutate(rr_upr = ifelse(rr_upr==3, Inf, rr_upr)) %>%
  filter(day %in% seq(14,84, by = 14)) %>%
  mutate(rr_ci = paste0(round(rr,2), " (", round(rr_lwr,2), ", ", round(rr_upr,2), ")")) %>%
  mutate(y=1)


p_c <- p_c +
  geom_text(aes(x=day, y=0, label = rr_ci), check_overlap = T, size=3, data=z_rr_wkly, inherit.aes = F)

p_c





