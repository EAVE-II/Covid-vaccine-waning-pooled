##### Matching modelling #####
# Carrying out cumulative incidence plots, poisson regression and cox regression


### Cumulative incidence plots ####
library(survival)
library(survminer)


#### Overall

png(file=paste0("./output/final/modelling/", z_event_endpoint, "/cumulative_risks/cumulative_overall.png"),
    width = 900, height=400)
par(mfrow=c(1,2))

z <- survfit(Surv(time_to_hosp, event) ~ vacc, data=df_cc_ps_matches)
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

#### By age group
png(file=paste0("./output/final/modelling/", z_event_endpoint, "/cumulative_risks/cumulative_age.png"),
    width = 900, height=900)

par(mfrow=c(3,2))
# 80+
z_age1 <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_ps_matches, subset=ageYear_vacc >= 80 & ageYear_vacc <= 110)
plot(z_age1, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "Age 80yrs +")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

z_age1 <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_ps_matches, subset=ageYear_vacc >= 80 & ageYear_vacc <= 110)
plot(z_age1, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "Age 80yrs +")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

# 65 to 79
z_age2 <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_ps_matches, subset=ageYear_vacc >= 65 & ageYear_vacc <= 79)
plot(z_age2, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "Age 65-79yrs")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

z_age2 <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_ps_matches, subset=ageYear_vacc >= 65 & ageYear_vacc <= 79)
plot(z_age2, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "Age 65-79yrs")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

# < 64
z_age3 <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_ps_matches, subset= ageYear_vacc < 65)
plot(z_age3, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="days from vaccination",ylab="cumulative risk",
      sub = "18-64yrs")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

z_age3 <- survfit(Surv(time_to_event14, event) ~ vacc  , data=df_cc_ps_matches, subset= ageYear_vacc < 65)
plot(z_age3, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "18-64yrs")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

dev.off()



### By Vaccine
png(file=paste0("./output/final/modelling/", z_event_endpoint, "/cumulative_risks/cumulative_vacc.png"),
    width = 900, height=600)


par(mfrow=c(2,2))

# PB
z_pb <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_ps_matches, subset = (vacc_type=="PB"))
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
z_az <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_ps_matches, subset = (vacc_type=="AZ"))
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
png(file=paste0("./output/final/modelling/", z_event_endpoint, "/cumulative_risks/cumulative_pb_age.png"),
    width = 900, height=900)

par(mfrow=c(3,2))
## PB
# 80+
z_pb <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_ps_matches, 
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
z_pb <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_ps_matches, 
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
z_pb <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_ps_matches, 
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
png(file=paste0("./output/final/modelling/", z_event_endpoint, "/cumulative_risks/cumulative_az_age.png"),
    width = 900, height=900)
par(mfrow=c(3,2))
# 80+
z_az <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_ps_matches, 
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
z_az <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_ps_matches, 
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
z_az <- survfit(Surv(time_to_hosp, event) ~ vacc  , data=df_cc_ps_matches, 
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




### In hospital####
df_cc_ps_matches <- df_cc_ps_matches %>%
  left_join(select(any_hosp_nov_distinct, 
                   EAVE_LINKNO, inhosp_vacc1), by = c("EAVE_LINKNO_vacc" = "EAVE_LINKNO"),
            suffix=c("", "_vacc")) %>%
  mutate(inhosp_vacc1 = replace_na(inhosp_vacc1, 0))


## In hospital at time of vaccination
par(mfrow=c(2,2))

## In hosp
z_hosp <- survfit(Surv(time_to_hosp, event) ~ vacc, data=df_cc_ps_matches, 
                  subset = inhosp_vacc1 == 1)
plot(z_hosp, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, xlab="days from vaccination",ylab="cumulative risk",
      sub = "In hospital at time of vaccination")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

# From 14 days
z_hosp <- survfit(Surv(time_to_event14, event) ~ vacc, data=df_cc_ps_matches, 
                  subset = inhosp_vacc1 == 1)
plot(z_hosp, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "In hospital at time of vaccination")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

## Not in hosp
z_hosp <- survfit(Surv(time_to_hosp, event) ~ vacc, data=df_cc_ps_matches, 
                  subset = inhosp_vacc1 == 0)
plot(z_hosp, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, xlab="days from vaccination",ylab="cumulative risk",
      sub = "Not in hospital at vaccination")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

# From 14 days
z_hosp <- survfit(Surv(time_to_event14, event) ~ vacc, data=df_cc_ps_matches, 
                  subset = inhosp_vacc1 == 0)
plot(z_hosp, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "Not in hospital at time of vaccination")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)



### Hospital status by vaccine 

## PB 
z_hosp <- survfit(Surv(time_to_hosp, event) ~ vacc, data=df_cc_ps_matches, 
                  subset = inhosp_vacc1 == 1 & vacc_type == "PB")
plot(z_hosp, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, xlab="days from vaccination",ylab="cumulative risk",
      sub = "PB: In hospital at time of vaccination")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

# From 14 days
z_hosp <- survfit(Surv(time_to_event14, event) ~ vacc, data=df_cc_ps_matches, 
                  subset = inhosp_vacc1 == 1 & vacc_type == "PB")
plot(z_hosp, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "PB: In hospital at time of vaccination")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

# Not in hosp
z_hosp <- survfit(Surv(time_to_hosp, event) ~ vacc, data=df_cc_ps_matches, 
                  subset = inhosp_vacc1 == 0 & vacc_type == "PB")
plot(z_hosp, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, xlab="days from vaccination",ylab="cumulative risk",
      sub = "PB: Not in hospital at vaccination")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)

# From 14 days
z_hosp <- survfit(Surv(time_to_event14, event) ~ vacc, data=df_cc_ps_matches, 
                  subset = inhosp_vacc1 == 0 & vacc_type == "PB")
plot(z_hosp, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk",
      sub = "PB: Not in hospital at time of vaccination")
legend("topleft",legend=levels(df_cc_ps_matches$vacc), lty=1, col=c(1,2), cex=0.8)







##### GLM - Overall ####

### GLM
## Time to event - all
#z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,14,21,28,35,42, max(df_cc_ps_matches$time_to_hosp) ) )
z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,13,20,27,34,41, 48, 55, 62, 69, 76, 83, max(df_cc_ps_matches$time_to_hosp) ) )

#aggregate for overall
z.agg <- pyears(Surv(time_to_hosp,event) ~ vacc + z.yr,
                data=df_cc_ps_matches,scale=1,  data.frame=TRUE)



df_res <- z.agg$data
df_res <- df_res %>% 
  mutate(pyears =sprintf("%.2f",round(pyears/365.25,1))) %>%
  mutate(Rate = round(event/pyears*1000,0)) %>% 
  #mutate(RR = Rate/first(Rate)) %>%
  mutate(vacc_z.yr = paste0("z.yr", z.yr, ":vacc",vacc)) %>%
  select(-n)
df_res



z_pois <- z.agg$data

# Check without 0-13 days
z_pois <- z_pois %>%
  filter(z.yr != "-1+ thru  13")

z_glm <- glm(event ~ offset(log(pyears)) + -1 + z.yr + z.yr:vacc, 
             family=poisson, data=z_pois)
summary(z_glm)

z_glm_est <- as.data.frame(round(exp(cbind(z_glm$coefficients, confint.default(z_glm) ) ), 3)) %>%
  rownames_to_column(var ="vacc_z.yr") %>%
  rename(est=2, lwr = 3, upr = 4) %>%
  mutate(est = sprintf("%.2f",round(est,2)),
         upr = sprintf("%.2f",round(upr,2)),
         lwr = sprintf("%.2f",round(lwr,2))) %>%
  mutate(RR_est = paste0(est, " (", lwr, ", ", upr, ")")) %>%  
  select(-c(est, lwr, upr))

z_glm_est

z_glm_output <- df_res %>%
  left_join(z_glm_est) %>%
  select(-vacc_z.yr)

z_glm_output

write.csv(z_glm_output, paste0("./output/final/modelling/", z_event_endpoint, "/poisson_overall.csv"))



#### GLM - Vaccine Type ####
#z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,13,20,27,34,41, 48, 55, 62, 69, 76, 83, max(df_cc_ps_matches$time_to_hosp) ) )
z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,13,20,27,34,41, 48, 55, 62, 69, max(df_cc_ps_matches$time_to_hosp) ) )


#aggregate for overall
z.agg <- pyears(Surv(time_to_hosp,event) ~ vacc + z.yr +vacc_type,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)


df_res <- z.agg$data
df_res <- df_res %>% 
  mutate(pyears =round(pyears/365.25,1)) %>%
  mutate(Rate = round(event/pyears*1000)) %>% 
  #mutate(RR = Rate/first(Rate)) %>%
  mutate(id = paste0("z.yr", z.yr, ":vacc",vacc, ":vacc_type", vacc_type)) %>%
  select(-n)
df_res

sum(df_res$event[which(df_res$vacc == "vacc")])

z_pois <- z.agg$data



# PB
z_glm_pb <- glm(event ~ offset(log(pyears)) + -1 + z.yr + z.yr:vacc, 
             family=poisson, data=z_pois, subset= vacc_type == "PB")
summary(z_glm_pb)

z_glm_pb_est <- data.frame(round(exp(cbind(z_glm_pb$coefficients, confint.default(z_glm_pb) ) ), 3)) %>%
  rownames_to_column(var ="id") %>%
  mutate(id = paste0(id, ":vacc_typePB")) %>%
  rename(est=2, lwr = 3, upr = 4) %>%
  mutate(est = sprintf("%.2f",round(est,2)),
         upr = sprintf("%.2f",round(upr,2)),
         lwr = sprintf("%.2f",round(lwr,2))) %>%
  mutate(RR_est = paste0(est, " (", lwr, ", ", upr, ")")) %>%  
  select(-c(est, lwr, upr))



# AZ
z_glm_az <- glm(event ~ offset(log(pyears)) + -1 + z.yr + z.yr:vacc, 
                family=poisson, data=z_pois, subset= vacc_type == "AZ")
summary(z_glm_az)

z_glm_az_est <- data.frame(round(exp(cbind(z_glm_az$coefficients, confint.default(z_glm_az) ) ), 3)) %>%
  rownames_to_column(var ="id") %>%
  mutate(id = paste0(id, ":vacc_typeAZ")) %>%
  rename(est=2, lwr = 3, upr = 4) %>%
  mutate(est = sprintf("%.2f",round(est,2)),
         upr = sprintf("%.2f",round(upr,2)),
         lwr = sprintf("%.2f",round(lwr,2))) %>%
  mutate(RR_est = paste0(est, " (", lwr, ", ", upr, ")")) %>%  
  select(-c(est, lwr, upr))


# Combine
z_glm_vacc_output <- df_res %>%
  left_join(full_join(z_glm_pb_est,z_glm_az_est)) %>%
  select(-id)


write.csv(z_glm_vacc_output, 
          paste0("./output/final/modelling/", z_event_endpoint, "/poisson_vacc.csv"),
          row.names = F)




#### GLM - Vaccine Type by age and sex ####
z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,13,20,27,34,41, 48, 55, 62, 69, 76, 83, max(df_cc_ps_matches$time_to_hosp) ) )

df_cc_ps_matches <- df_cc_ps_matches %>%
  mutate(age_grp2 = ifelse(age_grp == "18-64", "18-64", "65+")) %>%
  left_join(select(df_cohort, EAVE_LINKNO, Sex))

table(df_cc_ps_matches$age_grp, df_cc_ps_matches$age_grp2)
sum(table(df_cc_ps_matches$Sex))/nrow(df_cc_ps_matches)
  

## Funtions
rr_sex <- function(sex){
  z.agg <- pyears(Surv(time_to_hosp,event) ~ vacc + z.yr +vacc_type +Sex,
                  data=df_cc_ps_matches , scale=1, data.frame=TRUE)
  
  df_res <- z.agg$data
  df_res <- df_res %>% 
    mutate(pyears =round(pyears/365.25,1)) %>%
    mutate(Rate = round(event/pyears*1000)) %>% 
    #mutate(RR = Rate/first(Rate)) %>%
    mutate(id = paste0("z.yr", z.yr, ":vacc",vacc, ":vacc_type", vacc_type)) %>%
    select(-n)
  df_res
  
  z_pois <- z.agg$data
  
  
  # PB
  z_glm_pb <- glm(event ~ offset(log(pyears)) + -1 + z.yr + z.yr:vacc, 
                  family=poisson, data=z_pois, subset= vacc_type == "PB" &
                    Sex == sex)
  summary(z_glm_pb)
  
  z_glm_pb_est <- data.frame(round(exp(cbind(z_glm_pb$coefficients, confint.default(z_glm_pb) ) ), 3)) %>%
    rownames_to_column(var ="id") %>%
    mutate(id = paste0(id, ":vacc_typePB")) %>%
    rename(est=2, lwr = 3, upr = 4) %>%   
    mutate(RR_est = paste0(round(est,2), " (", round(lwr,2), ", ", round(upr, 2), ")")) %>%   
    select(-c(est, lwr, upr))
  
  
  
  # AZ
  z_glm_az <- glm(event ~ offset(log(pyears)) + -1 + z.yr + z.yr:vacc, 
                  family=poisson, data=z_pois, subset= vacc_type == "AZ" &
                    Sex == sex)
  summary(z_glm_az)
  
  z_glm_az_est <- data.frame(round(exp(cbind(z_glm_az$coefficients, confint.default(z_glm_az) ) ), 3)) %>%
    rownames_to_column(var ="id") %>%
    mutate(id = paste0(id, ":vacc_typeAZ")) %>%
    rename(est=2, lwr = 3, upr = 4) %>%   
    mutate(RR_est = paste0(round(est,2), " (", round(lwr,2), ", ", round(upr, 2), ")")) %>%  
    select(-c(est, lwr, upr))
  
  
  # Combine
  z_glm_vacc_output <- df_res %>%
    filter(Sex == sex) %>%
    left_join(full_join(z_glm_pb_est,z_glm_az_est)) %>%
    rename(group = Sex) %>%
    select(-id)
  
  z_glm_vacc_output
  
  
}

rr_age <- function(age){
  z.agg <- pyears(Surv(time_to_hosp,event) ~ vacc + z.yr +vacc_type +age_grp2,
                  data=df_cc_ps_matches , scale=1, data.frame=TRUE)
  
  df_res <- z.agg$data
  df_res <- df_res %>% 
    mutate(pyears =round(pyears/365.25,1)) %>%
    mutate(Rate = round(event/pyears*1000)) %>% 
    #mutate(RR = Rate/first(Rate)) %>%
    mutate(id = paste0("z.yr", z.yr, ":vacc",vacc, ":vacc_type", vacc_type)) %>%
    select(-n)
  df_res
  
  z_pois <- z.agg$data
  
  # PB
  z_glm_pb <- glm(event ~ offset(log(pyears)) + -1 + z.yr + z.yr:vacc, 
                  family=poisson, data=z_pois, subset= vacc_type == "PB" &
                    age_grp2 == age)
  summary(z_glm_pb)
  
  z_glm_pb_est <- data.frame(round(exp(cbind(z_glm_pb$coefficients, confint.default(z_glm_pb) ) ), 3)) %>%
    rownames_to_column(var ="id") %>%
    mutate(id = paste0(id, ":vacc_typePB")) %>%
    rename(est=2, lwr = 3, upr = 4) %>%   
    mutate(RR_est = paste0(round(est,2), " (", round(lwr,2), ", ", round(upr, 2), ")")) %>%   
    select(-c(est, lwr, upr))
  
  
  
  # AZ
  z_glm_az <- glm(event ~ offset(log(pyears)) + -1 + z.yr + z.yr:vacc, 
                  family=poisson, data=z_pois, subset= vacc_type == "AZ" &
                    age_grp2 == age)
  summary(z_glm_az)
  
  z_glm_az_est <- data.frame(round(exp(cbind(z_glm_az$coefficients, confint.default(z_glm_az) ) ), 3)) %>%
    rownames_to_column(var ="id") %>%
    mutate(id = paste0(id, ":vacc_typeAZ")) %>%
    rename(est=2, lwr = 3, upr = 4) %>%   
    mutate(RR_est = paste0(round(est,2), " (", round(lwr,2), ", ", round(upr, 2), ")")) %>%  
    select(-c(est, lwr, upr))
  
  
  # Combine
  z_glm_vacc_output <- df_res %>%
    filter(age_grp2 == age) %>%
    left_join(full_join(z_glm_pb_est,z_glm_az_est)) %>%
    rename(group = age_grp2) %>%
    select(-id)
  
  z_glm_vacc_output
  
  
}
    
age_list <- unique(df_cc_ps_matches$age_grp2)
sex_list <- unique(df_cc_ps_matches$Sex)

z_glm_list <- list()   
    
z_glm_list[[1]] <- rr_sex(sex_list[1])
z_glm_list[[2]] <- rr_sex(sex_list[2])
z_glm_list[[3]] <- rr_age(age_list[1])
z_glm_list[[4]] <- rr_age(age_list[2])
    
z_glm_outputs <- z_glm_list %>%
  reduce(full_join)
  
z_glm_outputs



write.csv(z_glm_outputs, paste0("./output/final/modelling/", z_event_endpoint, "/poisson_vacc_age_sex.csv"))

write.csv(z_glm_outputs, paste0("./output/final/modelling/", z_event_endpoint, "/poisson_vacc_age.csv"))

##### GAM Overall ####
library(mgcv)

# Simple first

# time_to_hosp =z.year as days rather than weeks
z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1:max(df_cc_ps_matches$time_to_hosp) ) )

#aggregate for overall
z.agg <- pyears(Surv(time_to_hosp,event) ~ vacc + z.yr,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)


z_pois <- z.agg$data


# Alternative
z_gam <- gam(event ~ offset(log(pyears)) + -1 + vacc + s(as.numeric(z.yr), by=vacc), 
             family=poisson, data=z_pois)


par(mfrow=c(1,2))
plot(z_gam)



# Predict risk of event
z_pois_pred <- z_pois
z_pred <- predict(z_gam, newdata=z_pois_pred, type="response", se =T) # lpmatrix
z_pois_pred$est <- z_pred$fit
z_pois_pred$mu <- (z_pred$fit)^2
z_pois_pred$var.fit <- (z_pred$se.fit)^2
z_pois_pred$upr <- z_pred$fit + 1.96*z_pred$se.fit
z_pois_pred$lwr <- z_pred$fit - 1.96*z_pred$se.fit

# Ratios
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

z_pois_pred_rr <- z_pois_pred %>%
  group_by(z.yr) %>%
  mutate(RR = est[vacc == "vacc"]/est[vacc == "uv"]) %>%
  mutate(var1 = (mu[vacc == "vacc"])/(mu[vacc == "uv"])) %>%
  mutate(var2 = (var.fit[vacc == "vacc"])/(mu[vacc == "vacc"]) + 
           (var.fit[vacc == "uv"])/(mu[vacc == "uv"])) %>%
  mutate(var = var1*var2) %>%
  select(-c(vacc, pyears, n, event, est, upr, lwr, var.fit, var1, var2)) %>%
  distinct() %>%
  mutate(upr = RR + 1.96*sqrt(var)) %>%
  mutate(lwr=RR - 1.96*sqrt(var))

z_pois_pred_rr

## GAMS
png(file=paste0("./output/final/modelling/", z_event_endpoint, "/gam_overall.png"),
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


## GAM RR
png(file=paste0("./output/final/modelling/", z_event_endpoint, "/gam_overall_rr.png"),
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


###  GAM Vaccine type ####

# time_to_hosp =z.year as days rather than weeks
z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1:max(df_cc_ps_matches$time_to_hosp) ) )

#aggregate for overall
z.agg <- pyears(Surv(time_to_hosp,event) ~ vacc + z.yr +vacc_type,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)


z_pois <- z.agg$data


# PB
z_gam_pb <- gam(event ~ offset(log(pyears)) + -1 + vacc + s(as.numeric(z.yr), by=vacc), 
             family=poisson, data=z_pois, subset= vacc_type == "PB")


summary(z_gam_pb)

par(mfrow=c(1,3))
plot(z_gam_pb)


# AZ
z_gam_az <- gam(event ~ offset(log(pyears)) + s(as.numeric(z.yr), by=vacc)+ vacc, 
                family=poisson, data=z_pois, subset= vacc_type == "AZ")

par(mfrow=c(1,2))
plot(z_gam_az)



# Predict risk of event
z_pois_pred <- z_pois %>%
  mutate(id = paste0("z.yr", z.yr, ":vacc",vacc, ":vacc_type", vacc_type))

z_pred_pb <- predict(z_gam_pb, newdata=z_pois_pred %>%
                       filter(vacc_type == "PB"), type="response", se =T)

z_pred_pb_df <- tibble(est = z_pred_pb$fit,
                              upr = z_pred_pb$fit + 1.96*z_pred_pb$se.fit,
                              lwr = z_pred_pb$fit - 1.96*z_pred_pb$se.fit,
                       id = z_pois_pred %>%
                         filter(vacc_type == "PB") %>%
                         pull(id))

z_pred_az <- predict(z_gam_az, newdata=z_pois_pred %>%
                       filter(vacc_type == "AZ"), type="response", se =T)

z_pred_az_df <- tibble(est = z_pred_az$fit,
                       upr = z_pred_az$fit + 1.96*z_pred_az$se.fit,
                       lwr = z_pred_az$fit - 1.96*z_pred_az$se.fit,
                       id = z_pois_pred %>%
                         filter(vacc_type == "AZ") %>%
                         pull(id))

z_pred_vacc <- full_join(z_pred_pb_df, z_pred_az_df)

z_pois_pred <- z_pois_pred %>%
  left_join(z_pred_vacc)


png(file=paste0("./output/final/modelling/", z_event_endpoint, "/gam_vacc.png"),
    width = 700, height=400)

ggplot(z_pois_pred) +
  geom_line(aes(x=as.numeric(z.yr), y=est, col=vacc)) +
  geom_ribbon(aes(x=as.numeric(z.yr), ymin=lwr, ymax=upr, fill = vacc), alpha = 0.5) +
  facet_wrap(~vacc_type) + 
  scale_color_manual(values = c(eave_blue, eave_green), labels=c("Control", "Exposed"))+
  scale_fill_manual(values = c(eave_blue, eave_green), labels=c("Control", "Exposed")) +
  labs(x="Days to event", fill="Vaccination group",col="Vaccination group", 
       subtitle = "") +
  geom_vline(xintercept = 14, linetype=2) +
  annotate("text", x=14.5, y=150, label = "14 days", hjust=0, size=3.5)+
  theme_light() 

dev.off()



###### Final figures by vacc type - cumulative risks and GAMs #########

z_vacc_type <- "PB"
z_vacc_type <- "AZ"

if(z_vacc_type == "PB"){
  z_vacc_title = "BNT162b2"
} else {
  z_vacc_title = "ChAdOx1"
}

model_outputs_vacc <- function(z_vacc_type){
  
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



png(file=paste0("./output/final/modelling/", z_event_endpoint, "/plots/cumulative_risks_",z_vacc_type,".png"),
    width = 600, height=900)
p_b
dev.off()



### C: GAM
z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,seq(0,max(df_cc_ps_matches$time_to_hosp),by=1) ))


#aggregate for overall

z.agg <- pyears(Surv(time_to_hosp,event) ~ vacc + z.yr +vacc_type,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)

z.agg$data %>% mutate(long = as.numeric(z.yr) >70) %>% group_by(vacc, vacc_type, long) %>% summarise(N=sum(event))


z_pois <- z.agg$data %>%
  mutate(day =  as.numeric(z.yr)-1) %>%
  filter(vacc_type == z_vacc_type) #%>%
  #filter(day %in% 1:70)


z_gam_vacc <- gam(event ~ offset(log(pyears)) + vacc + s(day, by=vacc), 
                  knots = list(day = seq(14,70, by =7)), 
                  family=poisson, data=z_pois, subset= vacc_type == z_vacc_type)

# Compare linearity
#z_gam_vacc2 <- gam(event ~ offset(log(pyears)) + vacc + s(day), 
                 # family=poisson, data=z_pois, subset= vacc_type == z_vacc_type)


#AIC(z_gam_vacc, z_gam_vacc2)


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
               day = 1:length(difference),
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
write.csv(z_rr, paste0("./output/final/modelling/", z_event_endpoint, "/poisson_gam_", z_vacc_type, ".csv"))

## P-value for quadratic term post 14 days
z_quad <- glm(event ~ offset(log(pyears)) +  day*vacc + I(day^2)*vacc, 
              family=poisson, data=z_pois, subset = vacc_type == z_vacc_type & day > 14)
z_linear <- glm(event ~ offset(log(pyears)) +  day*vacc, 
              family=poisson, data=z_pois, subset = vacc_type == z_vacc_type & day > 14)
summary(z_quad)
zz_p <-predict.glm(z_quad, newdata=z_nd, type="response")
z_nd$fit <- zz_p
z_nd$prate <- z_nd$fit/z_nd$pyears
z_nd %>% ggplot(aes(x=day, y=prate, colour=vacc)) + geom_line()

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
  geom_hline(yintercept = 0.5, linetype = 3) +
  annotate("text", x=80, y=0.6, label = "Waning threshold", size=3) +
  scale_x_continuous(breaks = seq(14,84, by = 7), 
                     limits = c(0,70))
  

# Weekly z_rr
z_rr_wkly <- z_rr %>%
  #mutate(rr_upr = ifelse(rr_upr==3, Inf, rr_upr)) %>%
  filter(day %in% seq(14,84, by = 14)) %>%
  mutate(rr_ci = paste0(round(rr,2), " (", round(rr_lwr,2), ", ", round(rr_upr,2), ")")) %>%
  mutate(y=1)


p_c <- p_c +
  geom_text(aes(x=day, y=0, label = rr_ci), check_overlap = T, size=3, data=z_rr_wkly, inherit.aes = F)

p_c
if(z_vacc_type == "PB"){
  p_c1 = p_c
}


png(file=paste0("./output/final/modelling/", z_event_endpoint, "/plots/gam_RR_",z_vacc_type,".png"),
    width = 800, height=400)

#gridExtra::grid.arrange(p_c, gam_rrs, ncol=1, heights = c(3,0.5))
p_c
dev.off() 


}



model_outputs_vacc("PB")
model_outputs_vacc("AZ")



## Plot GAMS together
png(file=paste0("./output/final/modelling/", z_event_endpoint, "/plots/gam_RR_both.png"),
    width = 800, height=600)

cowplot::plot_grid(p_c, p_c1, labels = "AUTO", ncol=1)

dev.off() 



##### GAM - Age and Sex ######

df_cc_ps_matches <- df_cc_ps_matches %>%
  mutate(age_grp2 = ifelse(age_grp == "18-64", "18-64", "65+")) %>%
  left_join(select(df_cohort, EAVE_LINKNO, Sex))


## Functions
# Age
gam_age_fn <- function(z_vacc_type, z_age){

if(z_vacc_type == "PB"){
  z_vacc_title = "BNT162b2"
} else {
  z_vacc_title = "ChAdOx1"
}

### C: GAM
z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,seq(0,max(df_cc_ps_matches$time_to_hosp),by=1) ))


#aggregate for overall

z.agg <- pyears(Surv(time_to_hosp,event) ~ vacc + z.yr +vacc_type + age_grp2,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)



z_pois <- z.agg$data %>%
  mutate(day =  as.numeric(z.yr)-1) %>%
  filter(vacc_type == z_vacc_type) %>%
  filter(age_grp2 == z_age)


z_gam_vacc <- gam(event ~ offset(log(pyears)) + vacc + s(day, by=vacc), 
                  knots = list(day = seq(14,70, by =7)), 
                  family=poisson, data=z_pois, subset= vacc_type == z_vacc_type)



# Create cut off so length of vacc and uv match
z <- z_pois %>%
  group_by(z.yr) %>%
  summarise(n=n()) %>%
  filter(n == 2)



## Plot GAM RR
#for using for predictions
z_nd <- z_pois %>%
  filter(z.yr %in% z$z.yr)

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
               day = 1:length(difference),
               rr_ci = paste0(round(rr,2), " (", round(rr_lwr,2), ", ", round(rr_upr,2), ")"))

z_rr <- z_rr %>%
  left_join(filter(z_pois, vacc == "vacc") %>%
              select(event, day) %>%
              rename(event_vacc =event))
z_rr <- z_rr %>%
  left_join(filter(z_pois, vacc == "uv") %>%
              select(event, day) %>%
              rename(event_uv =event))


## P-value for quadratic term post 14 days
z_quad <- glm(event ~ offset(log(pyears)) +  day*vacc + I(day^2)*vacc, 
              family=poisson, data=z_pois, subset = vacc_type == z_vacc_type & day > 14)

zz_p <-predict.glm(z_quad, newdata=z_nd, type="response")
z_nd$fit <- zz_p
z_nd$prate <- z_nd$fit/z_nd$pyears

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
       subtitle = paste0("Age group: ", z_age),
       caption = paste0("Vaccinated quadratic coefficient after 14 days: ",round(z_quad_sum[6,1],4),
                        #" (", round(upr,3), ", ",round(lwr,3),"), 
                        ", p.value = ", round(z_quad_sum[6,4],3))) +
  geom_vline(xintercept = 14, linetype=2) +
  annotate("text", x=14.5, y=2, label = "14 days", hjust=0, size=3.5)+
  theme_light() +
  lims(y=c(0,3)) +
  geom_hline(yintercept = 1, linetype = 1) +
  geom_hline(yintercept = 0.5, linetype = 3) +
  annotate("text", x=80, y=0.6, label = "Waning threshold", size=3) +
  scale_x_continuous(breaks = seq(14,84, by = 7), 
                     limits = c(0,70))


# Weekly z_rr
z_rr_wkly <- z_rr %>%
  #mutate(rr_upr = ifelse(rr_upr==3, Inf, rr_upr)) %>%
  filter(day %in% seq(14,84, by = 14)) %>%
  mutate(rr_ci = paste0(round(rr,2), " (", round(rr_lwr,2), ", ", round(rr_upr,2), ")")) %>%
  mutate(y=1)


p_c <- p_c +
  geom_text(aes(x=day, y=0, label = rr_ci), check_overlap = T, size=2, data=z_rr_wkly, inherit.aes = F)

p_c

}

# Sex
gam_sex_fn <- function(z_vacc_type, z_sex){
  
  if(z_vacc_type == "PB"){
    z_vacc_title = "BNT162b2"
  } else {
    z_vacc_title = "ChAdOx1"
  }
  
  ### C: GAM
  z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,seq(0,max(df_cc_ps_matches$time_to_hosp),by=1) ))
  
  
  #aggregate for overall
  
  z.agg <- pyears(Surv(time_to_hosp,event) ~ vacc + z.yr +vacc_type + Sex,
                  data=df_cc_ps_matches , scale=1, data.frame=TRUE)
  
  
  
  z_pois <- z.agg$data %>%
    mutate(day =  as.numeric(z.yr)-1) %>%
    filter(vacc_type == z_vacc_type) %>%
    filter(Sex == z_sex)
  
  
  z_gam_vacc <- gam(event ~ offset(log(pyears)) + vacc + s(day, by=vacc), 
                    knots = list(day = seq(14,70, by =7)), 
                    family=poisson, data=z_pois, subset= vacc_type == z_vacc_type)
  
  # Create cut off so length of vacc and uv match
  z <- z_pois %>%
    group_by(z.yr) %>%
    summarise(n=n()) %>%
    filter(n == 2)
  
  
  
  ## Plot GAM RR
  #for using for predictions
  z_nd <- z_pois %>%
    filter(z.yr %in% z$z.yr)
  
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
                 day = 1:length(difference),
                 rr_ci = paste0(round(rr,2), " (", round(rr_lwr,2), ", ", round(rr_upr,2), ")"))
  
  z_rr <- z_rr %>%
    left_join(filter(z_pois, vacc == "vacc") %>%
                select(event, day) %>%
                rename(event_vacc =event))
  z_rr <- z_rr %>%
    left_join(filter(z_pois, vacc == "uv") %>%
                select(event, day) %>%
                rename(event_uv =event))
  
  
  ## P-value for quadratic term post 14 days
  z_quad <- glm(event ~ offset(log(pyears)) +  day*vacc + I(day^2)*vacc, 
                family=poisson, data=z_pois, subset = vacc_type == z_vacc_type & day > 14)
  
  zz_p <-predict.glm(z_quad, newdata=z_nd, type="response")
  z_nd$fit <- zz_p
  z_nd$prate <- z_nd$fit/z_nd$pyears
  
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
         subtitle = paste0("Sex: ", z_sex),
         caption = paste0("Vaccinated quadratic coefficient after 14 days: ",round(z_quad_sum[6,1],4),
                          #" (", round(upr,3), ", ",round(lwr,3),"), 
                          ", p.value = ", round(z_quad_sum[6,4],3))) +
    geom_vline(xintercept = 14, linetype=2) +
    annotate("text", x=14.5, y=2, label = "14 days", hjust=0, size=3.5)+
    theme_light() +
    lims(y=c(0,3)) +
    geom_hline(yintercept = 1, linetype = 1) +
    geom_hline(yintercept = 0.5, linetype = 3) +
    annotate("text", x=80, y=0.6, label = "Waning threshold", size=3) +
    scale_x_continuous(breaks = seq(14,84, by = 7), 
                       limits = c(0,70))
  
  
  # Weekly z_rr
  z_rr_wkly <- z_rr %>%
    #mutate(rr_upr = ifelse(rr_upr==3, Inf, rr_upr)) %>%
    filter(day %in% seq(14,84, by = 14)) %>%
    mutate(rr_ci = paste0(round(rr,2), " (", round(rr_lwr,2), ", ", round(rr_upr,2), ")")) %>%
    mutate(y=1)
  
  
  p_c <- p_c +
    geom_text(aes(x=day, y=0, label = rr_ci), check_overlap = T, size=2, data=z_rr_wkly, inherit.aes = F)
  
  p_c
  
}



## Combine together
library("cowplot")

# AZ
p_title <- ggdraw() +
  draw_label(paste0("Rate ratios of ", z_title," for ChAdOx1 by age and sex"),
             size=10)
p1 <- gam_age_fn("AZ", "18-64")
p2 <- gam_age_fn("AZ", "65+")
p3 <- gam_sex_fn("AZ", "F")
p4 <- gam_sex_fn("AZ", "M")


png(file=paste0("./output/final/modelling/", z_event_endpoint, "/plots/gam_RR_agesex_AZ.png"),
    width = 1000, height=600)

plot_grid(p_title, plot_grid(p1, p2, p3, p4, labels = "AUTO", ncol=2), ncol=1,
          rel_heights = c(0.5,5))

dev.off() 


## PB
p_title <- ggdraw() +
  draw_label(paste0("Rate ratios of ", z_title," for BNT162b2 by age and sex"),
             size=10)
p1 <- gam_age_fn("PB", "18-64")
p2 <- gam_age_fn("PB", "65+")
p3 <- gam_sex_fn("PB", "F")
p4 <- gam_sex_fn("PB", "M")


png(file=paste0("./output/final/modelling/", z_event_endpoint, "/plots/gam_RR_agesex_PB.png"),
    width = 1000, height=500)

plot_grid(p_title, plot_grid(p1, p2, p3, p4, labels = "AUTO", ncol=2), ncol=1,
          rel_heights = c(0.5,6))

dev.off() 


##### GAMS with different knots #####
z_vacc_type <- "PB"
z_vacc_type <- "AZ"

if(z_vacc_type == "PB"){
  z_vacc_title = "BNT162b2"
} else {
  z_vacc_title = "ChAdOx1"
}

z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1,seq(0,max(df_cc_ps_matches$time_to_hosp),by=1) ))


#aggregate for overall

z.agg <- pyears(Surv(time_to_hosp,event) ~ vacc + z.yr +vacc_type,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)

#z.agg$data %>% mutate(long = as.numeric(z.yr) >70) %>% group_by(vacc, vacc_type, long) %>% summarise(N=sum(event))


z_pois <- z.agg$data %>%
  mutate(day =  as.numeric(z.yr)-1) %>%
  filter(vacc_type == z_vacc_type) %>%
filter(day %in% 14:56)

z_pois %>%
  filter(vacc == "vacc") %>%
  group_by(vacc_type) %>%
  summarise(events = sum(event))

z_gam_vacc <- gam(event ~ offset(log(pyears)) + vacc + s(day, by=vacc, k = 3), 
                 # knots = list(day = seq(14,70, by =7)), 
                  #knots = list(day = seq(14,56, by =14)), 
                  #knots = list(day = 14),
                  family=poisson, data=z_pois, subset= vacc_type == z_vacc_type)

# Compare linearity
#z_gam_vacc2 <- gam(event ~ offset(log(pyears)) + vacc + s(day), 
# family=poisson, data=z_pois, subset= vacc_type == z_vacc_type)


#AIC(z_gam_vacc, z_gam_vacc2)


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
#write.csv(z_rr, paste0("./output/final/modelling/", z_event_endpoint, "/poisson_gam_", z_vacc_type, ".csv"))

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
  geom_hline(yintercept = 0.5, linetype = 3) +
  annotate("text", x=80, y=0.6, label = "Waning threshold", size=3) +
  scale_x_continuous(breaks = seq(14,84, by = 7), 
                     limits = c(0,70))


# Weekly z_rr
z_rr_wkly <- z_rr %>%
  #mutate(rr_upr = ifelse(rr_upr==3, Inf, rr_upr)) %>%
  filter(day %in% seq(14,84, by = 14)) %>%
  mutate(rr_ci = paste0(round(rr,2), " (", round(rr_lwr,2), ", ", round(rr_upr,2), ")")) %>%
  mutate(y=1)


p_c <- p_c +
  geom_text(aes(x=day, y=0, label = rr_ci), check_overlap = T, size=3, data=z_rr_wkly, inherit.aes = F)

p_c







##### Previous GAM version (by hand) - do not use #####
# time_to_hosp =z.year as days rather than weeks
z.yr <- tcut(rep(0,nrow(df_cc_ps_matches)), c(-1:max(df_cc_ps_matches$time_to_hosp) ) )

#aggregate for overall
z.agg <- pyears(Surv(time_to_hosp,event) ~ vacc + z.yr +vacc_type,
                data=df_cc_ps_matches , scale=1, data.frame=TRUE)


z_pois <- z.agg$data %>%
  mutate(day =  as.numeric(z.yr)-1) %>%
  filter(!(day %in% c(1:13)) &
           vacc_type == z_vacc_type)

z_gam_vacc <- gam(event ~ offset(log(pyears)) + -1 + vacc + s(as.numeric(z.yr), by=vacc), 
                  family=poisson, data=z_pois, subset= vacc_type == z_vacc_type)

# Ratio of the two lines for un and vacc
# Predict risk of event
z_pois_pred <- z_pois %>%
  filter(vacc_type == z_vacc_type)
z_pred <- predict(z_gam_vacc, newdata=z_pois_pred, type="response", se =T) # lpmatrix
z_pois_pred$est <- z_pred$fit
z_pois_pred$mu <- (z_pred$fit)^2
z_pois_pred$var.fit <- (z_pred$se.fit)^2
z_pois_pred$upr <- z_pred$fit + 1.96*z_pred$se.fit
z_pois_pred$lwr <- z_pred$fit - 1.96*z_pred$se.fit

# Ratios

z_pois_pred_rr <- z_pois_pred %>%
  group_by(z.yr) %>%
  mutate(RR = est[vacc == "vacc"]/est[vacc == "uv"]) %>%
  mutate(var1 = (mu[vacc == "vacc"])/(mu[vacc == "uv"])) %>%
  mutate(var2 = (var.fit[vacc == "vacc"])/(mu[vacc == "vacc"]) + 
           (var.fit[vacc == "uv"])/(mu[vacc == "uv"])) %>%
  mutate(var = var1*var2) %>%
  select(-c(vacc, pyears, n, event, est, upr, lwr, var.fit, var1, var2)) %>%
  distinct() %>%
  mutate(upr = RR + 1.96*sqrt(var)) %>%
  mutate(lwr=RR - 1.96*sqrt(var))

z_pois_pred_rr



#truncate for plot
z_pois_pred_rr$upr[z_pois_pred_rr$upr>=3] <- 3
z_pois_pred_rr$lwr[z_pois_pred_rr$lwr<0] <- 0

ggplot(z_pois_pred_rr) +
  geom_line(aes(x=as.numeric(z.yr), y=RR), col=eave_blue) +
  geom_ribbon(aes(x=as.numeric(z.yr), ymin=lwr, ymax=upr), alpha = 0.25, fill=eave_blue) +
  #scale_color_manual(values = c(eave_blue, eave_green), labels=c("Control", "Exposed"))+
  #scale_fill_manual(values = c(eave_blue, eave_green), labels=c("Control", "Exposed")) +
  labs(x="Days to event", fill="Vaccination group", 
       subtitle = paste0("GAM Relative Risks for ", z_vacc_type)) +
  geom_vline(xintercept = 14, linetype=2) +
  #annotate("text", x=14.5, y=150, label = "14 days", hjust=0, size=3.5)+
  theme_light() +
  lims(y=c(0,3), x=c(0,84)) +
  geom_hline(yintercept = 1, linetype = 1) +
  geom_hline(yintercept = 0.5, linetype = 3) +
  annotate("text", x=80, y=0.6, label = "Waning threshold", size=3)


## Old
# No longer needed
gam_rrs <- ggplot(z_rr_wkly) +
  #geom_point(aes(x=days, y= 0.9), shape =15) +
  geom_text(aes(x=day, y=1, label = rr_ci), check_overlap = T, size=3.5) +
  #geom_text(aes(x=day, y=0.97, label = event), check_overlap = T, size=3.5) +
  labs(y="", x="") +
  lims(y=c(1.01, 0.95), x=c(14,84)) +
  #annotate("text", x=1, y=0.97, label="No. events:")+
  annotate("text", x=1, y=1, label="RR (95% CI):")+
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(angle = 0)) 