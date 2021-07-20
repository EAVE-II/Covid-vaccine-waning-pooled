###### Cohort analysis modelling #####


library(tidyverse)
#set-up directories
Location <- "/conf/"  # Server
project_path <- paste0(Location,"EAVE/GPanalysis/progs/RM/Vaccine/Vaccine_waning")

#df <- readRDS(paste0(project_path, "/output/df_long_death_hosp.RDS"))
#df <- readRDS(paste0(project_path, "/output/df_long_death_hosp_12wks.RDS"))
df <- readRDS(paste0("./output/df_long_death_hosp_10wks.RDS"))


## Eliminate care homes 
Cohort_Household <- readRDS("/conf/EAVE/GPanalysis/outputs/temp/Cohort_Household.rds") %>%
  mutate(n_hh_gp = cut(n_hh, breaks=c(0,1,2,5,10,30,100,max(n_hh)),
                       labels=c("1", "2", "3-5", "6-10", "11-30", "31-100", "101+")))%>% 
  mutate(ave_hh_age=if_else(is.na(ave_hh_age), mean(ave_hh_age, na.rm=T), ave_hh_age) )

# Add in cohort household information
df <- df %>%
  left_join(select(Cohort_Household, EAVE_LINKNO,
                   n_hh_gp, ave_hh_age, care_home_elderly)) %>%
  filter(care_home_elderly == 0)



print(z_event_endpoint)
print(table(df$event))
print(table(df$vacc_status, exclude=NULL))
print(table(df$event, df$vacc_type))
print(table(df$event, df$vacc_status))

z_combine_vacc_level <- FALSE  # to combine the first 1 or 2 levels of vacc_status with uv
output_list$combine_levels <- z_combine_vacc_level
z_adjustment <- "full"  #or "minimal" "full"
output_list$adjustment <- z_adjustment
output_list$prop_score <- "inverse propensity weighting"
#"no propensity weighting"  "inverse propensity weighting"
output_list$model <- "glm" 
#"cox"  "glm"

#data_selection_flags <- read.csv("02_data_selection_flags.csv")

library(readr)
data_selection_flags <- read.csv("./02_data_selection_flags.csv")

data_selection_flags <- data_selection_flags %>% filter(!grepl("omit", notes))

##### Run '03b_Functions.R' ######


source("./Eles_code/03b_Functions.R")


#rm(z_res)




#### Find estimates for data selection flags ####
i <- 1
i <- 6
i <- 10

z_out_list <- list()

for (i in c(1,6,10)) {
  
  z_df <- df
  

  ## Go through the conditions in the data selection flag row
  print(data_selection_flags[i,])
  
  ## Filter long data to vaccine type if selection requires it
  #vacc_type
  if (data_selection_flags[i,"vacc_type"] == "PB") z_df <- filter(z_df, vacc_type %in% c("uv","PB") )
  if (data_selection_flags[i,"vacc_type"] == "AZ") z_df <- filter(z_df, vacc_type %in% c("uv","AZ") )
  
  ## Filter long data to age group  if selection requires it
  #age
  #z_df <- filter(z_df, ageYear >= data_selection_flags[i,"age_lower"] & ageYear <= data_selection_flags[i,"age_upper"]) 
  z_age_centre <- data_selection_flags[i,"age_lower"]
  if (z_age_centre < 30) z_age_centre <- 30
  
  ## Filter long data to Sex if selection requires it
  #sex
  if (data_selection_flags[i,"sex"] %in% c("F","M")) z_df <- filter(z_df, Sex==data_selection_flags[i,"sex"])
  
  ## Filter long data to Previous positive tests if selection requires it
  #previous positive tests
  #if (data_selection_flags[i,"omit_prev_tests"] == "all") z_df <- filter(z_df, test_before_dec8 %in% c("no pos test","post-vacc"))
  #if (data_selection_flags[i,"omit_prev_tests"] == "4w") z_df <- filter(z_df, !(test_before_dec8 %in% c("1+m")) )
  #if (data_selection_flags[i,"omit_prev_tests"] == "2w") z_df <- filter(z_df, !(test_before_dec8 %in% c("1+m", "3w","4w")))
  #if (data_selection_flags[i,"omit_prev_tests"] == "1w") z_df <- filter(z_df, test_before_dec8 %in% c("0-6d","no pos test","post-vacc"))
  
  z_df$n_tests_gp <- cut(z_df$n_tests, breaks = c(-1,0,1,2,3,9,100), labels=c("0","1","2","3","4-9","10+"))
  #z_df$age_gp <- relevel(z_df$age_gp, ref="60-64")
  #z_df$vacc_status <- relevel(z_df$vacc_status, ref="v1_7:13")
  #drop to 10 controsl per case
  #z_ids_case <- unique(filter(z_df, event==1)$EAVE_LINKNO)
  #z_ids_cont <- unique(z_df$EAVE_LINKNO)
  #z_ids_cont <- z_ids_cont[!(z_ids_cont %in% z_ids_case)]
  #z_ids <- c(z_ids_case, sample(z_ids_cont, size=round(length(z_ids_cont)*0.20)))
  #z_df <- filter(z_df, EAVE_LINKNO %in% z_ids)
  
  #z_df <- filter(z_df, Q_DIAG_CKD_LEVEL==0)
  #z.yr <- tcut(rep(0,nrow(analysis.df)), c(-1,seq(7,as.numeric(a.analysis.to-a.analysis.from),by=7)))
 
  
  ## Get person years
   #aggregate for overall
  z.agg <- pyears(Surv(start,stop,event) ~ vacc_status, data=z_df , weight=weight, 
                  scale=365.25, data.frame=TRUE)
  
  df_res <- z.agg$data
  names(df_res) <- c("Vacc.Status","Person.Years","Count","Events")
  df_res <- df_res %>% mutate(Rate = Events/Person.Years*365.25,
                              Person.Years = round(Person.Years,0) ) %>% 
    select(-Count)
   # mutate(RR = Rate/first(Rate))
  df_res
  sum(df_res$Events) - df_res$Events[1]
  
  z_var <- "vacc_status"
  
  #plot(survfit(Surv(start, stop, event) ~  strata(period),data=z_df), fun="event")
  
  ## Age adjustment

    z.agg <- pyears(Surv(start,stop,event) ~ vacc_status + age_gp, data=z_df, weight=weight , 
                    scale=365.25, data.frame=TRUE)
    z_pois <- z.agg$data
    z_m <- glm(event ~ offset(log(pyears)) + vacc_status + age_gp, family=poisson, data=z_pois)
    #summary(z_m)
    z.estimates.1 <- fun_ve_glm(z_m, z_type="") %>%
      rename_at(vars(-vacc_status), function(x) paste0(x, "_AGE"))

  
    
    
  ## Multiple covariate adjustment
      #z.agg <- pyears(Surv(start,stop,event) ~ vacc_status + period + age_gp + simd2020_sc_quintile +
        #                Sex + n_risk_gps +n_tests_gp + Council, data=z_df, 
         #             weight=weight , scale=365.25, data.frame=TRUE)
    variables_hosp <- c("age_gp" , "EAVE_BP","Q_DIAG_DIABETES_2"   , "Q_DIAG_COPD"  ,       
                        "Q_DIAG_CKD_LEVEL","Q_DIAG_DEMENTIA","Q_DIAG_STROKE","Q_LEARN_CAT"   ,      
                        "Q_DIAG_FRACTURE","Q_DIAG_NEURO","Q_DIAG_CCF","Q_DIAG_ASTHMA"    ,   
                        "Q_DIAG_EPILEPSY","Q_DIAG_BLOOD_CANCER","Q_DIAG_VTE","Q_DIAG_CIRRHOSIS" ,   
                        "Q_DIAG_RA_SLE","Q_DIAG_PVD","Q_DIAG_AF","Q_DIAG_PULM_RARE" ,   
                        "Q_DIAG_PARKINSONS","Q_DIAG_DIABETES_1","Q_DIAG_PULM_HYPER",
                        "simd2020_sc_quintile","Sex", "n_tests_gp", "Council")
    
    z_pois <- z_df %>% mutate(pyears=stop-start)
    z.fmla <- as.formula(paste("event"," ~ offset(log(pyears)) + vacc_status + period +",
                               paste(variables_hosp, collapse= "+")))
    z_m <- glm(formula = z.fmla  , data=z_pois, family=poisson, weight = eave_weight)

      
      z.estimates.2 <- fun_ve_glm(z_m, z_type="") %>%
        rename_at(vars(-vacc_status), function(x) paste0(x, "_FULLADJ"))
      
      
      
      ## Adjusted for age and inhospital status
      z.agg <- pyears(Surv(start,stop,event) ~ vacc_status + period + age_gp + in_hosp, data=z_df, weight=weight , 
                      scale=365.25, data.frame=TRUE)
      z_pois <- z.agg$data
      z_m <- glm(event ~ offset(log(pyears)) + vacc_status + age_gp + in_hosp, family=poisson, data=z_pois)
      #summary(z_m)
      z.estimates.3 <- fun_ve_glm(z_m, z_type="") %>%
        rename_at(vars(-vacc_status), function(x) paste0(x, "_ADJ2"))
      
      
      #### Adjusted full and minimal?
      #variables_hosp2 <- c(variables_hosp, "simd2020_sc_quintile",
                      #     "Sex", "n_tests_gp", "Council")
      #z_pois <- z_df %>% mutate(pyears=stop-start)
      #z.fmla <- as.formula(paste("event"," ~ offset(log(pyears)) + vacc_status + period +",
             #                    paste(variables_hosp2, collapse= "+")))
      #z_m <- glm(formula = z.fmla  , data=z_pois, family=poisson, weight = eave_weight)
      
      
      #z.estimates.4 <- fun_ve_glm(z_m, z_type="") %>%
       # rename_at(vars(-vacc_status), function(x) paste0(x, "_FULLADJ2"))
      


   #end full adjustment 
  
  #merge the summary data and estimates
  z_names_or <- c("vacc_status","OR","OR_LCL","OR_UCL")
  z_out <- df_res %>% 
    left_join(z.estimates.1[,!(startsWith(colnames(z.estimates.1), "VE"))],
              by=c("Vacc.Status" = "vacc_status")) %>% 
    left_join(z.estimates.2[,!(startsWith(colnames(z.estimates.2), "VE"))], by=c("Vacc.Status" = "vacc_status")) %>%
    left_join(z.estimates.3[,!(startsWith(colnames(z.estimates.3), "VE"))], by=c("Vacc.Status" = "vacc_status"))# %>%
    #left_join(z.estimates.4[,!(startsWith(colnames(z.estimates.4), "VE"))], by=c("Vacc.Status" = "vacc_status"))
  
  
  
  z_out
  z <- data_selection_flags[i,]
  for (i_row in 2:nrow(z_out)) z <- bind_rows(z, data_selection_flags[i,])
  
  z_out <- bind_cols(z,z_out)
  z_out
  #z_out$selections <- paste(output_list[[1]] , output_list[[5]], output_list[[6]], output_list[[8]])
  
  #if (exists("z_res")) z_res <- bind_rows(z_res,z_out) else z_res <- z_out
  
  z_out_list[[i]] <- z_out
  
} # end loop through the models


z_out_results <- z_out_list[c(1,6,10)] %>%
  reduce(full_join)

#write.csv(z_out_results, paste0("./output/final/cohort_analysis/cohort_analysis_tbl.csv"))
#write.csv(z_out_results, paste0("./output/final/cohort_analysis/cohort_analysis_tbl_12wks.csv"))
write.csv(z_out_results, paste0("./output/final/cohort_analysis/cohort_analysis_tbl_10wks.csv"))
write.csv(z_out_results, paste0("./output/final/cohort_analysis/cohort_analysis_tbl_8wks_full.csv"))

