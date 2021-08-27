# EAVE II Analysis Repository
## COVID-19 vaccine waning 

### Scripts

**00_functions:** Functions used in the analysis 

**01_data_setup:** Initialises the baseline cohort (df_cohort), the vaccination data (df_vaccination) and the event datasets 

**02_descriptive:** Descriptive analysis on the characteristics of the baseline cohort 

**03_matching…:** Creates matched cohort, summarises their characteristics and models the cohort 

- ***03a_matching_time_varying:*** Performs time-varying propensity score matching between 1st dose and unvaccinated

- ***03b_matching_cc:*** Reformats matching cohort to stack matches and outputs the main dataframe for all analyses (df_cc_...)

- ***03c_matching_ps_summary:*** Descriptive analyses on matched cohort including descriptive tables and covariate balance plots

- ***03d_matching_modelling:*** Models data using cumulative risk plots, GLMs and GAMs

**04_cohort_analysis…:** Creates and models the cohort analysis 

**05_dataflow:** Calculates cohort sizes for data flow diagram

**06_checking_outputs:** Sense checks the matching to ensure it has worked

**07_meta_analysis:** Creates relevant tables for the pooled meta analysis

**08_2nd_dose...:** Creates matched cohort for the 2nd dose waning (same structure as 03_matching...)
