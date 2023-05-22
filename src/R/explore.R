#################################################################
# Copyright (c) 2021-2025 University of Missouri                   
# Author: Xing Song, xsm7f@umsystem.edu                            
# File: explore.R
# Description: exploratory analysis 
#################################################################
rm(list=ls())
pacman::p_load(
  tidyverse,
  magrittr,
  broom,
  survival,
  survminer,
  kableExtra,
  devtools,
  ggrepel,
  tidycmprsk
)

source_url("https://raw.githubusercontent.com/sxinger/utils/master/analysis_util.R")
path_to_data_folder<-file.path(
  getwd(),
  "osa-obesity-cpap-adherence",
  "data"
)

# overall summary
var_lst<-c(
  "AGE_AT_OSA_DX1","AGEGRP","SEXF","RACE_LABEL","LIS_DUAL_IND",
  "T2DM_HISTORY","HTN_HISTORY","OBESITY_HISTORY",
  "MI_HISTORY","STROKE_HISTORY","HF_HISTORY","REVASC_HISTORY","MACE_HISTORY",
  "bariatric_surgery","chronic_obstructive_pulmonary_disease","hypersomnia","insomina",
  "DEATH_time"
)
facvar_lst<-var_lst[!grepl("(age_at|time)",tolower(var_lst))]

var_lst2<-c(
  "AGE_AT_OSA_DX1","AGEGRP","SEXF","RACE_LABEL","LIS_DUAL_IND",
  "T2DM_HISTORY","HTN_HISTORY","OBESITY_HISTORY",
  "bariatric_surgery","chronic_obstructive_pulmonary_disease","hypersomnia","insomina",
  "MACE_time"
)
facvar_lst2<-var_lst2[!grepl("(age_at|time)",tolower(var_lst2))]

#===== exposure analysis ========================
#=== surv =======================================
df<-readRDS(file.path(path_to_data_folder,"cpap_exposure_aset.rda")) %>%
  filter(DEATH_time>0)

# overview
desc_cohort<-univar_analysis_mixed(
  df,
  id_col="PATID",
  grp=1,
  var_lst=var_lst,
  facvar_lst=facvar_lst,
  pretty=F
)
desc_cohort %>% View

# comparison
desc_case_ctrl<-univar_analysis_mixed(
  df,
  id_col="PATID",
  grp=df$CPAP_IND,
  var_lst=var_lst,
  facvar_lst=facvar_lst,
  pretty=F
)
desc_case_ctrl %>% View

#=== mace ========================================
df2<-df %>% filter(MACE_HISTORY==0&MACE_time>0)

desc_case_ctrl<-univar_analysis_mixed(
  df2,
  id_col="PATID",
  grp=df2$CPAP_IND,
  var_lst=var_lst2,
  facvar_lst=facvar_lst2,
  pretty=F
)
desc_case_ctrl %>% View


#===== adherence analysis ========================
#=== surv =======================================
df3<-readRDS(file.path(path_to_data_folder,"cpap_dose_response_aset.rda")) %>%
  filter(DEATH_time>0)

# overview
desc_cohort<-univar_analysis_mixed(
  df3,
  id_col="PATID",
  grp=1,
  var_lst=var_lst2,
  facvar_lst=facvar_lst2,
  pretty=F
)
desc_cohort %>% View

# comparison
desc_case_ctrl<-univar_analysis_mixed(
  df3,
  id_col="PATID",
  grp=df3$adherence_yr1_mttree,
  var_lst=var_lst2,
  facvar_lst=facvar_lst2,
  pretty=F
)
desc_case_ctrl %>% View

# unadjusted KM
risk_tbl<-summary(survfit(Surv(DEATH_time,DEATH_status) ~ adherence_yr1_mttree, data = df3),
                  times = 365*c(1:5))
km_mort_unadj<-ggsurvplot(
  fit = survfit(Surv(DEATH_time,DEATH_status) ~ adherence_yr1_mttree, data = df3),
  pval = TRUE, 
  conf.int = TRUE,
  risk.table = TRUE,
  linetype = "strata",
  break.x.by = 365,
  xlab = "Days", 
  ylab = "Mortality Endpoint")

km_mort_unadj$plot +
  geom_vline(xintercept=365*c(1:5),linetype=2)+
  geom_label_repel(data=data.frame(x=risk_tbl$time,
                                   y=risk_tbl$surv,
                                   label=round(risk_tbl$surv,2),
                                   label_int=paste0(round(risk_tbl$surv,2),"[",round(risk_tbl$lower,2),",",round(risk_tbl$upper,2),"]")),
                   aes(x=x,y=y,label=label),
                   max.overlaps = 40)

#=== mace =======================================
df4<-df3 %>% filter(MACE_HISTORY==0&MACE_time>0)

# overview
desc_cohort<-univar_analysis_mixed(
  df4,
  id_col="PATID",
  grp=1,
  var_lst=var_lst,
  facvar_lst=facvar_lst,
  pretty=F
)
desc_cohort %>% View

# comparison
desc_case_ctrl<-univar_analysis_mixed(
  df4,
  id_col="PATID",
  grp=df4$adherence_yr1_mctree,
  var_lst=var_lst,
  facvar_lst=facvar_lst,
  pretty=F
)
desc_case_ctrl %>% View


#unadjusted KM
risk_tbl<-summary(survfit(Surv(MACE_time,MACE_status) ~ adherence_yr1_tt, data = df4),
                  times = 365*c(1:5))
km_mort_unadj<-ggsurvplot(
  fit = survfit(Surv(MACE_time,MACE_status) ~ adherence_yr1_tt, data = df4),
  pval = TRUE, 
  conf.int = TRUE,
  risk.table = TRUE,
  linetype = "strata",
  break.x.by = 365,
  xlab = "Days", 
  ylab = "MACE Endpoint")

km_mort_unadj$plot +
  geom_vline(xintercept=365*c(1:5),linetype=2)+
  geom_label_repel(data=data.frame(x=risk_tbl$time,
                                   y=risk_tbl$surv,
                                   label=round(risk_tbl$surv,2),
                                   label_int=paste0(round(risk_tbl$surv,2),"[",round(risk_tbl$lower,2),",",round(risk_tbl$upper,2),"]")),
                   aes(x=x,y=y,label=label),
                   max.overlaps = 40)

