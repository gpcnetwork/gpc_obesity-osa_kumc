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

# webshot::install_phantomjs()

source_url("https://raw.githubusercontent.com/sxinger/utils/master/analysis_util.R")
path_to_data_folder<-file.path(
  getwd(),
  "data"
)

# overall summary
var_lst<-c(
   "AGEGRP"
  ,"SEXF"
  ,"RACE_LABEL"
  ,"LIS_DUAL_IND"
  ,'CCI_SCORE'
  ,'OBES'
  ,'T2DM'
  ,'CKD'
  ,'COPD'
  ,'AFIB'
  ,'HTN'
  ,'ND'
  ,'HYSM'
  ,'INSM'
  ,'ACG'
  ,'AHT'
  ,'ALP'
  ,'BGR'
  ,paste0(c('MI','HF','STROKE','REVASC',"MACE"),'_HIST')
  ,paste0(c('MI','HF','STROKE','REVASC',"MACE",'DEATH'),'_status')
  ,paste0(c('MI','HF','STROKE','REVASC',"MACE",'DEATH'),'_time')
  ,'DAYS_OSA_TO_CENSOR'
)
numvar_lst<-c(
   'CCI_SCORE'
  ,paste0(c('MI','HF','STROKE','REVASC',"MACE","DEATH"),'_time')
  ,'DAYS_OSA_TO_CENSOR'
)
facvar_lst<-var_lst[!var_lst %in% numvar_lst]

#===== exposure analysis ========================
#=== surv =======================================
df<-readRDS(file.path(
  path_to_data_folder,
  "cpap_exposure_final.rda"
))

# overview
desc_cohort<-univar_analysis_mixed(
  df = df,
  id_col="PATID",
  grp = 1,
  var_lst = var_lst,
  facvar_lst = facvar_lst,
  pretty = T
)

desc_cohort %>%
  save_kable(
    file.path(getwd(),"res/expos_all.pdf")
  )

# comparison
desc_case_ctrl<-univar_analysis_mixed(
  df = df,
  id_col = "PATID",
  grp = df$adherence_yr1_qt,
  var_lst = var_lst,
  facvar_lst = facvar_lst,
  pretty = T
)
desc_case_ctrl %>%
  save_kable(
    file.path(getwd(),"res/expos_qt4_pap.pdf")
  )

#=== mace ========================================
var_lst2<-var_lst[!grepl("(_HIST)+",var_lst)]
facvar_lst2<-facvar_lst[!grepl("(_HIST)+",facvar_lst)]

df2<-df %>% 
  filter(
    MACE_HIST==0 & 
    (DAYS_OSA_TO_MACE>DAYS_OSA_TO_CPAP_INIT |
     is.na(DAYS_OSA_TO_MACE) |
     is.na(DAYS_OSA_TO_CPAP_INIT)
     )
    )

desc_case_ctrl2<-univar_analysis_mixed(
  df = df2,
  id_col = "PATID",
  grp = df2$adherence_yr1_qt,
  var_lst = var_lst2,
  facvar_lst = facvar_lst2,
  pretty = T
)
desc_case_ctrl2 %>%
  save_kable(
    file.path(getwd(),"res/expos_mace_pap.pdf")
  )

#===== adherence analysis ========================
#=== surv =======================================
df<-readRDS(file.path(
  path_to_data_folder,
  "cpap_adherence_final.rda"
))

# overview
desc_cohort<-univar_analysis_mixed(
  df = df,
  id_col = "PATID",
  grp = 1,
  var_lst = var_lst,
  facvar_lst = facvar_lst,
  pretty = T
)
desc_cohort %>%
  save_kable(
    file.path(getwd(),"res/adhrn_all.pdf")
  )

# comparison
desc_case_ctrl<-univar_analysis_mixed(
  df = df,
  id_col = "PATID",
  grp = df3$adherence_yr1_mttree,
  var_lst = var_lst,
  facvar_lst = facvar_lst,
  pretty=F
)
desc_case_ctrl %>% 
  save_kable(
    file.path(getwd(),"res/adhrn_pap.pdf")
  )

#=== mace =======================================
var_lst2<-var_lst[!grepl("(_HIST)+",var_lst)]
facvar_lst2<-facvar_lst[!grepl("(_HIST)+",facvar_lst)]

df2<-df %>% 
  filter(
    MACE_HIST==0 & 
      (DAYS_OSA_TO_MACE>DAYS_OSA_TO_CPAP_INIT |
         is.na(DAYS_OSA_TO_MACE) |
         is.na(DAYS_OSA_TO_CPAP_INIT)
      )
  )

desc_case_ctrl2<-univar_analysis_mixed(
  df = df2,
  id_col = "PATID",
  grp = df2$adherence_yr1_qt,
  var_lst = var_lst2,
  facvar_lst = facvar_lst2,
  pretty = T
)
desc_case_ctrl2 %>%
  save_kable(
    file.path(getwd(),"res/expos_qt4_mace_pap.pdf")
  )

