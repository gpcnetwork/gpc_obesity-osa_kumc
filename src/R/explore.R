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
  tidycmprsk,
  webshot2
)

source_url("https://raw.githubusercontent.com/sxinger/utils/master/analysis_util.R")

path_to_datadir<-file.path(getwd(),"data")
path_to_outdir<-file.path(getwd(),"res")

# overall summary
var_lst<-c(
   "AGE_AT_OSA_DX1"
  ,"AGEGRP"
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
   'AGE_AT_OSA_DX1'
  ,'CCI_SCORE'
  ,paste0(c('MI','HF','STROKE','REVASC',"MACE","DEATH"),'_time')
  ,'DAYS_OSA_TO_CENSOR'
)
facvar_lst<-var_lst[!var_lst %in% numvar_lst]

#===== exposure analysis ========================
#=== surv =======================================
df<-readRDS(file.path(
  path_to_datadir,
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
    file.path(path_to_outdir,"expos_all.html")
  )


# comparison
desc_case_ctrl<-univar_analysis_mixed(
  df = df,
  id_col = "PATID",
  grp = df$CPAP_IND,
  var_lst = var_lst,
  facvar_lst = facvar_lst,
  pretty = T
)
desc_case_ctrl %>%
  save_kable(
    file.path(path_to_outdir,"expos_pap.html")
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
  grp = df2$CPAP_IND,
  var_lst = var_lst2,
  facvar_lst = facvar_lst2,
  pretty = T
)
desc_case_ctrl2 %>%
  save_kable(
    file.path(path_to_outdir,"expos_mace_pap.html")
  )

#===== adherence analysis ========================
#=== surv =======================================
df<-readRDS(file.path(
  path_to_datadir,
  "cpap_adherence_final2.rda"
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
    file.path(path_to_outdir,"adhrn_all.html")
  )

# comparison
desc_case_ctrl<-univar_analysis_mixed(
  df = df,
  id_col = "PATID",
  grp = df$adherence_yr1_qt,
  var_lst = var_lst,
  facvar_lst = facvar_lst,
  pretty=T
)
desc_case_ctrl %>% 
  save_kable(
    file.path(path_to_outdir,"adhrn_qt_pap.html")
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
    file.path(path_to_outdir,"adhrn_qt_mace_pap.html")
  )

