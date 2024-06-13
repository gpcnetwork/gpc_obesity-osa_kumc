#################################################################
# Copyright (c) 2021-2025 University of Missouri                   
# Author: Xing Song, xsm7f@umsystem.edu                            
# File: preproc.R
# Description: data preprocess to form final analytic set with
#              one patient per row
#################################################################

# install.packages("pacman")
pacman::p_load(
  tidyverse,
  magrittr,
  rpart,
  rpart.plot
)

rm(list=ls()); gc()
setwd("C:/repo/gpc-obesity-osa")

path_to_data_folder<-file.path(
  getwd(),
  "data"
)

##==== preprocess cohort for exposure analysis ========
# observation window: before osa onset
# routine data preprocess
dat_proc<-readRDS(file.path(path_to_data_folder,"pap_exposure_aset.rda")) %>%
  # survival object
  mutate(
    # MACE survival obj
    MACE_time = coalesce(DAYS_OSA_TO_MACE,pmin(DAYS_OSA_TO_DEATH,DAYS_OSA_TO_CENSOR,na.rm=T)),
    MACE_status = as.numeric(!is.na(DAYS_OSA_TO_MACE)),
    MACE_status_sub = case_when(
      !is.na(DAYS_OSA_TO_DEATH)&is.na(DAYS_OSA_TO_MACE) ~ 2,
      TRUE ~ as.numeric(!is.na(DAYS_OSA_TO_MACE))
    ),
    # MI survival obj
    MI_time = coalesce(DAYS_OSA_TO_MI,pmin(DAYS_OSA_TO_DEATH,DAYS_OSA_TO_CENSOR,na.rm=T)),
    MI_status = as.numeric(!is.na(DAYS_OSA_TO_MI)),
    MI_status_sub = case_when(
      !is.na(DAYS_OSA_TO_DEATH)&is.na(DAYS_OSA_TO_MI) ~ 2,
      TRUE ~ as.numeric(!is.na(DAYS_OSA_TO_MI))
    ),
    # STROKE survival obj
    STROKE_time = coalesce(DAYS_OSA_TO_STROKE,pmin(DAYS_OSA_TO_DEATH,DAYS_OSA_TO_CENSOR,na.rm=T)),
    STROKE_status = as.numeric(!is.na(DAYS_OSA_TO_STROKE)),
    STROKE_status_sub = case_when(
      !is.na(DAYS_OSA_TO_DEATH)&is.na(DAYS_OSA_TO_STROKE) ~ 2,
      TRUE ~ as.numeric(!is.na(DAYS_OSA_TO_STROKE))
    ),
    # HF survival obj
    HF_time = coalesce(DAYS_OSA_TO_HF,pmin(DAYS_OSA_TO_DEATH,DAYS_OSA_TO_CENSOR,na.rm=T)),
    HF_status = as.numeric(!is.na(DAYS_OSA_TO_HF)),
    HF_status_sub = case_when(
      !is.na(DAYS_OSA_TO_DEATH)&is.na(DAYS_OSA_TO_HF) ~ 2,
      TRUE ~ as.numeric(!is.na(DAYS_OSA_TO_HF))
    ),
    # REVASC survival obj
    REVASC_time = coalesce(DAYS_OSA_TO_REVASC,pmin(DAYS_OSA_TO_DEATH,DAYS_OSA_TO_CENSOR,na.rm=T)),
    REVASC_status = as.numeric(!is.na(DAYS_OSA_TO_REVASC)),
    REVASC_status_sub = case_when(
      !is.na(DAYS_OSA_TO_DEATH)&is.na(DAYS_OSA_TO_REVASC) ~ 2,
      TRUE ~ as.numeric(!is.na(DAYS_OSA_TO_REVASC))
    ),
    # all-cause DEATH survival obj
    DEATH_time = coalesce(DAYS_OSA_TO_DEATH,DAYS_OSA_TO_CENSOR),
    DEATH_status = as.numeric(!is.na(DAYS_OSA_TO_DEATH))
  ) %>%
  # pre-MACE event indicator
  mutate(
    MACE_HIST = as.numeric(!is.na(MACE_HIST)),
    MI_HIST = as.numeric(!is.na(MI_HIST)),
    STROKE_HIST = as.numeric(!is.na(STROKE_HIST)),
    HF_HIST = as.numeric(!is.na(HF_HIST)),
    REVASC_HIST = as.numeric(!is.na(REVASC_HIST))
  ) %>%
  # cleanup demographic fields
  mutate(
    RACE_LABEL = case_when(
      HISPANIC=='Y' ~ 'Hispanic',
      RACE=='05' ~ 'White',
      RACE=='03' ~ 'AA',
      RACE=='02' ~ 'Asian',
      RACE=='01' ~ 'AI',
      RACE=='NI' ~ 'Unknown',
      TRUE ~ 'Other'
    ),
    RACE_LABEL = relevel(as.factor(RACE_LABEL),ref="White"),
    RACE_LABEL_fac = paste0("RACE_",relevel(as.factor(RACE_LABEL),ref="White")),race_ind = 1,
    SEXF = as.numeric(SEX=='F'),
    AGEGRP = case_when(AGE_AT_OSA_DX1>=80 ~ 'agegrp4',
                       TRUE ~ paste0('agegrp',floor((AGE_AT_OSA_DX1-65)/5)+1)),
    AGEGRP_fac = paste0("AGEGRP_",relevel(as.factor(AGEGRP),ref="agegrp1")),agegrp_ind = 1,
  ) %>%
  spread(RACE_LABEL_fac,race_ind,fill=0) %>%
  spread(AGEGRP_fac,agegrp_ind,fill=0) %>%
  # add CCI classes
  mutate(
    CCI_CLASS = case_when(
      CCI_SCORE <= 2 & CCI_SCORE >= 1 ~ 'c1',
      CCI_SCORE <= 4 & CCI_SCORE >= 3 ~ 'c2',
      CCI_SCORE >= 5 ~ 'c3',
      TRUE ~ 'c0'
    )
  )

# fill na with 0 for covariate columns
all_colnm<-colnames(dat_proc)
cov_colnm<-c(
  all_colnm[grepl('^(AGEGRP)+',all_colnm)],
  all_colnm[grepl('^(RACE)+',all_colnm)],
  all_colnm[grepl('(_HIST)+?',all_colnm)],
  'SEXF',
  'CCI_SCORE',
  'OBES',
  'T2DM',
  'CKD',
  'COPD',
  'AFIB',
  'HTN',
  'ND',
  'HYSM',
  'INSM',
  'ACG',
  'AHT',
  'ALP',
  'BGR',
  'LIS_DUAL_IND'
)
dat_proc[,cov_colnm][is.na(dat_proc[,cov_colnm])]<-0

# create explicit interaction terms
dat_proc %<>%
  mutate(
    CPAP_x_INSM = CPAP_IND*INSM,
    CPAP_x_HTSM = CPAP_IND*HYSM,
    CPAP_x_COPD = CPAP_IND*COPD,
    CPAP_x_T2DM =  CPAP_IND*T2DM,
    CPAP_x_CKD =  CPAP_IND*CKD,
    CPAP_x_HTN = CPAP_IND*HTN,
    CPAP_x_OBES = CPAP_IND*OBES,
    CPAP_x_AFIB = CPAP_IND*AFIB,
    CPAP_x_ND = CPAP_IND*ND,
    CPAP_x_ACG = CPAP_IND*ACG,
    CPAP_x_AHT = CPAP_IND*AHT,
    CPAP_x_ALP = CPAP_IND*ALP,
    CPAP_x_BGR = CPAP_IND*BGR,
    CPAP_x_MACE = CPAP_IND*MACE_HIST,
    CPAP_x_SES = CPAP_IND*LIS_DUAL_IND
)

saveRDS(dat_proc,file.path(path_to_data_folder,"cpap_exposure_final.rda"))

##==== preprocess cohort for dose-response, adherence analysis ==========
# observation window: before 1-year mark since CPAP initialization
dat_proc<-readRDS(file.path(path_to_data_folder,"cpap_adherence_aset.rda")) %>%
  # survival object
  mutate(
    # MACE survival obj
    MACE_time = coalesce(DAYS_OSA_TO_MACE,pmin(DAYS_OSA_TO_DEATH,DAYS_OSA_TO_CENSOR,na.rm=T)),
    MACE_status = as.numeric(!is.na(DAYS_OSA_TO_MACE)),
    MACE_status_sub = case_when(
      !is.na(DAYS_OSA_TO_DEATH)&is.na(DAYS_OSA_TO_MACE) ~ 2,
      TRUE ~ as.numeric(!is.na(DAYS_OSA_TO_MACE))
    ),
    # MI survival obj
    MI_time = coalesce(DAYS_OSA_TO_MI,pmin(DAYS_OSA_TO_DEATH,DAYS_OSA_TO_CENSOR,na.rm=T)),
    MI_status = as.numeric(!is.na(DAYS_OSA_TO_MI)),
    MI_status_sub = case_when(
      !is.na(DAYS_OSA_TO_DEATH)&is.na(DAYS_OSA_TO_MI) ~ 2,
      TRUE ~ as.numeric(!is.na(DAYS_OSA_TO_MI))
    ),
    # STROKE survival obj
    STROKE_time = coalesce(DAYS_OSA_TO_STROKE,pmin(DAYS_OSA_TO_DEATH,DAYS_OSA_TO_CENSOR,na.rm=T)),
    STROKE_status = as.numeric(!is.na(DAYS_OSA_TO_STROKE)),
    STROKE_status_sub = case_when(
      !is.na(DAYS_OSA_TO_DEATH)&is.na(DAYS_OSA_TO_STROKE) ~ 2,
      TRUE ~ as.numeric(!is.na(DAYS_OSA_TO_STROKE))
    ),
    # HF survival obj
    HF_time = coalesce(DAYS_OSA_TO_HF,pmin(DAYS_OSA_TO_DEATH,DAYS_OSA_TO_CENSOR,na.rm=T)),
    HF_status = as.numeric(!is.na(DAYS_OSA_TO_HF)),
    HF_status_sub = case_when(
      !is.na(DAYS_OSA_TO_DEATH)&is.na(DAYS_OSA_TO_HF) ~ 2,
      TRUE ~ as.numeric(!is.na(DAYS_OSA_TO_HF))
    ),
    # REVASC survival obj
    REVASC_time = coalesce(DAYS_OSA_TO_REVASC,pmin(DAYS_OSA_TO_DEATH,DAYS_OSA_TO_CENSOR,na.rm=T)),
    REVASC_status = as.numeric(!is.na(DAYS_OSA_TO_REVASC)),
    REVASC_status_sub = case_when(
      !is.na(DAYS_OSA_TO_DEATH)&is.na(DAYS_OSA_TO_REVASC) ~ 2,
      TRUE ~ as.numeric(!is.na(DAYS_OSA_TO_REVASC))
    ),
    # all-cause DEATH survival obj
    DEATH_time = coalesce(DAYS_OSA_TO_DEATH,DAYS_OSA_TO_CENSOR),
    DEATH_status = as.numeric(!is.na(DAYS_OSA_TO_DEATH))
  ) %>%
  # pre-MACE event indicator
  mutate(
    MACE_HIST = as.numeric(!is.na(MACE_HIST)),
    MI_HIST = as.numeric(!is.na(MI_HIST)),
    STROKE_HIST = as.numeric(!is.na(STROKE_HIST)),
    HF_HIST = as.numeric(!is.na(HF_HIST)),
    REVASC_HIST = as.numeric(!is.na(REVASC_HIST))
  ) %>%
  mutate(
    # label demographic fields
    RACE_LABEL = case_when(
      HISPANIC=='Y' ~ 'Hispanic',
      RACE=='05' ~ 'White',
      RACE=='03' ~ 'AA',
      RACE=='02' ~ 'Asian',
      RACE=='01' ~ 'AI',
      RACE=='NI' ~ 'Unknown',
      TRUE ~ 'Other'
    ),
    RACE_LABEL = relevel(as.factor(RACE_LABEL),ref="White"),
    RACE_LABEL_fac = paste0("RACE_",relevel(as.factor(RACE_LABEL),ref="White")),race_ind = 1,
    SEXF = as.numeric(SEX=='F'),
    AGEGRP = case_when(AGE_AT_OSA_DX1>=80 ~ 'agegrp4',
                       TRUE ~ paste0('agegrp',floor((AGE_AT_OSA_DX1-65)/5)+1)),
    AGEGRP_fac = paste0("AGEGRP_",relevel(as.factor(AGEGRP),ref="agegrp1")),agegrp_ind = 1,
  ) %>%
  spread(RACE_LABEL_fac,race_ind,fill=0) %>%
  spread(AGEGRP_fac,agegrp_ind,fill=0) %>%
  # add CCI classes
  mutate(
    CCI_CLASS = case_when(
      CCI_SCORE <= 2 & CCI_SCORE >= 1 ~ 'c1',
      CCI_SCORE <= 4 & CCI_SCORE >= 3 ~ 'c2',
      CCI_SCORE >= 5 ~ 'c3',
      TRUE ~ 'c0'
    )
  )

# fill na with 0 for covariate columns
all_colnm<-colnames(dat_proc)
cov_colnm<-c(
  all_colnm[grepl('^(AGEGRP)+',all_colnm)],
  all_colnm[grepl('^(RACE)+',all_colnm)],
  all_colnm[grepl('(_HIST)+?',all_colnm)],
  'SEXF',
  'CCI_SCORE',
  'OBES',
  'T2DM',
  'CKD',
  'COPD',
  'AFIB',
  'HTN',
  'ND',
  'HYSM',
  'INSM',
  'ACG',
  'AHT',
  'ALP',
  'BGR',
  'LIS_DUAL_IND'
)
dat_proc[,cov_colnm][is.na(dat_proc[,cov_colnm])]<-0

# ggplot(dat_proc,aes(x=CPAP_YR1))+
#   geom_histogram(binwidth = 1,color="red",fill="white")+
#   stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-1.5)+
#   scale_x_continuous(breaks = 1:37,labels = 1:37)

# create adherence metrics based on year 1 charges
cut_yr1_med<-quantile(dat_proc$CPAP_YR1, probs = seq(0, 1, length.out = 3))
cut_yr1_tt<-quantile(dat_proc$CPAP_YR1, probs = seq(0, 1, length.out = 4))
cut_yr1_qt<-quantile(dat_proc$CPAP_YR1, probs = seq(0, 1, length.out = 5))
cut_yr1_inc4<-c(seq(min(dat_proc$CPAP_YR1),26,by=4),max(dat_proc$CPAP_YR1))
cut_yr1_inc8<-c(seq(min(dat_proc$CPAP_YR1),26,by=8),max(dat_proc$CPAP_YR1))
dat_proc %<>%
  mutate(
    adherence_yr1_med = paste0('tt_',cut(CPAP_YR1,breaks = cut_yr1_med, right=FALSE, include.lowest=TRUE,labels=FALSE)),
    adherence_yr1_tt = paste0('tt_',cut(CPAP_YR1,breaks = cut_yr1_tt, right=FALSE, include.lowest=TRUE,labels=FALSE)),
    adherence_yr1_tt = relevel(factor(adherence_yr1_tt),ref="tt_1"),
    adherence_yr1_qt = paste0('qt_',cut(CPAP_YR1,breaks = cut_yr1_qt, right=FALSE, include.lowest=TRUE,labels=FALSE)),
    adherence_yr1_qt = relevel(factor(adherence_yr1_qt),ref="qt_1")
  ) %>%
  mutate(
    adherence_yr1_inc4 = paste0('inc4_',cut(CPAP_YR1,breaks = cut_yr1_inc4, right=FALSE, include.lowest=TRUE,labels=FALSE)),
    adherence_yr1_inc4 = relevel(factor(adherence_yr1_inc4),ref="inc4_1"),
    adherence_yr1_inc8 = paste0('inc8_',cut(CPAP_YR1,breaks = cut_yr1_inc8, right=FALSE, include.lowest=TRUE,labels=FALSE)),
    adherence_yr1_inc8 = relevel(factor(adherence_yr1_inc8),ref="inc8_1")
  ) %>%
  mutate(
    adherence_yr1_ww = case_when(CPAP_YR1 > 12 ~ 'adh_yr1_ww_3',
                                 CPAP_YR1 <= 12 & CPAP_YR1 >= 4 ~ 'adh_yr1_ww_2',
                                 TRUE ~ 'adh_yr1_ww_1'),
    adherence_yr1_ww = relevel(factor(adherence_yr1_ww),ref = "adh_yr1_ww_1")
  )

# based on smoothed mixture model?
# mixed<-normalmixEM(dat2_proc$CPAP_YR1,k=2) # cluster 1 biased upward
# plot(mixed,which=2)
# 
# mixed<-spEMsymloc(dat2_proc$CPAP_YR1,mu0=c(4,15)) # very slow
# plot(mixed)

# inductive 
fit_tree <- rpart(
  DEATH_status ~ CPAP_YR1, 
  data = dat_proc, 
  method = "class", 
  minsplit = 2, 
  minbucket = 1, 
  cp = -1
)
pr_tree<-prune(fit_tree,cp=fit_tree$cptable[which.min(fit_tree$cptable[-1,"xerror"])+1,"CP"])
rpart.plot(pr_tree,cex=1)
mort_cut<-c(min(dat_proc$CPAP_YR1),8,15,25,max(dat_proc$CPAP_YR1))
dat_proc %<>% 
  mutate(
     adherence_yr1_mttree = paste0('mttree_',cut(CPAP_YR1,breaks = mort_cut, right=FALSE, include.lowest=TRUE,labels=FALSE))
    ,adherence_yr1_mttree = relevel(factor(adherence_yr1_mttree),ref="mttree_1")
  )

fit_tree2 <- rpart(
  MACE_status ~ CPAP_YR1, 
  data = dat_proc %>% filter(MACE_HIST==0&MACE_time>0), 
  method = "class", 
  minsplit = 2, 
  minbucket = 1, 
  cp = -1
)
pr_tree2<-prune(fit_tree2,cp=fit_tree2$cptable[which.min(fit_tree2$cptable[-1,"xerror"])+1,"CP"])
rpart.plot(pr_tree2,cex=1)
mace_cut<-c(min(dat_proc$CPAP_YR1),12,15,26,max(dat_proc$CPAP_YR1))
dat_proc %<>% 
  mutate(
     adherence_yr1_mctree = paste0('mctree_',cut(CPAP_YR1,breaks = mace_cut, right=FALSE, include.lowest=TRUE,labels=FALSE))
    ,adherence_yr1_mctree = relevel(factor(adherence_yr1_mctree),ref="mctree_1")
  )

# save results
# saveRDS(dat_proc,file.path(path_to_data_folder,"cpap_adherence_final.rda"))
saveRDS(dat_proc,file.path(path_to_data_folder,"cpap_adherence_final2.rda"))

