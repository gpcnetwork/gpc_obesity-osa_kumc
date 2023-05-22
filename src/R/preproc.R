#################################################################
# Copyright (c) 2021-2025 University of Missouri                   
# Author: Xing Song, xsm7f@umsystem.edu                            
# File: preproc.R
# Description: data preprocess to form final analytic set with
#              one patient per row
#################################################################

rm(list=ls()); gc()
setwd("C:/repo/GPC-Analytics-Weight-Cohort")

# install.packages("pacman")
pacman::p_load(
  tidyverse,
  magrittr,
  devtools,
  flexmix,
  mixtools,
  rpart,rpart.plot
)

path_to_data_folder<-file.path(
  getwd(),
  "osa-obesity-cpap-adherence",
  "data"
)

##==== preprocess cohort for exposure analysis ========
# observation window: covariates no later than first OSA diagnosis
dat_cov_hist<-readRDS(file.path(path_to_data_folder,"cpap_exposure_cov.rds")) %>%
  filter(DAYS_OSA_TO_COV <= 0) %>% 
  group_by(PATID) %>% 
  arrange(DAYS_OSA_TO_COV) %>% slice(1:1) %>% ungroup %>%
  mutate(DAYS_OSA_TO_COV = 1) %>%
  spread(COV,DAYS_OSA_TO_COV, fill = 0)

# routine data preprocess
dat_proc<-readRDS(file.path(path_to_data_folder,"cpap_exposure_tbl1.rds")) %>%
  mutate(
    CPAP_IND = as.numeric(!is.na(DAYS_OSA_TO_CPAP_INIT)), 
    # MACE survival obj
    MACE_time = coalesce(DAYS_OSA_TO_MACE,DAYS_OSA_TO_CENSOR),
    MACE_status = as.numeric(!is.na(DAYS_OSA_TO_MACE)),
    # MI survival obj
    MI_time = coalesce(DAYS_OSA_TO_MI,DAYS_OSA_TO_CENSOR),
    MI_status = as.numeric(!is.na(DAYS_OSA_TO_MI)),
    # STROKE survival obj
    STROKE_time = coalesce(DAYS_OSA_TO_STROKE,DAYS_OSA_TO_CENSOR),
    STROKE_status = as.numeric(!is.na(DAYS_OSA_TO_STROKE)),
    # HF survival obj
    HF_time = coalesce(DAYS_OSA_TO_HF,DAYS_OSA_TO_CENSOR),
    HF_status = as.numeric(!is.na(DAYS_OSA_TO_HF)),
    # REVASC survival obj
    REVASC_time = coalesce(DAYS_OSA_TO_REVASC,DAYS_OSA_TO_CENSOR),
    REVASC_status = as.numeric(!is.na(DAYS_OSA_TO_REVASC)),
    # all-cause DEATH survival obj
    DEATH_time = coalesce(DAYS_OSA_TO_DEATH,DAYS_OSA_TO_CENSOR),
    DEATH_status = as.numeric(!is.na(DAYS_OSA_TO_DEATH)),
    # label demographic fields
    RACE_LABEL = case_when(RACE=='05' ~ 'White',
                           RACE=='03' ~ 'AA',
                           RACE=='02' ~ 'Asian',
                           RACE=='01' ~ 'AI',
                           RACE=='OT' & HISPANIC=='Y' ~ 'Hispanic',
                           RACE=='NI' ~ 'Unknown',
                           TRUE ~ 'Other'),
    RACE_LABEL = relevel(as.factor(RACE_LABEL),ref="White"),
    RACE_LABEL_fac = paste0("RACE_",relevel(as.factor(RACE_LABEL),ref="White")),race_ind = 1,
    SEXF = as.numeric(SEX=='F'),
    AGEGRP = case_when(AGE_AT_OSA_DX1>=80 ~ 'agegrp4',
                       TRUE ~ paste0('agegrp',floor((AGE_AT_OSA_DX1-65)/5)+1)),
    AGEGRP_fac = paste0("AGEGRP_",relevel(as.factor(AGEGRP),ref="agegrp1")),agegrp_ind = 1,
  ) %>%
  spread(RACE_LABEL_fac,race_ind,fill=0) %>%
  spread(AGEGRP_fac,agegrp_ind,fill=0) %>%
  left_join(dat_cov_hist,by="PATID")

# fill na with 0 for covariate columns
cov_colnm<-colnames(dat_cov_hist)[!grepl("PATID",colnames(dat_cov_hist))]
dat_proc[,cov_colnm][is.na(dat_proc[,cov_colnm])]<-0

# create explicit interaction terms
dat_proc %<>%
  mutate(
    CPAP_x_insomnia = CPAP_IND*insomina,
    CPAP_x_hypersomnia = CPAP_IND*hypersomnia,
    CPAP_x_obesity = CPAP_IND*OBESITY_HISTORY,
    CPAP_x_COPD = CPAP_IND*`chronic-obstructive-pulmonary-disease`,
    CPAP_x_T2DM =  CPAP_IND*T2DM_HISTORY,
    CPAP_x_HTN = CPAP_IND*HTN_HISTORY,
    CPAP_x_MACE = CPAP_IND*MACE_HISTORY,
    CPAP_x_LowSES = CPAP_IND*LIS_DUAL_IND
)

colnames(dat_proc)<-gsub("-","_",colnames(dat_proc)) # so the column name is easier to work with in R
saveRDS(dat_proc,file.path(path_to_data_folder,"cpap_exposure_aset.rda"))

##==== preprocess cohort for dose-response, adherence analysis ==========
dose_resp_tbl1<-readRDS(file.path(path_to_data_folder,"cpap_dose_response_tbl1.rds")) %>%
  filter(coalesce(DAYS_CPAP_INIT_TO_DEATH,DAYS_CPAP_INIT_TO_CENSOR) > 365)

# observation window: before 1-year mark since CPAP initialization
dat2_cov_hist<-readRDS(file.path(path_to_data_folder,"cpap_dose_response_cov.rds")) %>%
  left_join(dose_resp_tbl1 %>% select(PATID,DAYS_OSA_TO_CPAP_INIT),by=c("PATID"),multiple = "all") %>%
  mutate(DAYS_CPAP_INIT_TO_COV = DAYS_OSA_TO_COV - DAYS_OSA_TO_CPAP_INIT) %>%
  filter(DAYS_CPAP_INIT_TO_COV+365 <= 0) %>% 
  select(-DAYS_OSA_TO_COV) %>%
  group_by(PATID) %>% 
  arrange(DAYS_CPAP_INIT_TO_COV) %>% slice(1:1) %>% ungroup %>%
  mutate(DAYS_CPAP_INIT_TO_COV = 1) %>%
  spread(COV,DAYS_CPAP_INIT_TO_COV, fill = 0)

dat2_proc<-dose_resp_tbl1 %>%
  mutate(
    # MACE survival obj
    MACE_time = coalesce(DAYS_CPAP_INIT_TO_MACE,DAYS_CPAP_INIT_TO_CENSOR) - 365,
    MACE_status = as.numeric(!is.na(DAYS_CPAP_INIT_TO_MACE)),
    # MI survival obj
    MI_time = coalesce(DAYS_CPAP_INIT_TO_MI,DAYS_CPAP_INIT_TO_CENSOR) - 365,
    MI_status = as.numeric(!is.na(DAYS_CPAP_INIT_TO_MI)),
    # STROKE survival obj
    STROKE_time = coalesce(DAYS_CPAP_INIT_TO_STROKE,DAYS_CPAP_INIT_TO_CENSOR) - 365,
    STROKE_status = as.numeric(!is.na(DAYS_CPAP_INIT_TO_STROKE)),
    # HF survival obj
    HF_time = coalesce(DAYS_CPAP_INIT_TO_HF,DAYS_CPAP_INIT_TO_CENSOR) - 365,
    HF_status = as.numeric(!is.na(DAYS_CPAP_INIT_TO_HF)),
    # REVASC survival obj
    REVASC_time = coalesce(DAYS_CPAP_INIT_TO_REVASC,DAYS_CPAP_INIT_TO_CENSOR) - 365,
    REVASC_status = as.numeric(!is.na(DAYS_CPAP_INIT_TO_REVASC)),
    # all-cause DEATH survival obj
    DEATH_time = coalesce(DAYS_CPAP_INIT_TO_DEATH,DAYS_CPAP_INIT_TO_CENSOR) - 365,
    DEATH_status = as.numeric(!is.na(DAYS_CPAP_INIT_TO_DEATH)),
    # label demographic fields
    RACE_LABEL = case_when(RACE=='05' ~ 'White',
                           RACE=='03' ~ 'AA',
                           RACE=='02' ~ 'Asian',
                           RACE=='01' ~ 'AI',
                           RACE=='OT' & HISPANIC=='Y' ~ 'Hispanic',
                           RACE=='NI' ~ 'Unknown',
                           TRUE ~ 'Other'),
    RACE_LABEL = relevel(as.factor(RACE_LABEL),ref="White"),
    RACE_LABEL_fac = paste0("RACE_",relevel(as.factor(RACE_LABEL),ref="White")),race_ind = 1,
    SEXF = as.numeric(SEX=='F'),
    AGEGRP = case_when(AGE_AT_OSA_DX1>=80 ~ 'agegrp4',
                       TRUE ~ paste0('agegrp',floor((AGE_AT_OSA_DX1-65)/5)+1)),
    AGEGRP_fac = paste0("AGEGRP_",relevel(as.factor(AGEGRP),ref="agegrp1")),agegrp_ind = 1
  ) %>%
  spread(RACE_LABEL_fac,race_ind,fill=0) %>%
  spread(AGEGRP_fac,agegrp_ind,fill=0) %>%
  # gather cpap supply annual summary
  group_by(PATID) %>%
  mutate(
    cpap_total = length(DAYS_CPAP_INIT_TO_SUPPLY)+1,
    cpap_total_aft_yr1 = sum(DAYS_CPAP_INIT_TO_SUPPLY>365),
    cpap_yr1 = sum((DAYS_CPAP_INIT_TO_SUPPLY<=365))+1,
    cpap_yr2 = sum((DAYS_CPAP_INIT_TO_SUPPLY>365 & DAYS_CPAP_INIT_TO_SUPPLY<=365*2)),
    cpap_yr3 = sum((DAYS_CPAP_INIT_TO_SUPPLY>365*2 & DAYS_CPAP_INIT_TO_SUPPLY<=365*3)),
    cpap_yr4 = sum((DAYS_CPAP_INIT_TO_SUPPLY>365*3 & DAYS_CPAP_INIT_TO_SUPPLY<=365*4)),
    cpap_yr5 = sum((DAYS_CPAP_INIT_TO_SUPPLY>365*4 & DAYS_CPAP_INIT_TO_SUPPLY<=365*5)),
    cpap_yr6 = sum((DAYS_CPAP_INIT_TO_SUPPLY>365*5 & DAYS_CPAP_INIT_TO_SUPPLY<=365*6)),
    cpap_yr7 = sum((DAYS_CPAP_INIT_TO_SUPPLY>365*6 & DAYS_CPAP_INIT_TO_SUPPLY<=365*7)),
    cpap_yrs = pmax(round(max(DAYS_CPAP_INIT_TO_SUPPLY)/365.25),1),
    cpap_per_yr = cpap_total/cpap_yrs,
    cpap_yrs_aft_yr1 = pmax(round(max(DAYS_CPAP_INIT_TO_SUPPLY)/365.25-1),1),
    cpap_per_yr_aft_yr1 = cpap_total_aft_yr1/cpap_yrs_aft_yr1
  ) %>%
  ungroup

# # calculate adherence metrics
# #-- Wickwire(2021): total cpap charges
# dat2_proc %<>%
#   mutate(
#     adherence_ww = case_when(cpap_total > 12 ~ 'adh_ww_high',
#                              cpap_total <= 12 & cpap_total >= 4 ~ 'adh_ww_partial',
#                              TRUE ~ 'adh_ww_low'),
#     adherence_ww = relevel(factor(adherence_ww),ref = "adh_ww_low")
#   )
#          
# #-- empirical1: overall total, tertile, quartile
# cut_all_tt<-quantile(dat2_proc$cpap_total[dat2_proc$cpap_total>0], probs = seq(0, 1, length.out = 4))
# cut_all_qt<-quantile(dat2_proc$cpap_total[dat2_proc$cpap_total>0], probs = seq(0, 1, length.out = 5))
# dat2_proc %<>%
#   mutate(
#     adherence_all_tt = paste0('adh_all_tt_',cut(cpap_total,breaks = cut_all_tt, include.lowest=TRUE,labels=FALSE)),
#     adherence_all_qt = paste0('adh_all_qt_',cut(cpap_total,breaks = cut_all_qt, include.lowest=TRUE,labels=FALSE))
#   )
# 
# #-- empirical2: year1 total, tertile, quartile
# cut_yr1_tt<-quantile(dat2_proc$cpap_yr1[dat2_proc$cpap_yr1>0], probs = seq(0, 1, length.out = 4))
# cut_yr1_qt<-quantile(dat2_proc$cpap_yr1[dat2_proc$cpap_yr1>0], probs = seq(0, 1, length.out = 5))
# dat2_proc %<>%
#   mutate(
#     adherence_yr1_tt = paste0('adh_yr1_tt_',cut(cpap_yr1,breaks = cut_yr1_tt, right=FALSE, include.lowest=TRUE,labels=FALSE)),
#     adherence_yr1_qt = paste0('adh_yr1_qt_',cut(cpap_yr1,breaks = cut_yr1_qt, right=FALSE, include.lowest=TRUE,labels=FALSE))
#   ) %>%
#   mutate(
#     adherence_yr1_ww = case_when(cpap_yr1 > 12 ~ 'adh_yr1_ww_high',
#                                  cpap_yr1 <= 12 & cpap_yr1 >= 4 ~ 'adh_yr1_ww_partial',
#                                  TRUE ~ 'adh_yr1_ww_low'),
#     adherence_yr1_ww = relevel(factor(adherence_yr1_ww),ref = "adh_yr1_ww_low"),
#     adherence_yr1_emp = case_when(cpap_yr1 > 15 ~ 'adh_yr1_emp_high',
#                                  cpap_yr1 <= 15 & cpap_yr1 > 4 ~ 'adh_yr1_emp_partial',
#                                  TRUE ~ 'adh_yr1_emp_low'),
#     adherence_yr1_emp = relevel(factor(adherence_yr1_emp),ref = "adh_yr1_emp_low")
#   )
# 
# #-- empirical3: yearly, tertile, quartile
# cut_peryr_tt<-quantile(dat2_proc$cpap_per_yr[dat2_proc$cpap_per_yr>0], probs = seq(0, 1, length.out = 4))
# cut_peryr_qt<-quantile(dat2_proc$cpap_per_yr[dat2_proc$cpap_per_yr>0], probs = seq(0, 1, length.out = 5))
# dat2_proc %<>%
#   mutate(
#     adherence_peryr_tt = case_when(cpap_per_yr>0 ~ paste0('adh_peryr_tt_',cut(cpap_per_yr,breaks = cut_peryr_tt, include.lowest=TRUE,labels=FALSE)),
#                                    TRUE ~ 'adh_peryr_tt_0'),
#     adherence_peryr_qt = case_when(cpap_per_yr>0 ~ paste0('adh_peryr_qt_',cut(cpap_per_yr,breaks = cut_peryr_qt, include.lowest=TRUE,labels=FALSE)),
#                                    TRUE ~ 'adh_peryr_qt_0')
#   )
# 
# #-- empirical4: subsequent yearly, tertile, quartile
# cut_subperyr_tt<-quantile(dat2_proc$cpap_per_yr_aft_yr1[dat2_proc$cpap_per_yr_aft_yr1>0], probs = seq(0, 1, length.out = 4))
# cut_subperyr_qt<-quantile(dat2_proc$cpap_per_yr_aft_yr1[dat2_proc$cpap_per_yr_aft_yr1>0], probs = seq(0, 1, length.out = 5))
# dat2_proc %<>%
#   mutate(
#     adherence_subperyr_tt = case_when(cpap_per_yr_aft_yr1>0 ~ paste0('adh_subperyr_tt_',cut(cpap_per_yr_aft_yr1,breaks = cut_subperyr_tt, include.lowest=TRUE,labels=FALSE)),
#                                    TRUE ~ 'adh_subperyr_tt_0'),
#     adherence_subperyr_qt = case_when(cpap_per_yr_aft_yr1>0 ~ paste0('adh_subperyr_qt_',cut(cpap_per_yr_aft_yr1,breaks = cut_subperyr_qt, include.lowest=TRUE,labels=FALSE)),
#                                    TRUE ~ 'adh_subperyr_qt_0')
#   )

# only use year 1 adherence metrics for outcome model
dat2_proc %<>%
  filter(cpap_yr1 > 0) %>%
  mutate(
    MACE_HISTORY = pmax(MACE_HISTORY,as.numeric(DAYS_CPAP_INIT_TO_MACE <= 365)),
    MI_HISTORY = pmax(MI_HISTORY,as.numeric(DAYS_CPAP_INIT_TO_MI <= 365)),
    STROKE_HISTORY = pmax(STROKE_HISTORY,as.numeric(DAYS_CPAP_INIT_TO_STROKE <= 365)),
    HF_HISTORY = pmax(HF_HISTORY,as.numeric(DAYS_CPAP_INIT_TO_HF <= 365)),
    REVASC_HISTORY = pmax(REVASC_HISTORY,as.numeric(DAYS_CPAP_INIT_TO_REVASC <= 365)),
    T2DM_HISTORY = pmax(T2DM_HISTORY,as.numeric(T2DM_AFT_MIN < DAYS_OSA_TO_CPAP_INIT + 365)),
    HTN_HISTORY = pmax(HTN_HISTORY,as.numeric(HTN_AFT_MIN < DAYS_OSA_TO_CPAP_INIT + 365)),
    OBESITY_HISTORY = pmax(OBESITY_HISTORY,as.numeric(OBESITY_AFT_MIN < DAYS_OSA_TO_CPAP_INIT + 365))
  ) %>%
  select(
    PATID,DAYS_OSA_TO_CPAP_INIT,
    AGE_AT_OSA_DX1, SEX, SEXF, RACE_LABEL, AGEGRP, LIS_DUAL_IND,
    AGEGRP_agegrp1,AGEGRP_agegrp2,AGEGRP_agegrp3,AGEGRP_agegrp4,
    RACE_AA,RACE_AI,RACE_Asian,RACE_Other,RACE_Unknown,RACE_White,
    MACE_HISTORY,MI_HISTORY,STROKE_HISTORY,HF_HISTORY,REVASC_HISTORY,
    T2DM_HISTORY, HTN_HISTORY,OBESITY_HISTORY,
    cpap_yr1,
    MACE_time, MACE_status,
    MI_time,MI_status,
    STROKE_time, STROKE_status,
    HF_time, HF_status,
    REVASC_time, REVASC_status,
    DEATH_time, DEATH_status
  ) %>% unique %>%
  left_join(dat2_cov_hist,
            by=c("PATID","DAYS_OSA_TO_CPAP_INIT"))

# fill na with 0 for covariate columns
dat2_proc[is.na(dat2_proc)]<-0

# ggplot(dat2_proc,aes(x=cpap_yr1))+
#   geom_histogram(binwidth = 1,color="red",fill="white")+
#   stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-1.5)+ 
#   scale_x_continuous(breaks = 1:37,labels = 1:37)

# create adherence metrics based on year 1 charges
cut_yr1_med<-quantile(dat2_proc$cpap_yr1, probs = seq(0, 1, length.out = 3))
cut_yr1_tt<-quantile(dat2_proc$cpap_yr1, probs = seq(0, 1, length.out = 4))
cut_yr1_qt<-quantile(dat2_proc$cpap_yr1, probs = seq(0, 1, length.out = 5))
cut_yr1_inc4<-c(seq(min(dat2_proc$cpap_yr1),26,by=4),max(dat2_proc$cpap_yr1))
cut_yr1_inc8<-c(seq(min(dat2_proc$cpap_yr1),26,by=8),max(dat2_proc$cpap_yr1))
dat2_proc %<>%
  mutate(
    adherence_yr1_med = paste0('tt_',cut(cpap_yr1,breaks = cut_yr1_med, right=FALSE, include.lowest=TRUE,labels=FALSE)),
    adherence_yr1_tt = paste0('tt_',cut(cpap_yr1,breaks = cut_yr1_tt, right=FALSE, include.lowest=TRUE,labels=FALSE)),
    adherence_yr1_tt = relevel(factor(adherence_yr1_tt),ref="tt_1"),
    adherence_yr1_qt = paste0('qt_',cut(cpap_yr1,breaks = cut_yr1_qt, right=FALSE, include.lowest=TRUE,labels=FALSE)),
    adherence_yr1_qt = relevel(factor(adherence_yr1_qt),ref="qt_1")
  ) %>%
  mutate(
    adherence_yr1_inc4 = paste0('inc4_',cut(cpap_yr1,breaks = cut_yr1_inc4, right=FALSE, include.lowest=TRUE,labels=FALSE)),
    adherence_yr1_inc4 = relevel(factor(adherence_yr1_inc4),ref="inc4_1"),
    adherence_yr1_inc8 = paste0('inc8_',cut(cpap_yr1,breaks = cut_yr1_inc8, right=FALSE, include.lowest=TRUE,labels=FALSE)),
    adherence_yr1_inc8 = relevel(factor(adherence_yr1_inc8),ref="inc8_1")
  ) %>%
  mutate(
    adherence_yr1_ww = case_when(cpap_yr1 > 12 ~ 'adh_yr1_ww_3',
                                 cpap_yr1 <= 12 & cpap_yr1 >= 4 ~ 'adh_yr1_ww_2',
                                 TRUE ~ 'adh_yr1_ww_1'),
    adherence_yr1_ww = relevel(factor(adherence_yr1_ww),ref = "adh_yr1_ww_1")
  )

# based on smoothed mixture model?
# mixed<-normalmixEM(dat2_proc$cpap_yr1,k=2) # cluster 1 biased upward
# plot(mixed,which=2)
# 
# mixed<-spEMsymloc(dat2_proc$cpap_yr1,mu0=c(4,15)) # very slow
# plot(mixed)

# inductive 
fit_tree <- rpart(
  DEATH_status ~ cpap_yr1, 
  data = dat2_proc, 
  method = "class", 
  minsplit = 2, 
  minbucket = 1, 
  cp = -1
)
pr_tree<-prune(fit_tree,cp=fit_tree$cptable[which.min(fit_tree$cptable[-1,"xerror"])+1,"CP"])
rpart.plot(pr_tree)
mort_cut<-c(min(dat2_proc$cpap_yr1),12,16,max(dat2_proc$cpap_yr1))
dat2_proc %<>% 
  mutate(
    adherence_yr1_mttree = paste0('mttree_',cut(cpap_yr1,breaks = mort_cut, right=FALSE, include.lowest=TRUE,labels=FALSE))
    ,adherence_yr1_mttree = relevel(factor(adherence_yr1_mttree),ref="mttree_1")
  )


fit_tree2 <- rpart(
  MACE_status ~ cpap_yr1, 
  data = dat2_proc %>% filter(MACE_HISTORY==0&MACE_time>0), 
  method = "class", 
  minsplit = 2, 
  minbucket = 1, 
  cp = -1
)
pr_tree2<-prune(fit_tree2,cp=fit_tree2$cptable[which.min(fit_tree2$cptable[-1,"xerror"])+1,"CP"])
rpart.plot(pr_tree2)
mace_cut<-c(min(dat2_proc$cpap_yr1),15,17,24,max(dat2_proc$cpap_yr1))
dat2_proc %<>% 
  mutate(
    adherence_yr1_mctree = paste0('mctree_',cut(cpap_yr1,breaks = mace_cut, right=FALSE, include.lowest=TRUE,labels=FALSE))
    ,adherence_yr1_mctree = relevel(factor(adherence_yr1_mctree),ref="mctree_1")
  )

# save results
colnames(dat2_proc)<-gsub("-","_",colnames(dat2_proc)) # so the column name is easier to work with in R
saveRDS(dat2_proc,file.path(path_to_data_folder,"cpap_dose_response_aset.rda"))

# clean again
rm(list=ls()); gc()
