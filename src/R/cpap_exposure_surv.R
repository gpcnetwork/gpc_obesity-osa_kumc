#################################################################
# Copyright (c) 2021-2025 University of Missouri                   
# Author: Xing Song, xsm7f@umsystem.edu                            
# File: cpap_exposure_surv.R
# Description: model how cpap exposure (binary) affect survival 
#################################################################

rm(list=ls())
pacman::p_load(tidyverse,
               magrittr,
               broom,
               survival,
               survminer,
               kableExtra,
               devtools,
               cmprsk,
               Matrix,
               glmnet,
               scales,
               islasso,
               ranger
)

source_url("https://raw.githubusercontent.com/sxinger/utils/master/sample_util.R")
source_url("https://raw.githubusercontent.com/sxinger/utils/master/analysis_util.R")
source_url("https://raw.githubusercontent.com/sxinger/utils/master/model_util.R")

boots<-2
boots_iter<-1:boots

#==== load data
path_to_data_folder<-file.path(
  getwd(),
  "osa-obesity-cpap-adherence",
  "data"
)
df<-readRDS(file.path(path_to_data_folder,"cpap_exposure_aset.rda"))

#==== matched sampling
path_to_file<-file.path(
  path_to_data_folder,
  paste0("cpap_exposure_surv_adjset_bts",boots,".rda"
  ))

if(!file.exists(path_to_file)){
  adj_set<-matched_sample.ptdm(
    ref_dat = df %>% filter(CPAP_IND==1),
    match_dat = df %>% filter(CPAP_IND==0),
    id_col = 'PATID',
    update_ref = 'DAYS_OSA_TO_CPAP_INIT',
    update_col = 'DEATH_time',
    boots=boots,
    replace = TRUE
  )
  saveRDS(adj_set,file=path_to_file)
}else{
  adj_set<-readRDS(path_to_file)
}

#==== specify pre-selected variables
var_demo_lst<-c(
  "AGEGRP_agegrp2", "AGEGRP_agegrp3","AGEGRP_agegrp4",
  "SEXF",
  "RACE_AA","RACE_AI","RACE_Asian","RACE_Other","RACE_Unknown",
  "LIS_DUAL_IND"  
)

var_comorb_lst<-c(
  "T2DM_HISTORY",
  "HTN_HISTORY",
  "OBESITY_HISTORY",
  "MACE_HISTORY",
  "chronic_obstructive_pulmonary_disease",
  "hypersomnia",
  "insomina"
)

var_interact<-c(
  "CPAP_x_insomnia",
  "CPAP_x_hypersomnia",
  "CPAP_x_obesity",
  "CPAP_x_COPD",
  "CPAP_x_T2DM",
  "CPAP_x_HTN",
  "CPAP_x_MACE",
  "CPAP_x_LowSES"
)

var_outcome<-c(
  "MI_status",
  "STROKE_status",
  "HF_status",
  "REVASC_status",
  "MACE_status",
  "DEATH_status"
)

# repeat experiments over bootstrapped samples
for(boot_i in boots_iter){
  #==== construct matching sample
  df_adj<-adj_set[[boot_i]]$pos %>% select(PATID) %>% 
    left_join(df,by="PATID") %>%
    bind_rows(adj_set[[boot_i]]$neg %>% select(PATID,adj_time) %>%
                left_join(df,by="PATID") %>% mutate(DEATH_time = adj_time) %>% select(-adj_time))
  
  #==== propensity score
  x<-data.matrix(df[,c(var_demo_lst,var_comorb_lst)])
  y<-df$CPAP_IND
  fit_ps<-cv.glmnet(x=x,y=y)
  # plot(fit_ps)
  ps_cpap<-predict(fit_ps, newx = x, s = "lambda.min", type="response")
  ps_cpap_smth<-data.frame(
    PATID=df$PATID,ps_cpap=ps_cpap[,1], y=y
  ) %>% 
    mutate(idx=row_number()) %>%
    mutate(cpap_wt=case_when(y==1 ~ ps_cpap,
                             TRUE ~ 1-ps_cpap)) %>%
    select(PATID,cpap_wt)
  
  #==== decompose
  # out<-islasso(y~x,data=df,lambda=fit_tw$lambda.min)
  # summary(out,pval=0.05)
  
  #===== IPTW-adjusted, Main Effect ===========
  #--create subdir
  path_to_dir<-file.path(path_to_data_folder,
                         "results","surv")
  if(!dir.exists(path_to_dir)) dir.create(path_to_dir)
  #--create subdir/subdir
  path_to_dir<-file.path(path_to_data_folder,
                         "results","surv",
                         paste0("boot",boot_i))
  if(!dir.exists(path_to_dir)) dir.create(path_to_dir)

  #--build IPW model
  path_to_file<-file.path(path_to_dir,"coxph_iptw_main.rda")
  if(!file.exists(path_to_file)){
    df_iptw<-df %>% select(PATID) %>% left_join(ps_cpap_smth,by="PATID")
    fit_frm<-formula(
      paste0("Surv(DEATH_time, DEATH_status) ~ ",
             paste(c(var_demo_lst,var_comorb_lst,"CPAP_IND"),collapse = "+"))
    )
    fit_mort_msm<-coxph(fit_frm, data = df, weights = 1/df_iptw$cpap_wt, model=TRUE)
    # save model - takes about 5min to train a msm model
    saveRDS(fit_mort_msm, file=path_to_file)
  }
  
  #===== IPTW-adjusted, with Interaction ==========
  path_to_file<-file.path(path_to_dir,"coxph_iptw_intx.rda")
  if(file.exists(path_to_file)){
    df_iptw<-df %>% select(PATID) %>% left_join(ps_cpap_smth,by="PATID")
    fit_frm<-formula(
      paste0("Surv(DEATH_time, DEATH_status) ~ ",
             paste(c(var_demo_lst,var_comorb_lst,"CPAP_IND",var_interact),collapse = "+"))
    )
    fit_mort_msm<-coxph(fit_frm, data = df, weights = 1/df_iptw$cpap_wt, model=TRUE)
    # save model - takes about 5min to train a msm model
    saveRDS(fit_mort_msm, file=path_to_file)
  }
  
  #===== Stratified Models ==========
  #-- define strata
  strata<-data.frame(
    val=c(unique(df$AGEGRP),
          unique(df$SEX),
          unique(as.character(df$RACE_LABEL)),
          unique(df$insomina),
          unique(df$hypersomnia),
          unique(df$LIS_DUAL_IND),
          unique(df$OBESITY_HISTORY),
          unique(df$chronic_obstructive_pulmonary_disease),
          unique(df$T2DM_HISTORY),
          unique(df$HTN_HISTORY),
          unique(df$MACE_HISTORY)
    ),
    var=c(rep("AGEGRP",length(unique(df$AGEGRP))),
          rep("SEX",length(unique(df$SEX))),
          rep("RACE_LABEL",length(unique(df$RACE_LABEL))),
          rep("insomina",length(unique(df$insomina))),
          rep("hypersomnia",length(unique(df$hypersomnia))),
          rep("LIS_DUAL_IND",length(unique(df$LIS_DUAL_IND))),
          rep("OBESITY_HISTORY",length(unique(df$OBESITY_HISTORY))),
          rep("chronic_obstructive_pulmonary_disease",length(unique(df$chronic_obstructive_pulmonary_disease))),
          rep("T2DM_HISTORY",length(unique(df$T2DM_HISTORY))),
          rep("HTN_HISTORY",length(unique(df$HTN_HISTORY))),
          rep("MACE_HISTORY",length(unique(df$MACE_HISTORY)))
    ),
    excld=c(rep("AGE",length(unique(df$AGEGRP))),
            rep("SEXF",length(unique(df$SEX))),
            rep("RACE",length(unique(df$RACE_LABEL))),
            rep("insomina",length(unique(df$insomina))),
            rep("hypersomnia",length(unique(df$hypersomnia))),
            rep("(LIS_DUAL_IND|LowSES)",length(unique(df$LIS_DUAL_IND))),
            rep("OBESITY",length(unique(df$OBESITY_HISTORY))),
            rep("(chronic_obstructive_pulmonary_disease|COPD)",length(unique(df$chronic_obstructive_pulmonary_disease))),
            rep("T2DM",length(unique(df$T2DM_HISTORY))),
            rep("HTN",length(unique(df$HTN_HISTORY))),
            rep("MACE",length(unique(df$MACE_HISTORY)))
    )
  )
  
  #==== IPW-adjusted, Main Effects ==========
  #-- built IPW-adjusted model within each strata
  result<-c()
  path_to_file<-file.path(path_to_dir,"coxph_strata_main.csv")
  if(!file.exists(path_to_file)){
    for(i in seq_len(nrow(strata))){
      #- model training
      var_filter<-c(var_demo_lst,var_comorb_lst)[!grepl(strata$excld[i],c(var_demo_lst,var_comorb_lst))]
      fit_frm<-formula(paste0("Surv(DEATH_time, DEATH_status) ~ ",
                              paste(c(var_filter,"CPAP_IND"),collapse = "+")))
      df_str<-df %>% filter(.data[[strata$var[i]]]==strata$val[i])
      df_iptw_str<-df_str %>% select(PATID) %>% left_join(ps_cpap_smth,by="PATID")
      fit_mort_cov<-coxph(fit_frm,data = df_str, weights = 1/df_iptw_str$cpap_wt)
      
      #- get all coefficients
      fit_summ<-summary(fit_mort_cov)$coefficients
      fit_var<-rownames(fit_summ)
      rownames(fit_summ)<-NULL
      result<-rbind(result,
                    cbind(stratum_var=strata$var[i],
                          stratum_val=strata$val[i],
                          fit_var=fit_var,
                          fit_summ))
      # print(paste0("stratified by:",strata$var[i],strata$val[i]))
    }
  }
  write.csv(result,
            file=file.path(path_to_dir,"coxph_strata_main.csv"),
            row.names = F)
  
  #==== IPW-adjusted, with interactions ==========
  #-- built IPW-adjusted model within each strata
  result<-c()
  path_to_file<-file.path(path_to_dir,"coxph_strata_main.csv")
  if(!file.exists(path_to_file)){
    for(i in seq_len(nrow(strata))){
      #- model training
      var_filter<-c(var_demo_lst,var_comorb_lst,var_interact)[!grepl(strata$excld[i],c(var_demo_lst,var_comorb_lst,var_interact))]
      fit_frm<-formula(paste0("Surv(DEATH_time, DEATH_status) ~ ",
                              paste(c(var_filter,"CPAP_IND"),collapse = "+")))
      df_str<-df %>% filter(.data[[strata$var[i]]]==strata$val[i])
      df_iptw_str<-df_str %>% select(PATID) %>% left_join(ps_cpap_smth,by="PATID")
      fit_mort_cov<-coxph(fit_frm,data = df_str, weights = 1/df_iptw_str$cpap_wt)
      
      #- get all coefficients
      fit_summ<-summary(fit_mort_cov)$coefficients
      fit_var<-rownames(fit_summ)
      rownames(fit_summ)<-NULL
      result<-rbind(result,
                    cbind(stratum_var=strata$var[i],
                          stratum_val=strata$val[i],
                          fit_var=fit_var,
                          fit_summ))
      # print(paste0("stratified by:",strata$var[i],strata$val[i]))
    }
  }
  write.csv(result,
            file=file.path(path_to_dir,"coxph_strata_intx.csv"),
            row.names = F)
}





