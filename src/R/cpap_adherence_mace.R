#################################################################
# Copyright (c) 2021-2025 University of Missouri                   
# Author: Xing Song, xsm7f@umsystem.edu                            
# File: cpap_exposure_surv.R
# Description: model how cpap exposure (binary) affect survival 
#################################################################

rm(list=ls())
knitr::opts_chunk$set(echo = FALSE,
                      message = FALSE,
                      warning = FALSE,
                      fig.height=5)

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
               glmnetcr,
               inflection
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
df<-readRDS(file.path(path_to_data_folder,"cpap_dose_response_aset.rda")) %>%
  filter(MACE_HISTORY == 0&MACE_time>0)

#==== specify pre-selected variables
var_demo_lst<-c(
  "AGEGRP_agegrp2", "AGEGRP_agegrp3","AGEGRP_agegrp4", # ref = "AGEGRP_agegrp4"
  "SEXF",
  "RACE_AA","RACE_AI","RACE_Asian","RACE_Other","RACE_Unknown", # ref = "RACE_White"
  "LIS_DUAL_IND"  
)

var_comorb_lst<-c(
  "T2DM_HISTORY",
  "HTN_HISTORY",
  "OBESITY_HISTORY",
  "chronic_obstructive_pulmonary_disease",
  "hypersomnia",
  "insomina"
)

adh_metric<-c(
   "cpap_yr1"
  ,"adherence_yr1_med"
  ,"adherence_yr1_tt"
  ,"adherence_yr1_qt"
  ,"adherence_yr1_inc4"
  ,"adherence_yr1_inc8"
  ,"adherence_yr1_ww"
  ,"adherence_yr1_mctree"
)

mace_endpts<-c(
  'MACE',
  'HF',
  'MI',
  'STROKE',
  'REVASC'
)

for(mace_endpt in mace_endpts){
  for(adh in adh_metric){
    for(boot_i in boots_iter){
      # mace_endpt<-mace_endpts[1]
      # adh<-adh_metric[1]
      # boot_i<-1
      #==== propensity score
      df_ps<-df[,c(adh,var_demo_lst,var_comorb_lst)] %>% 
        filter(!is.na(get(adh)))
      x<-data.matrix(df_ps[,-1])
      y<-unlist(df_ps[,1]) 
      if(adh == 'cpap_yr1'){
        fit_ps<-cv.glmnet(x=x,y=y,family = "poisson")
        # plot(fit_ps)
        ps_cpap<-predict(fit_ps, newx = x, s = "lambda.min", type="response")
        ps_cpap_smth<-data.frame(
          PATID=df$PATID,
          lambda=ps_cpap[,1],
          k = y
        ) %>%
          mutate(
            explambda = exp(-lambda),
            lambdak=lambda^k,
            kfac=factorial(k),
            ps_wt = explambda*lambdak/kfac,
            cpap_wt=rescale(ps_wt,to=c(0.01,0.99))
            )
      }else{
        fit_ps<-cv.glmnet(x=x,y=y,family = "multinomial", type.multinomial = "grouped")
        # plot(fit_ps)
        ps_cpap<-predict(fit_ps, newx = x, s = "lambda.min", type="response")
        ps_cpap_smth<-data.frame(
          PATID=df$PATID,
          lambda=ps_cpap[,,1],
          actual = y
        ) %>%
          gather(pred,ps_wt,-actual,-PATID) %>%
          mutate(pred = gsub("lambda\\.","",pred)) %>%
          filter(pred==actual) %>%
          mutate(ps_wt = as.numeric(ps_wt),
                 cpap_wt=rescale(ps_wt,to=c(0.01,0.99))
          )
      }
      
      #==== decompose
      # out<-islasso(y~x,data=df,lambda=fit_tw$lambda.min)
      # summary(out,pval=0.05)
      
      #===== IPTW-adjusted, Main Effect ===========
      #--create subdir
      path_to_dir<-file.path(path_to_data_folder,
                             "results",paste0("adh_",tolower(mace_endpt)))
      if(!dir.exists(path_to_dir)) dir.create(path_to_dir)
      #--create subdir/subdir
      path_to_dir<-file.path(path_to_dir,adh)
      if(!dir.exists(path_to_dir)) dir.create(path_to_dir)
      #--create subdir/subdir
      path_to_dir<-file.path(path_to_dir,paste0("boot",boot_i))
      if(!dir.exists(path_to_dir)) dir.create(path_to_dir)
      
      #--build IPW model
      path_to_file<-file.path(path_to_dir,"coxph_iptw_main_acm.rda")
      df_iptw<-df %>% select(PATID) %>% left_join(ps_cpap_smth,by="PATID")
      if(!file.exists(path_to_file)){
        if(adh == 'cpap_yr1'){
          fit_frm<-formula(
            paste0("Surv(",mace_endpt,"_time, ",mace_endpt,"_status) ~ ",
                   paste(c(var_demo_lst,var_comorb_lst),collapse = "+"),
                   "+ pspline(",adh,", df = 4)")
          )
          fit_mort_msm<-coxph(fit_frm, data = df, weights = 1/df_iptw$cpap_wt)
          # tp<-termplot(fit_mort_msm, term=18, se=TRUE, col.term=1, col.se=1,plot=F)
        }else{
          fit_frm<-formula(
            paste0("Surv(",mace_endpt,"_time, ",mace_endpt,"_status) ~ ",
                   paste(c(var_demo_lst,var_comorb_lst,adh),collapse = "+"))
          )
          fit_mort_msm<-coxph(fit_frm, data = df, weights = 1/df_iptw$cpap_wt)
          # summary(fit_mort_msm)
        }
        #----------------------------------
        print(paste0(mace_endpt,"...",adh,"...",boot_i,"...","IPTW model developed"))
        saveRDS(fit_mort_msm, file=path_to_file)
      }else{
        #--post-hoc: extract nonlinear termplot
        if(adh == 'cpap_yr1'){
          if(!file.exists(file.path(path_to_dir,"pspline_termplot.rda"))){
            fit_mort_msm<-readRDS(path_to_file)
            tp<-termplot(fit_mort_msm, term=18, se=TRUE, col.term=1, col.se=1,plot=F)
            #----------------------------------
            print(paste0(mace_endpt,"...",adh,"...",boot_i,"...","spline model extracted"))
            saveRDS(tp,file=file.path(path_to_dir,"pspline_termplot.rda"))
          }
        }
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
      path_to_file<-file.path(path_to_dir,"coxph_strata_main_acm.csv")
      if(!file.exists(path_to_file)){
        result<-list()
        for(i in seq_len(nrow(strata))){
          #- model training
          var_filter<-c(var_demo_lst,var_comorb_lst)[!grepl(strata$excld[i],c(var_demo_lst,var_comorb_lst))]
          if(adh == 'cpap_yr1'){
            fit_frm<-formula(
              paste0("Surv(",mace_endpt,"_time, ",mace_endpt,"_status) ~ ",
                     paste(var_filter,collapse = "+"),
                     "+ pspline(",adh,", df = 4)")
            )
            fit_mort_msm<-coxph(fit_frm, data = df, weights = 1/df_iptw$cpap_wt)
            # tp<-termplot(fit_mort_msm, term=18, se=TRUE, col.term=1, col.se=1,plot=F)
          }else{
            fit_frm<-formula(
              paste0("Surv(",mace_endpt,"_time, ",mace_endpt,"_status) ~ ",
                     paste(c(var_filter,adh),collapse = "+"))
            )
            fit_mort_msm<-coxph(fit_frm, data = df, weights = 1/df_iptw$cpap_wt)
            # ggforest(fit_mort_msm)
          }
          df_str<-df %>% filter(.data[[strata$var[i]]]==strata$val[i])
          df_iptw_str<-df_str %>% select(PATID) %>% left_join(ps_cpap_smth,by="PATID")
          fit_mort_cov<-coxph(fit_frm,data = df_str, weights = 1/df_iptw_str$cpap_wt)
          
          #- get all coefficients
          fit_summ<-summary(fit_mort_cov)$coefficients
          fit_var<-rownames(fit_summ)
          rownames(fit_summ)<-NULL
          result[[paste0(strata$var[i],"_",strata$val[i])]]<-list(
            summ = cbind(stratum_var=strata$var[i],
                         stratum_val=strata$val[i],
                         fit_var=fit_var,
                         fit_summ),
            model = fit_mort_cov)
          #----------------------------------              
          print(paste0(mace_endpt,"...",adh,"...",boot_i,"...","stratified by:",
                       strata$var[i]," - ",strata$val[i]))
        }
        #----------------------------------
        print(paste0(mace_endpt,"...",adh,"...",boot_i,"...","stratified models collected"))
        saveRDS(result,file=file.path(path_to_dir,"coxph_strata_main_acm.rda"))
      }
    }
  }
}
