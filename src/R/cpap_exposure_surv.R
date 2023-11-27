#################################################################
# Copyright (c) 2021-2025 University of Missouri                   
# Author: Xing Song, xsm7f@umsystem.edu                            
# File: cpap_exposure_surv.R
# Description: model how cpap exposure (binary) affect survival 
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

#==== load data
path_to_data_folder<-file.path(
  getwd(),
  "data"
)
df<-readRDS(file.path(
  path_to_data_folder,
  "cpap_exposure_final.rda"
))

#==== matched sampling
boots<-5
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
var_demo<-c(
  "AGEGRP_agegrp2", "AGEGRP_agegrp3","AGEGRP_agegrp4", #ref=AGEGRP_agegrp1
  "SEXF",
  "RACE_AA","RACE_AI","RACE_Asian","RACE_Other","RACE_Unknown", #ref=RACE_WH
  "LIS_DUAL_IND"  
)

var_comorb<-c(
   'CCI_SCORE'
  ,'OBES'
  ,'T2DM'
  ,'CKD'
  ,'COPD'
  ,'AFIB'
  ,'HTN'
  ,'ND'
  ,'HYSM'
  ,'INSM'
  ,paste0(c('MI','HF','STROKE','REVASC',"MACE"),'_HIST')
)

var_med<-c(
   'ACG'
  ,'AHT'
  ,'ALP'
  ,'BGR'
)

var_all<-c(
  var_demo,
  var_comorb,
  var_med
)

# model formulation
fit_main<-formula(
  paste0(
   "Surv(DEATH_time, DEATH_status) ~ ",
    paste(c(var_all,"CPAP_IND"),collapse = "+")
  )
)

# define strata
strata<-data.frame(
  val=c(
    unique(df$AGEGRP),
    unique(df$SEX),
    unique(as.character(df$RACE_LABEL)),
    unique(df$INSM),
    unique(df$HYSM),
    unique(df$LIS_DUAL_IND),
    unique(df$CCI_CLASS),
    unique(df$OBES),
    unique(df$T2DM),
    unique(df$CKD),
    unique(df$AFIB),
    unique(df$HTN),
    unique(df$HTN),
    unique(df$ND),
    unique(df$MACE_HIST),
    unique(df$ACG),
    unique(df$AHT),
    unique(df$ALP),
    unique(df$BGR) 
  ),
  var=c(
    rep("AGEGRP",length(unique(df$AGEGRP))),
    rep("SEX",length(unique(df$SEX))),
    rep("RACE_LABEL",length(unique(df$RACE_LABEL))),
    rep("INSM",length(unique(df$INSM))),
    rep("HYSM",length(unique(df$HYSM))),
    rep("LIS_DUAL_IND",length(unique(df$LIS_DUAL_IND))),
    rep("CCI_CLASS",length(unique(df$CCI_CLASS))),
    rep("OBES",length(unique(df$OBES))),
    rep("T2DM",length(unique(df$T2DM))),
    rep("CKD",length(unique(df$CKD))),
    rep("COPD",length(unique(df$COPD))),
    rep("AFIB",length(unique(df$AFIB))),
    rep("HTN",length(unique(df$HTN))),
    rep("ND",length(unique(df$ND))),
    rep("MACE_HIST",length(unique(df$MACE_HIST))),
    rep("ACG",length(unique(df$ACG))),
    rep("AHT",length(unique(df$AHT))),
    rep("ALP",length(unique(df$ALP))),
    rep("BGR",length(unique(df$BGR)))
  ),
  excld=c(
    rep("AGE",length(unique(df$AGEGRP))),
    rep("SEXF",length(unique(df$SEX))),
    rep("RACE",length(unique(df$RACE_LABEL))),
    rep("INSM",length(unique(df$INSM))),
    rep("HYSM",length(unique(df$HYSM))),
    rep("LIS_DUAL_IND",length(unique(df$LIS_DUAL_IND))),
    rep("CCI",length(unique(df$CCI_CLASS))),
    rep("OBES",length(unique(df$OBES))),
    rep("T2DM",length(unique(df$T2DM))),
    rep("CKD",length(unique(df$CKD))),
    rep("COPD",length(unique(df$COPD))),
    rep("AFIB",length(unique(df$AFIB))),
    rep("HTN",length(unique(df$HTN))),
    rep("ND",length(unique(df$ND))),
    rep("MACE",length(unique(df$MACE_HIST))),
    rep("ACG",length(unique(df$ACG))),
    rep("AHT",length(unique(df$AHT))),
    rep("ALP",length(unique(df$ALP))),
    rep("BGR",length(unique(df$BGR)))
  )
)

#==== main effect models ====
boots_iter<-1:2
for(boot_i in boots_iter){
  #==== create subfolder structure
  #--create subdir
  path_to_dir<-file.path(path_to_data_folder,"surv")
  if(!dir.exists(path_to_dir)) dir.create(path_to_dir)
  #--create subdir/subdir
  path_to_dir<-file.path(path_to_data_folder,"surv",paste0("boot",boot_i))
  if(!dir.exists(path_to_dir)) dir.create(path_to_dir)
  
  #==== construct matching sample
  df_adj<-adj_set[[boot_i]]$pos %>% 
    select(PATID) %>% 
    left_join(df,by="PATID") %>%
    bind_rows(
      adj_set[[boot_i]]$neg %>% 
        select(PATID,adj_time) %>%
        left_join(df,by="PATID") %>% 
        mutate(DEATH_time = adj_time) %>% 
        select(-adj_time)
    )
  
  #==== propensity score
  path_to_file<-file.path(path_to_dir,"ps_pap_elasticnet.rda")
  if(!file.exists(path_to_file)){
    x<-data.matrix(df_adj[,var_all])
    y<-df_adj$CPAP_IND
    fit_ps<-cv.glmnet(x=x,y=y,family="binomial")
    # plot(fit_ps)
    ps_cpap<-predict(
      fit_ps, 
      newx = x, 
      s = "lambda.min", 
      type="response"
    )
    wt_long<-data.frame(
      PATID=df_adj$PATID,
      time=1,
      wt_num=1,
      p=ps_cpap[,1],
      tgt=y
    ) %>%
      mutate(
        wt_den = tgt*p + (1-tgt)*(1-p)
      )
    ipw_df<-ipw.naive(
      wt_long = wt_long,
      id_col = 'PATID',
      ot_cols = 'tgt',
      truncate = TRUE,
      truncate_lower = 0.0001,
      truncate_upper = 0.99
    ) %>% ungroup
    
    # save ps model
    ps_out<-list(
      fit_ps = fit_ps,
      ps_score = ps_cpap,
      ipw = ipw_df
    )
    saveRDS(ps_out, file=path_to_file)
  }else{
    ipw_df<-readRDS(path_to_file)$ipw
  }

  #===== IPTW-adjusted, Main Effect 
  path_to_file<-file.path(path_to_dir,"coxph_iptw_main.rda")
  if(!file.exists(path_to_file)){
    # align
    df_iptw<-df_adj %>% select(PATID) %>% 
      left_join(
        ipw_df %>% select(PATID,iptw),
        by="PATID")
    
    # fit
    fit_cox<-coxph(
      fit_main,
      data = df_adj, 
      weights = df_iptw$iptw,
      model=TRUE
    )
    # ggforest(fit_cox) #quick print
    
    # proportionality check
    sch_res<-cox.zph(fit_cox)
    # ggcoxzph(sch_res) #quick print - caution! may kill the R session
    
    # save model
    out<-list(
      fit_cox = fit_cox,
      sch_res = sch_res
    )
    saveRDS(out, file=path_to_file)
    
    #----------------------------------
    print("coxph: global model")
  }
  
  #==== IPTW-adjusted, Main Effect, Stratified
  result<-c()
  path_to_file<-file.path(path_to_dir,"coxph_strata_main.csv")
  if(!file.exists(path_to_file)){
    for(i in seq_len(nrow(strata))){
      # align
      df_str<-df_adj %>% filter(.data[[strata$var[i]]]==strata$val[i])
      df_iptw_str<-df_str %>% select(PATID) %>% left_join(df_iptw,by="PATID")
      
      #fit
      var_filter<-var_all[!grepl(strata$excld[i],var_all)]
      fit_frm<-formula(
        paste0(
          "Surv(DEATH_time, DEATH_status) ~ ",
          paste(c(var_filter,"CPAP_IND"),collapse = "+"))
      )
      fit_cox_str<-coxph(
        fit_frm,
        data = df_str, 
        weights = df_iptw_str$iptw,
        model=TRUE
      )

      # get all coefficients
      fit_summ<-summary(fit_cox_str)$coefficients
      fit_var<-rownames(fit_summ)
      rownames(fit_summ)<-NULL
      result<-rbind(
        result,
        cbind(stratum_var=strata$var[i],
              stratum_val=strata$val[i],
              fit_var=fit_var,
              fit_summ)
      )
      #----------------------------------
      print(paste0("stratified by:",strata$var[i],":",strata$val[i]))
    }
    # save results
    saveRDS(result,file=path_to_file)
  }
}


#==== interaction effect models ====
# var_interact<-c(
#   "CPAP_x_INSM",
#   "CPAP_x_HTSM",
#   "CPAP_x_COPD",
#   "CPAP_x_T2DM",
#   "CPAP_x_CKD",
#   "CPAP_x_HTN",
#   "CPAP_x_OBES",
#   "CPAP_x_AFIB",
#   "CPAP_x_ND",
#   "CPAP_x_ACG",
#   "CPAP_x_AHT",
#   "CPAP_x_ALP",
#   "CPAP_x_BGR",
#   "CPAP_x_MACE",
#   "CPAP_x_SES"
# )
# 
# fit_intx<-formula(
#   paste0(
#     "Surv(DEATH_time, DEATH_status) ~ ",
#     paste(
#       c(
#         var_all,
#         "CPAP_IND",
#         var_interact
#         ),
#       collapse = "+")
#   )
# )
# 
# for(boot_i in boots_iter){
#   #==== create subfolder structure
#   #--create subdir
#   path_to_dir<-file.path(path_to_data_folder,"surv")
#   if(!dir.exists(path_to_dir)) dir.create(path_to_dir)
#   #--create subdir/subdir
#   path_to_dir<-file.path(path_to_data_folder,"surv",paste0("boot",boot_i))
#   if(!dir.exists(path_to_dir)) dir.create(path_to_dir)
#   
#   #==== construct matching sample
#   df_adj<-adj_set[[boot_i]]$pos %>% 
#     select(PATID) %>% 
#     left_join(df,by="PATID") %>%
#     bind_rows(
#       adj_set[[boot_i]]$neg %>% 
#         select(PATID,adj_time) %>%
#         left_join(df,by="PATID") %>% 
#         mutate(DEATH_time = adj_time) %>% 
#         select(-adj_time)
#     )
#   
#   #==== propensity score
#   #-- same as the ps model for main-effect models
#   path_to_file<-file.path(path_to_dir,"ps_pap_elasticnet.rda")
#   if(!file.exists(path_to_file)){
#     x<-data.matrix(df_adj[,c(var_demo,var_comorb,var_med)])
#     y<-df_adj$CPAP_IND
#     fit_ps<-cv.glmnet(x=x,y=y,family="binomial")
#     # plot(fit_ps)
#     ps_cpap<-predict(
#       fit_ps, 
#       newx = x, 
#       s = "lambda.min", 
#       type="response"
#     )
#     ipw_df<-ipw.naive(
#       wt_long = data.frame(
#         PATID=df_adj$PATID,
#         time=1,
#         wt_den=ps_cpap[,1],
#         wt_num=1,
#         tgt=y
#       ),
#       id_col = 'PATID',
#       ot_cols = 'tgt',
#       truncate = TRUE,
#       truncate_lower = 0.0001,
#       truncate_upper = 0.99
#     ) %>% ungroup
#     
#     # save ps model
#     ps_out<-list(
#       fit_ps = fit_ps,
#       ps_score = ps_cpap,
#       ipw = ipw_df
#     )
#     saveRDS(ps_out, file=path_to_file)
#   }else{
#     ipw_df<-readRDS(path_to_file)$ipw
#   }
#   
#   #===== IPTW-adjusted, Interactive Effect
#   path_to_file<-file.path(path_to_dir,"coxph_iptw_intx.rda")
#   if(!file.exists(path_to_file)){
#     # align
#     df_iptw<-df_adj %>% select(PATID) %>%
#       left_join(ipw_df %>% select(PATID,iptw),
#                 by="PATID")
#     # fit
#     fit_cox<-coxph(
#       fit_intx,
#       data = df_adj,
#       weights = df_iptw$iptw,
#       model=TRUE
#     )
#     # result quick print
#     # ggforest(fit_cox)
#     # save model - takes about 5min to train a msm model
#     saveRDS(fit_cox, file=path_to_file)
#     gc()
#   }
# }
