#################################################################
# Copyright (c) 2021-2025 University of Missouri                   
# Author: Xing Song, xsm7f@umsystem.edu                            
# File: cpap_exposure_mace.R
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
  gtsummary,
  tableone,
  tidycmprsk
)

source_url("https://raw.githubusercontent.com/sxinger/utils/master/sample_util.R")
source_url("https://raw.githubusercontent.com/sxinger/utils/master/analysis_util.R")
source_url("https://raw.githubusercontent.com/sxinger/utils/master/model_util.R")

#==== load data ======================
path_to_data_folder<-file.path(
  getwd(),
  "data"
)
df<-readRDS(file.path(path_to_data_folder,"cpap_exposure_final.rda")) %>%
  filter(MACE_HIST==0&MACE_time>0)

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
    rep("ACG",length(unique(df$ACG))),
    rep("AHT",length(unique(df$AHT))),
    rep("ALP",length(unique(df$ALP))),
    rep("BGR",length(unique(df$BGR)))
  )
)

mace_endpts<-c(
  'MACE',
  'HF',
  'MI',
  'STROKE',
  'REVASC'
)

#==== MACE models ====
boots_iter<-1:1
for(mace_endpt in mace_endpts){
  # mace_endpt<-'MACE' # only for testing
  #==== matched sampling ====
  boots<-5
  path_to_file<-file.path(
    path_to_data_folder,
    paste0("cpap_exposure_",tolower(mace_endpt),"_adjset_bts",boots,".rda"
    ))
  if(!file.exists(path_to_file)){
    adj_set<-matched_sample.ptdm(
      ref_dat = df %>% filter(CPAP_IND==1),
      match_dat = df %>% filter(CPAP_IND==0),
      id_col = 'PATID',
      update_ref = 'DAYS_OSA_TO_CPAP_INIT',
      update_col = paste0(mace_endpt,'_time'),
      boots = boots,
      replace = TRUE
    )
    saveRDS(adj_set,file=path_to_file)
  }else{
    adj_set<-readRDS(path_to_file)
  }
  #----------------------------------
  print("prep:sample debiasing")
  
  for(boot_i in boots_iter){
    #==== create subfolder structure
    #--create subdir
    path_to_dir<-file.path(path_to_data_folder,mace_endpt)
    if(!dir.exists(path_to_dir)) dir.create(path_to_dir)
    #--create subdir/subdir
    path_to_dir<-file.path(path_to_data_folder,mace_endpt,paste0("boot",boot_i))
    if(!dir.exists(path_to_dir)) dir.create(path_to_dir)
    
    #==== construct matching sample
    adj_col<-paste0(mace_endpt,"_time")
    df_adj<-adj_set[[boot_i]]$pos %>% 
      select(PATID) %>% 
      left_join(df,by="PATID") %>%
      bind_rows(
        adj_set[[boot_i]]$neg %>% 
          select(PATID,adj_time) %>%
          left_join(df,by="PATID") %>% 
          mutate(!!adj_col:= adj_time) %>% 
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
    
    gc()
    #===== IPTW-adjusted, Main Effect, Cause-specific
    path_to_file<-file.path(path_to_dir,"coxph_iptw_main.rda")
    if(!file.exists(path_to_file)){
      # align
      df_iptw<-df_adj %>% 
        select(PATID) %>% 
        left_join(ipw_df,by="PATID")
      
      # fit cause-specific
      cox_frm<-formula(
        paste0(
          "Surv(",mace_endpt,"_time, ",mace_endpt,"_status) ~ ",
          paste(c(var_all,"CPAP_IND"),collapse = "+"))
      )
      fit_cox<-coxph(
        cox_frm,
        data = df_adj,
        weights = df_iptw$iptw,
        model = TRUE
      )
      # ggforest(fit_cox) #quick print
      out<-list(
        fit_cox = fit_cox,
        cox_res = cox.zph(fit_cox)
      )
      saveRDS(out, file=path_to_file)
      
      #----------------------------------
      print("coxph: global model")
    }
    
    gc()
    #==== IPTW-adjusted, Main Effect, Stratified, Cause-specific
    path_to_file<-file.path(path_to_dir,"coxph_strata_main.rda")
    if(!file.exists(path_to_file)){
      result<-c()
      for(i in seq_len(nrow(strata))){
        # align 
        df_str<-df_adj %>% filter(.data[[strata$var[i]]]==strata$val[i])
        df_iptw_str<-df_str %>% select(PATID) %>% left_join(df_iptw,by="PATID")
        var_filter<-var_all[!grepl(strata$excld[i],var_all)]
        
        # fit cause-specific
        cox_frm<-formula(
          paste0(
            "Surv(",mace_endpt,"_time, ",mace_endpt,"_status) ~ ",
            paste(c(var_all,"CPAP_IND"),collapse = "+"))
        )
        fit_cox_str<-coxph(
          cox_frm,
          data = df_str, 
          weights = df_iptw_str$iptw,
          model = TRUE
        )
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
        print(paste0("coxph: stratified by:",strata$var[i],":",strata$val[i]))
      }
      
      # save results
      saveRDS(result,file=path_to_file)
    }
    
    gc()
    # #===== IPTW-adjusted, Main Effect, Fine-Gray Regression
    # path_to_file<-file.path(path_to_dir,"cmprsk_iptw_main.rda")
    # if(!file.exists(path_to_file)){
    #   # align
    #   df_iptw<-df_adj %>% 
    #     select(PATID) %>% 
    #     left_join(ipw_df,by="PATID")
    # 
    #   # fit fine-gray
    #   fgr_frm<-formula(
    #     paste0(
    #       "Surv(",mace_endpt,"_time, ",mace_endpt,"_status_sub) ~ ",
    #       paste(c(var_all,"CPAP_IND"),collapse = "+"))
    #   )
    #   # crr requires status column to be factor
    #   df_adj %<>% 
    #     mutate(!!paste0(mace_endpt,"_status_sub"):= as.factor(get(paste0(mace_endpt,"_status_sub"))))
    #   # https://stats.stackexchange.com/questions/13151/competing-risk-regression-with-crr-slow-on-large-datasets
    #   fit_fgr<-crr(
    #     fgr_frm,
    #     data = df_adj,
    #     weights = df_iptw$iptw,
    #     model = TRUE,
    #     variance = FALSE
    #   )
    #   # ggforest(fit_fgr) #quick print
    #   
    #   # save model
    #   out<-list(
    #     fit_fgr = fit_fgr,
    #     fgr_res = cox.zph(fit_fgr)
    #   )
    #   saveRDS(out, file=path_to_file)
    #   
    #   #----------------------------------
    #   print("fine-gray: global model")
    # }
    # 
    # gc()
    # #==== IPTW-adjusted, Main Effect, Stratified, Fine-Gray Regression
    # path_to_file<-file.path(path_to_dir,"cmprsk_strata_main.rda")
    # if(!file.exists(path_to_file)){
    #   result<-c()
    #   for(i in seq_len(nrow(strata))){
    #     # align 
    #     df_str<-df_adj %>% filter(.data[[strata$var[i]]]==strata$val[i])
    #     df_iptw_str<-df_str %>% select(PATID) %>% left_join(df_iptw,by="PATID")
    #     var_filter<-var_all[!grepl(strata$excld[i],var_all)]
    #     
    #     # fit fine-gray
    #     fgr_frm<-formula(
    #       paste0(
    #         "Surv(",mace_endpt,"_time, ",mace_endpt,"_status_sub) ~ ",
    #         paste(c(var_all,"CPAP_IND"),collapse = "+"))
    #     )
    #     fit_fgr_str<-crr(
    #       fgr_frm,
    #       data = df_str,
    #       weights = df_iptw_str$iptw,
    #       model = TRUE,
    #       variance = FALSE
    #     )
    #     fit_summ<-summary(fit_fgr_str)$coefficients
    #     fit_var<-rownames(fit_summ)
    #     rownames(fit_summ)<-NULL
    #     result<-rbind(
    #       result,
    #       cbind(stratum_var=strata$var[i],
    #             stratum_val=strata$val[i],
    #             fit_var=fit_var,
    #             fit_summ)
    #     )
    #     
    #     #----------------------------------
    #     print(paste0("fine-gray: stratified by:",strata$var[i],":",strata$val[i]))
    #   }
    #   
    #   # save results
    #   saveRDS(result,file=path_to_file)
    # }
  }
}
