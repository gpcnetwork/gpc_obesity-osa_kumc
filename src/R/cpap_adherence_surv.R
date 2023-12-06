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
  glmnetcr,
  inflection
)

source_url("https://raw.githubusercontent.com/sxinger/utils/master/sample_util.R")
source_url("https://raw.githubusercontent.com/sxinger/utils/master/analysis_util.R")
source_url("https://raw.githubusercontent.com/sxinger/utils/master/model_util.R")

#==== load data
path_to_data_folder<-file.path(
  getwd(),
  "data"
)
df<-readRDS(file.path(path_to_data_folder,"cpap_adherence_final.rda"))

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

# repeat experiments over bootstrapped samples
adh_metric<-c(
   "CPAP_YR1"
  ,"adherence_yr1_med"
  ,"adherence_yr1_tt"
  ,"adherence_yr1_qt"
  ,"adherence_yr1_inc4"
  ,"adherence_yr1_inc8"
  ,"adherence_yr1_ww"
  ,"adherence_yr1_mttree"
)
boots_iter<-1:1
for(adh in adh_metric){
  for(boot_i in boots_iter){
    # adh<-adh_metric[1]
    boot_i<-1
    #==== create subfolder structure
    #--create subdir
    path_to_dir<-file.path(path_to_data_folder,"ADHRN_ACM")
    if(!dir.exists(path_to_dir)) dir.create(path_to_dir)
    #--create subdir/subdir
    path_to_dir<-file.path(path_to_dir,tolower(adh))
    if(!dir.exists(path_to_dir)) dir.create(path_to_dir)
    #--create subdir/subdir
    path_to_dir<-file.path(path_to_dir,paste0("boot",boot_i))
    if(!dir.exists(path_to_dir)) dir.create(path_to_dir)
    
    #==== propensity score
    path_to_file<-file.path(path_to_dir,"ps_papadh_elasticnet.rda")
    if(!file.exists(path_to_file)){
      df_ps<-df[,c(adh,var_all)] %>% 
        filter(!is.na(get(adh)))
      x<-data.matrix(df_ps[,-1])
      y<-unlist(df_ps[,1]) 
      if(adh == 'CPAP_YR1'){
        fit_ps<-cv.glmnet(x=x,y=y,family = "poisson")
        # plot(fit_ps)
        ps_cpap<-predict(fit_ps, newx = x, s = "lambda.min", type="response")
        wt_long<-data.frame(
          PATID=df$PATID,
          lambda=ps_cpap[,1],
          tgt = y
        ) %>%
          mutate(
            explambda = exp(-lambda),
            lambdak=(lambda)^tgt,
            kfac=factorial(tgt),
            wt_den = explambda*lambdak/kfac
          ) %>%
          group_by(tgt) %>%
          mutate(wt_num = length(unique(PATID))/length(unique(df$PATID))) %>%
          ungroup
      }else{
        fit_ps<-cv.glmnet(x=x,y=y,family = "multinomial", type.multinomial = "grouped")
        # plot(fit_ps)
        ps_cpap<-predict(fit_ps, newx = x, s = "lambda.min", type="response")
        wt_long<-data.frame(
          PATID=df$PATID,
          lambda=ps_cpap[,,1],
          tgt = y
        ) %>%
          gather(pred,ps,-tgt,-PATID) %>%
          mutate(pred = gsub("lambda\\.","",pred)) %>%
          filter(pred==tgt) %>%
          mutate(wt_den = as.numeric(ps)) %>%
          group_by(tgt) %>%
          mutate(wt_num = length(unique(PATID))/length(unique(df$PATID))) %>%
          ungroup
      }
      ipw_df<-ipw.naive(
        wt_long = wt_long %>% 
          mutate(time=1),
        id_col = 'PATID',
        ot_cols = 'tgt',
        truncate = TRUE,
        truncate_lower = 0.1,
        truncate_upper = 0.9
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
    #--build IPW model
    path_to_file<-file.path(path_to_dir,"coxph_iptw_main.rda")
    if(!file.exists(path_to_file)){
      # align
      df_iptw<-df %>% select(PATID) %>% 
        left_join(
          ipw_df %>% select(PATID,iptw),
          by="PATID")
      
      # fit
      if(adh == 'CPAP_YR1'){
        fit_frm<-formula(
          paste0(
            "Surv(DEATH_time, DEATH_status) ~ ",
            paste(var_all,collapse = "+"),
            "+ pspline(",adh,", df = 4)"
          )
        )
        fit_cox<-coxph(
          fit_frm,
          data = df, 
          weights = df_iptw$iptw,
          model = TRUE
        )
        tm<-names(fit_cox$pterms)
        tm_idx<-which(grepl(adh,tm), arr.ind = FALSE)
        tp<-termplot(
          fit_cox, 
          data = df,
          term=tm_idx, 
          se=TRUE, 
          col.term=1, 
          col.se=1,plot=F
        )
        saveRDS(tp,file=file.path(path_to_dir,"pspline_termplot.rda"))
      }else{
        # model formulation
        fit_frm<-formula(
          paste0(
            "Surv(DEATH_time, DEATH_status) ~ ",
            paste(c(var_all,adh),collapse = "+")
          )
        )
        fit_cox<-coxph(
          fit_frm,
          data = df, 
          weights = df_iptw$iptw,
          model=TRUE
        )
      }
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
      print(paste0(adh,":coxph: global model"))
    }

    #==== IPTW-adjusted, Main Effect, Stratified
    result<-c()
    path_to_file<-file.path(path_to_dir,"coxph_strata_main.rda")
    if(!file.exists(path_to_file)){
      for(i in seq_len(nrow(strata))){
        # align
        df_str<-df %>% filter(.data[[strata$var[i]]]==strata$val[i])
        df_iptw_str<-df_str %>% select(PATID) %>% left_join(df_iptw,by="PATID")
        
        #fit
        var_filter<-var_all[!grepl(strata$excld[i],var_all)]
        if(adh == 'CPAP_YR1'){
          fit_frm<-formula(
            paste0(
              "Surv(DEATH_time, DEATH_status) ~ ",
              paste(var_filter,collapse = "+"),
              "+ pspline(",adh,", df = 4)"
              )
            )
        }else{
          fit_frm<-formula(
            paste0(
              "Surv(DEATH_time, DEATH_status) ~ ",
              paste(c(var_filter,adh),collapse = "+")
              )
            )
        }
        fit_cox_str<-coxph(
          fit_frm,
          data = df_str, 
          weights = df_iptw_str$iptw,
          model = TRUE
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
        print(paste0(adh,":stratified by:",strata$var[i],":",strata$val[i]))
      }
      # save results
      saveRDS(result,file=path_to_file)
    }
  }
  
  #------------------------------------------------
  print(paste0("finish adherence measure:",adh))
}
