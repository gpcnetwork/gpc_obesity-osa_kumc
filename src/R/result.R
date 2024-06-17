rm(list=ls()); gc()

pacman::p_load(
  survival,
  survminer,
  tidyverse,
  magrittr,
  kableExtra,
  gridExtra,
  devtools,
  grid,
  forestploter,
  ggrepel,
  ggpubr,
  svglite
  # adjustedCurves
)

source_url("https://raw.githubusercontent.com/sxinger/utils/master/plot_util.R")

path_to_datadir<-file.path(getwd(),"data")
path_to_outdir<-file.path(getwd(),"res")

#==== KM plots ====
#------ ACM,KM summary
pap_use<-readRDS(file.path(path_to_datadir,"cpap_exposure_final.rda"))
summary(survfit(Surv(DEATH_time,DEATH_status) ~ 1, data = pap_use),times = 365*c(1:5))

#--- ACM,unadj
sfit_obj<-survfit(Surv(DEATH_time,DEATH_status) ~ CPAP_IND, data = pap_use)
risk_tbl<-summary(sfit_obj,times = 365*c(1:5))
km_mort_unadj<-ggsurvplot(
  fit = sfit_obj,
  # pval = TRUE,
  conf.int = TRUE,
  legend.labs=c("wo/ PAP", "w/ PAP"),
  # risk.table = TRUE,
  linetype = "strata",
  break.x.by = 365,
  xlim = c(0, 1825),
  xlab = "Days since Index Date", 
  ylab = "Unadjusted Survival Probability"
)

km_mort_unadj2<-km_mort_unadj$plot +
  geom_vline(xintercept=365*c(1:5),linetype=2)+
  geom_label_repel(
    data=data.frame(
      x=risk_tbl$time,
      y=risk_tbl$surv,
      label=round(risk_tbl$surv,3),
      label_int=paste0(round(risk_tbl$surv,3),"[",round(risk_tbl$lower,3),",",round(risk_tbl$upper,3),"]")),
    aes(x=x,y=y,label=label)) +
  geom_text(aes(x=150, y=0.2, label = "p < 0.0001", fontface=0))

saveRDS(
  km_mort_unadj2,
  file=file.path(path_to_datadir,"output","expos_ACM_unadjKM.rda")
)
rm(km_mort_unadj,km_mort_unadj2,risk_tbl,sfit_obj); gc()

#--median follow-up
survfit(Surv(DEATH_time,1-DEATH_status) ~ 1, data = pap_use)

#--- ACM, adj
fitcox<-readRDS(file.path(path_to_datadir,".archive","ACM","boot1","coxph_iptw_main.rda"))$fit_cox
pval<-summary(fitcox)$coefficients["CPAP_IND",6]
pred_df<-pap_use %>% select(attr(fitcox$means,"names")) %>% select(-CPAP_IND) %>% unique %>% sample_frac(0.1)

sfit0<-summary(survfit(fitcox,newdata=pred_df %>% mutate(CPAP_IND=0),conf.int = T)); gc() # predicted risks
time_scale<-sfit0$time
adjkm_df<-data.frame(
  time = time_scale,
  surv = apply(sfit0$surv,1,mean),
  std = apply(sfit0$std.err,1,mean),
  # lower = apply(sfit0$surv,1,function(x) quantile(x,0.025)),
  # upper = apply(sfit0$surv,1,function(x) quantile(x,0.975)),
  CPAP_IND = rep(0,length(time_scale))
)
rm(sfit0); gc()

sfit1<-summary(survfit(fitcox,newdata=pred_df %>% mutate(CPAP_IND=1),conf.int = T)); gc() # predicted risks
adjkm_df %<>%
  bind_rows(
    data.frame(
      time = sfit1$time,
      surv = apply(sfit1$surv,1,mean),
      std = apply(sfit1$std.err,1,mean),
      # lower = apply(sfit1$surv,1,function(x) quantile(x,0.025)),
      # upper = apply(sfit1$surv,1,function(x) quantile(x,0.975)),
      CPAP_IND = rep(1,length(time_scale))
    )
  ) %>%
  filter(time > 0) %>%
  mutate(
    lower = surv - 1.96*std,
    upper = surv + 1.96*std
  ) %>%
  mutate(
    CPAP_IND_lbl = case_when(
      CPAP_IND==1 ~ 'Evidence of PAP initiation',
      TRUE ~ 'No Evidence of PAP initiation'
    )
  )
rm(sfit1); gc()

risk_tbl2<-adjkm_df %>%  filter(time %in% c(365,730,1095,1460,1825))
km_mort_adj<-ggplot(adjkm_df,aes(x=time,y=surv)) +
  geom_step(aes(group=CPAP_IND,color = CPAP_IND_lbl),linewidth=2.5) +
  geom_ribbon(aes(fill = CPAP_IND_lbl, ymin = lower, ymax = upper),alpha=0.5) + 
  geom_vline(xintercept=365*c(1:5),linetype=2) +
  # geom_label_repel(
  #   data=data.frame(
  #     x=risk_tbl2$time,
  #     y=risk_tbl2$surv,
  #     label=round(risk_tbl2$surv,3),
  #     label_int=paste0(round(risk_tbl2$surv,3),"[",round(risk_tbl2$lower,3),",",round(risk_tbl2$upper,3),"]")),
  #   aes(x=x,y=y,label=label)
  # ) + 
  geom_text(aes(x=170, y=0.2, label = ifelse(pval<0.0001, "p < 0.0001",pval)), fontface=0) +
  scale_x_continuous(
    breaks = c(0,365,730,1095,1460,1825),
    labels = c(0,365,730,1095,1460,1825),
    limits = c(0, 1825)
  ) +
  ylim(0,1) + 
  ylab("Adjusted Survival Probability") + 
  xlab("Days since Index Date") +
  theme_survminer()+
  scale_color_discrete(name="")+
  scale_fill_discrete(name="")

saveRDS(
  km_mort_adj,
  file=file.path(path_to_datadir,"output","expos_ACM_adjKM.rda")
)

ggsave(
  file=file.path(path_to_outdir,"expos_ACM_adjKM.svg"), 
  plot=km_mort_adj, 
  width=8, height=4
)

rm(km_mort_adj,risk_tbl2,fitcox,pred_df,adjkm_df); gc()

#--- MACE, KM summary
pap_use2<-pap_use %>% filter(MACE_HIST==0&MACE_time>0)
rm(pap_use);gc()
summary(survfit(Surv(MACE_time,MACE_status) ~ 1, data = pap_use2),times = 365*c(1:5))

#--- MACE, unadj
sfit_obj2<-survfit(Surv(MACE_time,MACE_status) ~ CPAP_IND, data = pap_use2)
risk_tbl<-summary(sfit_obj2, data = pap_use2,times = 365*c(1:5))
km_mace_unadj<-ggsurvplot(
  fit = sfit_obj2,
  # pval = TRUE,
  conf.int = TRUE,
  legend.labs=c("wo/ PAP", "w/ PAP"),
  # risk.table = TRUE,
  linetype = "strata",
  break.x.by = 365,
  xlim = c(0, 1825),
  xlab = "Days since Index Date", 
  ylab = "Unadjusted MACE-free Probability")

km_mace_unadj2<-km_mace_unadj$plot +
  geom_vline(xintercept=365*c(1:5),linetype=2)+
  # geom_label_repel(
  #   data=data.frame(
  #     x=risk_tbl$time,
  #     y=risk_tbl$surv,
  #     label=round(risk_tbl$surv,3),
  #     label_int=paste0(round(risk_tbl$surv,3),"[",round(risk_tbl$lower,3),",",round(risk_tbl$upper,3),"]")
  #   ),
  #   aes(x=x,y=y,label=label)
  # ) +
  geom_text(aes(x=150, y=0.2, label = "p < 0.0001", fontface=0))

saveRDS(
  km_mace_unadj2,
  file=file.path(path_to_datadir,"output","expos_MACE_unadjKM.rda")
)
rm(km_mace_unadj,km_mace_unadj2,risk_tbl,sfit_obj2); gc()

#--- MACE, adj 
fitcox<-readRDS(file.path(path_to_datadir,"MACE","boot1","coxph_iptw_main.rda"))$fit_cox
pred_df<-pap_use2 %>% select(attr(fitcox$means,"names")) %>% select(-CPAP_IND) %>% unique %>% sample_frac(0.1)
pval<-summary(fitcox)$coefficients["CPAP_IND",6]

sfit0<-summary(survfit(fitcox,newdata=pred_df %>% mutate(CPAP_IND=0),conf.int = T)); gc() # predicted risks
time_scale<-sfit0$time
adjkm_df<-data.frame(
  time = time_scale,
  surv = apply(sfit0$surv,1,mean),
  std = apply(sfit0$std.err,1,mean),
  # lower = apply(sfit0$surv,1,function(x) quantile(x,0.025)),
  # upper = apply(sfit0$surv,1,function(x) quantile(x,0.975)),
  CPAP_IND = rep(0,length(time_scale))
) 
rm(sfit0); gc()

sfit1<-summary(survfit(fitcox,newdata=pred_df %>% mutate(CPAP_IND=1),conf.int = T)); gc() # predicted risks
adjkm_df %<>%
  bind_rows(
    data.frame(
      time = sfit1$time,
      surv = apply(sfit1$surv,1,mean),
      std = apply(sfit1$std.err,1,mean),
      # lower = apply(sfit1$surv,1,function(x) quantile(x,0.025)),
      # upper = apply(sfit1$surv,1,function(x) quantile(x,0.975)),
      CPAP_IND = rep(1,length(time_scale))
    )
  ) %>%
  filter(time > 0) %>%
  mutate(
    lower = surv - 1.96*std,
    upper = surv + 1.96*std
  ) %>%
  mutate(
    CPAP_IND_lbl = case_when(
      CPAP_IND==1 ~ 'Evidence of PAP initiation',
      TRUE ~ 'No Evidence of PAP initiation'
    )
  )
rm(sfit1); gc()

risk_tbl2<-adjkm_df %>% filter(time %in% c(365,730,1095,1460,1825))
km_mace_adj<-ggplot(adjkm_df,aes(x=time,y=surv)) +
  geom_step(aes(group=CPAP_IND,color = CPAP_IND_lbl),linewidth=2.5) +
  geom_ribbon(aes(fill = CPAP_IND_lbl, ymin = lower, ymax = upper),alpha=0.5) + 
  geom_vline(xintercept=365*c(1:5),linetype=2) +
  # geom_label_repel(
  #   data=data.frame(
  #     x=risk_tbl2$time,
  #     y=risk_tbl2$surv,
  #     label=round(risk_tbl2$surv,3),
  #     label_int=paste0(round(risk_tbl2$surv,3),"[",round(risk_tbl2$lower,3),",",round(risk_tbl2$upper,3),"]")),
  #   aes(x=x,y=y,label=label)) + 
  geom_text(aes(x=170, y=0.2, label = ifelse(pval<0.0001, "p < 0.0001",pval)), fontface=0) +
  scale_x_continuous(
    breaks = c(0, 365,730,1095,1460,1825),
    labels = c(0, 365,730,1095,1460,1825),
    limits = c(0, 1825)
  ) +
  ylim(0,1) + 
  ylab("Adjusted MACE-free Probability") + 
  xlab("Days since Index Date") +
  theme_survminer()+
  scale_color_discrete(name="")+
  scale_fill_discrete(name="")

saveRDS(
  km_mace_adj,
  file=file.path(path_to_datadir,"output","expos_MACE_adjKM.rda")
)

ggsave(
  file=file.path(path_to_outdir,"expos_MACE_adjKM.svg"), 
  plot=km_mace_adj, 
  width=8, height=4
)

rm(km_mace_adj,risk_tbl2,fitcox,pred_df,adjkm_df); gc()

#--- put the plots together
km_adj<-ggarrange(
  readRDS(file.path(path_to_datadir,"output","expos_ACM_adjKM.rda")),
  readRDS(file.path(path_to_datadir,"output","expos_MACE_adjKM.rda")),
  ncol = 2, common.legend = TRUE
)

svglite(file.path(path_to_outdir,"expos_adjKM2.svg"),width=10, height=4)
plot(km_adj)
dev.off()

ggsave(
  file=file.path(path_to_outdir,"expos_adjKM.svg"), 
  plot=km_adj, 
  width=10, height=4
)

#==== forestplots =====
endpts<-c(
   "ACM"
  ,"MACE"
  ,"MI"
  ,"HF"
  ,"STROKE"
  ,"REVASC"
)

nboots<-1

dat<-c()
for(endpt_i in endpts){
  # endpt_i<-endpts[2]
  # boot_i<-1
  adh_df<-c()
  for(boot_i in 1:nboots){
    # whole-data model
    #--main-effect
    boot<-readRDS(
      file.path(path_to_datadir,endpt_i,paste0("boot",boot_i),
                "coxph_iptw_main.rda"))[["fit_cox"]]
    adh_df %<>%
      bind_rows(
        data.frame(summary(boot)$coefficients) %>%
          cbind(data.frame(confint(boot))) %>%
          rownames_to_column(var = "vari") %>%
          filter(grepl("CPAP_IND",vari)) %>%
          gather(summ_var,summ_val,-vari) %>%
          dplyr::filter(summ_var %in% c("coef","X2.5..","X97.5..","Pr...z..")) %>%
          mutate(model = "main",
                 stratum_var = "None",
                 stratum_val = " ")
      )
    
    # stratified model
    #--main-effect
    boot<-readRDS(
      file.path(path_to_datadir,endpt_i,paste0("boot",boot_i),
                "coxph_strata_main.rda")) %>% as.data.frame()
    adh_df %<>%
      bind_rows(
        boot %>%
          mutate(
            vari = fit_var,
            coef = as.numeric(coef),
            `se(coef)` = as.numeric(`se(coef)`),
            `Pr...z..` = `Pr(>|z|)`
          ) %>%
          filter(grepl("CPAP_IND",vari)) %>%
          mutate(
            `X2.5..`= coef - 1.96*`se(coef)`,
            `X97.5..`= coef + 1.96*`se(coef)`
          ) %>%
          gather(summ_var,summ_val,-vari,-stratum_var,-stratum_val) %>%
          dplyr::filter(summ_var %in% c("coef","X2.5..","X97.5..","Pr...z..")) %>%
          mutate(summ_val = as.numeric(summ_val),
                 model = "main")
      )
  }
  dat %<>%
    dplyr::bind_rows(
      adh_df %>%
        dplyr::group_by(stratum_var,stratum_val,model,vari,summ_var) %>%
        dplyr::summarize(summ_val_mu=mean(summ_val),
                         summ_val_sd=sd(summ_val),
                         .groups = "drop") %>%
        dplyr::mutate(endpt = endpt_i)
    )
}

dat %<>%
  pivot_wider(
    id_cols = c("endpt","stratum_var","stratum_val","model","vari"),
    names_from = "summ_var",
    values_from = "summ_val_mu"
  )
colnames(dat)<-c("endpt","stratum_var","stratum_val","model","vari","pval","ci_lb","ci_ub","coef")

dat %<>%
  # fit forestplot default keys
  dplyr::mutate(
    pvallabel = case_when(
      pval < 0.001 ~ "<0.001",
      TRUE ~ as.character(round(pval,3))
    ),
    mean = exp(coef),
    lower = exp(ci_lb),
    upper = exp(ci_ub)
  )

dat_all<-dat %>% filter(stratum_var=="None") %>%
  mutate(stratum_var=" ")
fplt1<-forestplot.HR(
  df = dat_all %>% filter(endpt %in% c("ACM","MACE")),
  x_idx1="stratum_var", # 1st layer index
  x_idx2="stratum_val", # 2nd layer index
  y_idx="endpt", # 1st layer y index
  est="mean", # estimates
  lower="lower", # 95% CI lower bound
  upper="upper", # 95% CI upper bound
  pval="pval", # p value
  plt_par = list(
    ci_column = c(2,4),
    ref_line = c(0.55,0.95),
    xlim = list(c(0.45, 0.55),c(0.85,0.95)),
    ticks_at = list(c(0.45,0.5,0.55),c(0.85,0.9,0.95))
  ), 
  ny = 2, # number of y groups (must be the same as groups of y_idx)
  idx_display = "Full\nModel"
)
# save figure
ggsave(
  file.path(path_to_outdir,"pap_expo_on_surv_mace_full.svg"),
  plot = fplt1,
  dpi = 100,
  width = 8, 
  height = 4, 
  units = "in",
  device = "svg"
)

dat %<>%
  filter(stratum_var!="None") %>% 
  mutate(
    stratum_var_lbl = recode(
      stratum_var,
      "AGEGRP" = "01.Age",
      "SEX" = "02.Sex",
      "RACE_LABEL" = "03.Race",
      "LIS_DUAL_IND" = "04.Low-income-subsidy/Dual Eligibility",
      "HYSM" = "05.Hypersomnia",
      "INSM" = "06.Insomina",
      "MACE_HIST" = "07.MACE History",
      "OBES" = "08.Obesity",
      "COPD" = "09.COPD",
      "HTN" = "10.Hypertension",
      "T2DM" = "11.T2DM", 
      "ND" = "12.Anxiety Disorder",
      "AFIB" = "13.Atrial Fibrilation",
      "CKD" = "14.Chronic Kidney Disease",
      "CCI_CLASS" = "15.Charlson Comorbidity Index",
      "ACG" = "16.Anti-coagulant",
      "AHT" = "17.Anti-hypertensive",
      "ALP" = "18.Anti-lipemic",
      "BGR" = "19.Blood glucose regulator",
      
    ),
    stratum_val_lbl = recode(
      stratum_val,
      "agegrp1" = "65-69 years",
      "agegrp2" = "70-74 years",
      "agegrp3" = "75-79 years",
      "agegrp4" = "80+ years",
      "0" = "No",
      "1" = "Yes",
      " " = " ",
      "AA" = "Black",
      "AI" = "Native American",
      "Asian" = "Asian",
      "Other" = "Other",
      "Unknown" = "Unknown",
      "White" = "White",
      "F" = "Female",
      "M" = "Male"
    ),
    endpt_lbl = endpt
  )

fplt1<-forestplot.HR(
  df = dat %>% 
    filter(
      endpt %in% c("ACM","MACE") & 
        model =="main" & 
        grepl("^(0[1-9])+",stratum_var_lbl)
    ),  # long table
  x_idx1="stratum_var_lbl", # 1st layer index
  x_idx2="stratum_val_lbl", # 2nd layer index
  y_idx="endpt_lbl", # 1st layer y index
  est="mean", # estimates
  lower="lower", # 95% CI lower bound
  upper="upper", # 95% CI upper bound
  pval="pval", # p value
  plt_par = list(
    xlim = rep(list(c(0, 1.4)),2),
    vert_line = rep(list(c(0.3, 1.2)),2),
    ticks_at = rep(list(c(0.1, 0.5, 1, 1.2)),2)
  ), 
  ny = 2, # number of y groups (must be the same as groups of y_idx)
  idx_display = "Stratification",
)
# save figure
ggsave(
  file.path(path_to_outdir,"pap_expo_on_surv_mace_str_p1.tiff"),
  plot = fplt1,
  dpi = 150,
  width = 12, 
  height = 12, 
  units = "in",
  device = "tiff"
)

fplt2<-forestplot.HR(
  df = dat %>% 
    filter(
      endpt %in% c("ACM","MACE") & 
        model =="main" & 
        !grepl("^(0[1-9])+",stratum_var_lbl)
    ),  # long table
  x_idx1="stratum_var_lbl", # 1st layer index
  x_idx2="stratum_val_lbl", # 2nd layer index
  y_idx="endpt_lbl", # 1st layer y index
  est="mean", # estimates
  lower="lower", # 95% CI lower bound
  upper="upper", # 95% CI upper bound
  pval="pval", # p value
  plt_par = list(
    xlim = rep(list(c(0, 1.4)),2),
    vert_line = rep(list(c(0.3, 1.2)),2),
    ticks_at = rep(list(c(0.1, 0.5, 1, 1.2)),2)
  ), 
  ny = 2, # number of y groups (must be the same as groups of y_idx)
  idx_display = "Stratification",
)
# save figure
ggsave(
  file.path(path_to_outdir,"pap_expo_on_surv_mace_str_p2.tiff"),
  plot = fplt2,
  dpi = 150,
  width = 12, 
  height = 12, 
  units = "in",
  device = "tiff"
)

fplt1<-forestplot.HR(
  df = dat_all %>% filter(!endpt %in% c("ACM","MACE")),
  x_idx1="stratum_var", # 1st layer index
  x_idx2="stratum_val", # 2nd layer index
  y_idx="endpt", # 1st layer y index
  est="mean", # estimates
  lower="lower", # 95% CI lower bound
  upper="upper", # 95% CI upper bound
  pval="pval", # p value
  plt_par = list(
    xlim = rep(list(c(0, 1.4)),4),
    vert_line = rep(list(c(0.3, 1.2)),4),
    ticks_at = rep(list(c(0.1, 0.5, 1, 1.2)),4)
  ), 
  ny = 4, # number of y groups (must be the same as groups of y_idx)
  idx_display = "Full\nModel"
)
# save figure
ggsave(
  file.path(path_to_outdir,"pap_expo_on_mace_events_full.tiff"),
  plot = fplt1,
  dpi = 100,
  width = 18, 
  height = 5, 
  units = "in",
  device = "tiff"
)

fplt1<-forestplot.HR(
  df = dat %>% 
    filter(
      !endpt %in% c("ACM","MACE") &
        grepl("^(0[1-9])+",stratum_var_lbl)
    ),  # long table
  x_idx1="stratum_var_lbl", # 1st layer index
  x_idx2="stratum_val_lbl", # 2nd layer index
  y_idx="endpt_lbl", # 1st layer y index
  est="mean", # estimates
  lower="lower", # 95% CI lower bound
  upper="upper", # 95% CI upper bound
  pval="pval", # p value
  plt_par = list(
    xlim = rep(list(c(0, 1.4)),4),
    vert_line = rep(list(c(0.3, 1.2)),4),
    ticks_at = rep(list(c(0.1, 0.5, 1, 1.2)),4)
  ), 
  ny = 4, # number of y groups (must be the same as groups of y_idx)
  idx_display = "Stratification",
)
# save figure
ggsave(
  file.path(path_to_outdir,"pap_expo_on_mace_events_str_p1.tiff"),
  plot = fplt1,
  dpi = 150,
  width = 20,
  height = 12, 
  units = "in",
  device = "tiff"
)

fplt2<-forestplot.HR(
  df = dat %>% 
    filter(
      !endpt %in% c("ACM","MACE") &
        !grepl("^(0[1-9])+",stratum_var_lbl)
    ),  # long table
  x_idx1="stratum_var_lbl", # 1st layer index
  x_idx2="stratum_val_lbl", # 2nd layer index
  y_idx="endpt_lbl", # 1st layer y index
  est="mean", # estimates
  lower="lower", # 95% CI lower bound
  upper="upper", # 95% CI upper bound
  pval="pval", # p value
  plt_par = list(
    xlim = rep(list(c(0, 1.4)),4),
    vert_line = rep(list(c(0.3, 1.2)),4),
    ticks_at = rep(list(c(0.1, 0.5, 1, 1.2)),4)
  ), 
  ny = 4, # number of y groups (must be the same as groups of y_idx)
  idx_display = "Stratification",
)
# save figure
ggsave(
  file.path(path_to_outdir,"pap_expo_on_mace_events_str_p2.tiff"),
  plot = fplt2,
  dpi = 150,
  width = 20,
  height = 12, 
  units = "in",
  device = "tiff"
)

adh_metrics<-c(
  "CPAP_YR1"
  ,"adherence_yr1_med"
  ,"adherence_yr1_tt"
  ,"adherence_yr1_qt"
  ,"adherence_yr1_ww"
  ,"adherence_yr1_inc4"
  ,"adherence_yr1_inc8"
  ,"adherence_yr1_mttree"
  ,"adherence_yr1_mctree"
)

cpap_yr1_use<-readRDS(file.path(path_to_datadir,"cpap_adherence_final.rda")) %>%
  select(all_of(c("PATID",adh_metrics)))
N<-nrow(cpap_yr1_use)

ggplot(cpap_yr1_use,aes(x=CPAP_YR1))+
  geom_histogram(
    aes(y=..density..),
    binwidth = 1,
    fill="grey",
    color="blue")+
  geom_density(
    aes(y = ..density..),
    color = "red",
    stat = 'density',
    bw = 2,
    size = 1.5) +
  stat_bin(
    aes(y = ..density..,
        label=scales::percent(round(..density..,3))),
    geom='text',
    binwidth=1,
    vjust = 1.5,
    fontface = "bold"
  )+
  scale_x_continuous(
    name = "CPAP Total Charges in Year 1",
    breaks = 1:37,labels = 1:37)+
  scale_y_continuous(
    name = "Percentage (%)",
    breaks = seq(0,0.2,by=0.02),
    labels = seq(0,0.2,by=0.02)*100,
    sec.axis = sec_axis(
      trans = ~ . * N,
      name = "Count",
      breaks = seq(0,40000,by=2000),
      labels = seq(0,40000,by=2000))
  )+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"),
    panel.grid.major.y = element_line(colour = "grey"),
    text = element_text(face="bold",size=14)
  )
# save figure
ggsave(
  file.path(path_to_outdir,"pap_yr1_dist.tiff"),
  dpi = 100,
  width = 16, 
  height = 6, 
  units = "in",
  device = 'tiff'
)

adh_metrics_label<-c(
  "Raw",
  "Quantile_Median",
  "Quantile_Tertile",
  "Quantile_Quartile",
  "Rule_Based",
  "EqualSpaced_By4",
  "EqualSpaced_By8",
  "Empirical_All-cause-Mortality",
  "Empirical_MACE-composite"
)

plot_lst<-list()
# https://stackoverflow.com/questions/48212824/shaded-area-under-density-curve-in-ggplot2
for(i in seq_along(adh_metrics[-1])){
  adh_m<-adh_metrics[-1][i]
  adh_m_lbl<-adh_metrics_label[-1][i]
  bound<-cpap_yr1_use %>%
    group_by_at(vars(adh_m)) %>%
    summarise(lb=min(CPAP_YR1),
              ub=max(CPAP_YR1),
              .groups = "drop")
  xlabs <- c(bound$lb[1],bound$ub)
  ncut <- length(xlabs)
  annot <- paste0("lev",seq_along(xlabs)[-ncut])
  annot_x <- (xlabs[-ncut] + xlabs[-1])/2
  density_df <- density(cpap_yr1_use$CPAP_YR1,bw=2) %$% 
    data.frame(x = x, y = y) %>% 
    cross_join(bound) %>%
    filter(x<=ub & x>=lb-1) %>%
    rename("adh_cat"=adh_m)
  
  plot_lst[[adh_m]]<-ggplot(density_df,
                            aes(x = x, ymin = 0, ymax = y, fill = adh_cat))+
    geom_ribbon() +
    geom_line(aes(y = y)) +
    geom_vline(
      xintercept = xlabs,
      color="red", size=1,linetype=2) +
    annotate(
      geom = "text",
      label = xlabs,
      x = xlabs, 
      y = seq(0.01,0.08,length.out = ncut),
      fontface = "bold"
    ) +
    annotate(
      geom = "label",
      label = annot,
      x = annot_x, 
      y = rep(0.1,length(annot_x)),
      fontface = "bold"
    ) +
    theme(
      panel.background = element_rect(fill = "white", colour = "grey50"),
      panel.grid.major.y = element_line(colour = "grey"),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )+
    ggtitle(adh_m_lbl)
}
grid.arrange(grobs = plot_lst, ncol = 3)
# save figure
ggsave(
  file.path(path_to_outdir,"pap_yr1_discrete.tiff"),
  arrangeGrob(grobs = plot_lst, ncol = 3),
  dpi = 100,
  width = 12, 
  height = 8, 
  units = "in",
  device = 'tiff'
)

whole_dat<-c()
for(endpt_i in endpts){
  for(adh_i in adh_metrics){
    if(endpt_i=="ACM"&grepl("mctree",adh_i)) next
    if(endpt_i!="ACM"&grepl("mttree",adh_i)) next
    # endpt_i<-endpts[1]
    # adh_i<-adh_metrics[1]
    # boot_i<-1
    adh_df<-c()
    for(boot_i in 1:nboots){
      adh_boot<-readRDS(
        file.path(
          path_to_datadir,paste0("ADHRN_",endpt_i),adh_i,paste0("boot",boot_i),
          "coxph_iptw_main.rda")
      )$fit_cox
      
      if(adh_i=="CPAP_YR1"){
        adh_df %<>%
          bind_rows(
            data.frame(summary(adh_boot)$coefficients) %>%
              filter(!is.na(coef)) %>%
              mutate(`X2.5..`= coef-1.96*`se.coef.`,
                     `X97.5..`= coef+1.96*`se.coef.`,
                     `Pr...z..`= p) %>%
              rownames_to_column(var = "vari") %>%
              filter(grepl(adh_i,vari)) %>%
              gather(summ_var,summ_val,-vari)
          )
      }else{
        adh_df %<>%
          bind_rows(
            data.frame(summary(adh_boot)$coefficients) %>%
              cbind(data.frame(confint(adh_boot))) %>%
              rownames_to_column(var = "vari") %>%
              filter(grepl(adh_i,vari)) %>%
              gather(summ_var,summ_val,-vari)
          )
      }
    }
    whole_dat %<>%
      dplyr::bind_rows(
        adh_df %>%
          dplyr::group_by(vari,summ_var) %>%
          dplyr::summarize(summ_val_mu=mean(summ_val),
                           summ_val_sd=sd(summ_val),
                           .groups = "drop") %>%
          dplyr::filter(summ_var %in% c("coef","X2.5..","X97.5..","Pr...z..")) %>%
          dplyr::mutate(vari_cat = gsub(adh_i,"",vari),
                        vari = adh_i,
                        endpt = endpt_i)
      )
  }
}

whole_dat %<>%
  pivot_wider(
    id_cols = c("endpt","vari","vari_cat"),
    names_from = "summ_var",
    values_from = "summ_val_mu"
  )
colnames(whole_dat)<-c("endpt","vari","vari_cat","pval","ci_lb","ci_ub","coef")

whole_dat %<>%
  # fixed labels
  dplyr::mutate(
    vari_cat = case_when(grepl("med",vari) ~ "med_2",
                         TRUE ~ vari_cat)
  ) %>%
  # fit forestplot default keys
  dplyr::mutate(
    pvallabel = case_when(pval < 0.001 ~ "<0.001",
                          TRUE ~ as.character(round(pval,3))),
    mean = exp(coef),
    lower = exp(ci_lb),
    upper = exp(ci_ub),
    group =  case_when(grepl("pspline",vari_cat) ~ "cnt",
                       TRUE ~ gsub(".*_","lev",vari_cat)),
    proxy = recode(
      vari,
      CPAP_YR1 = "Charges_Raw",             
      adherence_yr1_med = "Quantile_Median",
      adherence_yr1_tt = "Quantile_Tertile",
      adherence_yr1_qt = "Quantile_Quartile",
      adherence_yr1_ww = "Rule_Based",
      adherence_yr1_inc4 = "EqualSpaced_By4",
      adherence_yr1_inc8 = "EqualSpaced_By8",
      adherence_yr1_mttree = "Empirical",
      adherence_yr1_mctree = "Empirical"
    ),
    endpt_lbl = endpt
  )

# ggplot(whole_dat %>% filter(endpt %in% c("ACM","mace")),
#        aes(color=vari,y=vari_cat))+
#   geom_point(aes(x=coef))+
#   geom_errorbar(aes(xmin=ci_lb,xmax=ci_ub))+
#   geom_vline(xintercept=1,linetype=2,color = "red") +
#   theme(
#         panel.background = element_rect(fill = "white", colour = "grey50"),
#         panel.grid.major.y = element_line(colour = "grey"),
#         legend.position = "none",
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank()
#       )+
#   facet_wrap(~ endpt, ncol = 2)

fplt<-forestplot.HR(
  df = whole_dat %>% filter(endpt %in% c("ACM","MACE")),  # long table
  x_idx1="proxy", # 1st layer index
  x_idx2="group", # 2nd layer index
  y_idx="endpt_lbl", # 1st layer y index
  est="mean", # estimates
  lower="lower", # 95% CI lower bound
  upper="upper", # 95% CI upper bound
  pval="pval", # p value
  plt_par = list(
    ci_column = c(2,4),
    ref_line = c(1,1),
    xlim = rep(list(c(0.6, 1.1)),2),
    ticks_at = rep(list(c(0.6,0.8,1)),2)
  ), 
  ny = 2, # number of y groups (must be the same as groups of y_idx)
  idx_display = "AdherenceProxy",
)
# save figure
ggsave(
  file.path(path_to_outdir,"pap_yr1_discrete_on_surv_mace.svg"),
  plot = fplt,
  dpi = 150,
  width = 12, 
  height = 8, 
  units = "in",
  device = 'svg'
)

# ggplot(whole_dat %>% filter(endpt %in% c("ACM","mace")),
#        aes(color=vari,y=vari_cat))+
#   geom_point(aes(x=coef))+http://127.0.0.1:23643/graphics/plot_zoom_png?width=1920&height=1017
#   geom_errorbar(aes(xmin=ci_lb,xmax=ci_ub))+
#   geom_vline(xintercept=1,linetype=2,color = "red") +
#   theme(
#         panel.background = element_rect(fill = "white", colour = "grey50"),
#         panel.grid.major.y = element_line(colour = "grey"),
#         legend.position = "none",
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank()
#       )+
#   facet_wrap(~ endpt, ncol = 2)

# fplt<-forestplot.HR(
#   df = whole_dat %>% filter(endpt %in% c("ACM","MACE") & proxy=="Quantile_Quartile"),  # long table
#   x_idx1="proxy", # 1st layer index
#   x_idx2="group", # 2nd layer index
#   y_idx="endpt_lbl", # 1st layer y index
#   est="mean", # estimates
#   lower="lower", # 95% CI lower bound
#   upper="upper", # 95% CI upper bound
#   pval="pval", # p value
#   plt_par = list(), # other plotting parameters passed in forest function
#   ny = 2, # number of y groups (must be the same as groups of y_idx)
#   idx_display = "AdherenceProxy",
# )
# # save figure
# ggsave(
#   file.path(path_to_outdir,"pap_yr1_quartile_on_surv_mace.tiff"),
#   plot = fplt,
#   dpi = 100,
#   width = 12, 
#   height = 8, 
#   units = "in",
#   device = 'tiff'
# )
# 
# fplt<-forestplot.HR(
#   df = whole_dat %>% filter(!endpt %in% c("ACM","MACE")),  # long table
#   x_idx1="proxy", # 1st layer index
#   x_idx2="group", # 2nd layer index
#   y_idx="endpt_lbl", # 1st layer y index
#   est="mean", # estimates
#   lower="lower", # 95% CI lower bound
#   upper="upper", # 95% CI upper bound
#   pval="pval", # p value
#   plt_par = list(), # other plotting parameters passed in forest function
#   ny = 4, # number of y groups (must be the same as groups of y_idx)
#   idx_display = "AdherenceProxy",
# )
# # save figure
# ggsave(
#   file.path(path_to_dir,"res","cpap_yr1_discrete_on_mace_events.tiff"),
#   plot = fplt,
#   dpi = 150,
#   width = 20, 
#   height = 8, 
#   units = "in",
#   device = 'tiff'
# )
# 
# fplt<-forestplot.HR(
#   df = whole_dat %>% filter(!endpt %in% c("ACM","MACE") & proxy=="Quantile_Quartile"),  # long table
#   x_idx1="proxy", # 1st layer index
#   x_idx2="group", # 2nd layer index
#   y_idx="endpt_lbl", # 1st layer y index
#   est="mean", # estimates
#   lower="lower", # 95% CI lower bound
#   upper="upper", # 95% CI upper bound
#   pval="pval", # p value
#   plt_par = list(), # other plotting parameters passed in forest function
#   ny = 4, # number of y groups (must be the same as groups of y_idx)
#   idx_display = "AdherenceProxy",
# )
# # save figure
# ggsave(
#   file.path(path_to_dir,"res","cpap_yr1_quartile_on_mace_events.tiff"),
#   plot = fplt,
#   dpi = 100,
#   width = 20, 
#   height = 8, 
#   units = "in",
#   device = 'tiff'
# )

endpt_lbl<-c(
  "1.All Cause Mortality",
  "2.MACE Composite",
  "3.Myocardial Infarction",
  "4.Heart Failure",
  "5.Stroke",
  "6.Revascularization"
)
raw_plt_lst<-list()
for(i in seq_along(endpts)){
  # i<-1
  endpt_i<-endpts[i]
  whole_dat_raw<-c()
  turnpt<-c()
  altpt<-c()
  for(boot_i in 1:nboots){
    # spline plot
    tp<-readRDS(
      file.path(path_to_datadir,paste0("ADHRN_",endpt_i),"CPAP_YR1",paste0("boot",boot_i),
                "pspline_termplot.rda"))     
    whole_dat_raw %<>%
      bind_rows(tp %>% mutate(boot = boot_i,endpt = endpt_i))
    
    # turning points
    delta_y<-diff(tp$y)
    turns_ind <- which(delta_y[-1] * delta_y[-length(delta_y)] < 0) + 1
    x_turns <- tp$x[turns_ind]
    turnpt<-unique(c(turnpt,x_turns))
    # sign-altering point
    signalt_ind<-which(tp$y[-1] * tp$y[-length(tp$y)] < 0) + 1
    x_alts <- tp$x[signalt_ind]
    altpt<-unique(c(altpt,x_alts))
  }
  # plot
  whole_dat_raw %<>%
    group_by(x) %>% summarise(y=mean(y),se=mean(se),.groups="drop")
  raw_plt_lst[[endpt_lbl[i]]]<-ggplot(whole_dat_raw,aes(x=x,y=exp(y))) +
    geom_ribbon(
      aes(ymin = exp(y - 1.96*se),ymax = exp(y + 1.96*se)),
      size = 1.5,fill = "grey70") +
    geom_line(size = 1.5) +
    geom_vline(
      data=data.frame(x=altpt),
      aes(xintercept=x),
      color = "blue",linetype=2,size=1
    ) + 
    annotate(
      geom = "label",
      label = turnpt,
      x = turnpt,
      y = rep(2,length(turnpt)),
      size = 5,
      fontface="bold"
    ) +
    annotate(
      geom = "label",
      label = altpt,
      x = altpt, 
      y = rep(0.5,length(altpt)),
      size = 5,
      fontface="bold",
      color = "red"
    ) +
    geom_hline(aes(yintercept=1),
               color = "red",linetype=2,size=1) + 
    labs(x = "Raw Year-1 Total Charges",
         y = "HR") +
    scale_y_continuous(
      limits = c(0,3),
      breaks = seq(0,3,0.5),
      labels = seq(0,3,0.5)) + 
    theme(
      panel.background = element_rect(fill = "white", colour = "grey50"),
      panel.grid.major.y = element_line(colour = "grey"),
      legend.position = "none",
      text = element_text(face ="bold")
    ) +
    ggtitle(endpt_lbl[[i]])
}
grid.arrange(grobs = raw_plt_lst, ncol = 2)
# save figure
ggsave(
  file.path(path_to_outdir,"pap_yr1_raw_on_surv_mace.tiff"),
  arrangeGrob(grobs = raw_plt_lst, ncol = 2),
  dpi = 100,
  width = 10, 
  height = 8, 
  units = "in",
  device = 'tiff'
)

strat_dat<-c()
# stratefied model
for(endpt_i in endpts){
  for(adh_i in adh_metrics){
    if(endpt_i=="ACM"&grepl("mctree",adh_i)) next
    if(endpt_i!="ACM"&grepl("mttree",adh_i)) next
    # endpt_i<-endpts[1]
    adh_i<-"adherence_yr1_qt"
    # boot_i<-1
    for(boot_i in 1:nboots){
      # non-stratified model
      adh_boot<-readRDS(
        file.path(path_to_datadir,paste0("ADHRN_",endpt_i),adh_i,paste0("boot",boot_i),
                  "coxph_iptw_main.rda")
      )$fit_cox 
      strat_dat %<>%
        bind_rows(
          data.frame(summary(adh_boot)$coefficients) %>%
            filter(!is.na(coef)) %>%
            mutate(across(everything(), as.character)) %>%
            rownames_to_column(var = "fit_var") %>%
            mutate(stratum_var = "None",
                   stratum_val = " ",
                   endpt = endpt_i,
                   boot = boot_i)
        )
      # stratified model
      adh_boot<-readRDS(
        file.path(path_to_datadir,paste0("ADHRN_",endpt_i),adh_i,paste0("boot",boot_i),
                  "coxph_strata_main.rda") 
      ) %>% as.data.frame()
      strat_dat %<>%
        bind_rows(
          data.frame(adh_boot) %>% 
            filter(!is.na(coef)) %>%
            mutate(endpt = endpt_i,
                   boot = boot_i)
        )
    }
  }
}

strat_dat %<>%
  # filter(stratum_var != "MACE_HISTORY") %>%
  group_by(stratum_var,stratum_val,fit_var,endpt) %>%
  summarize(
    mean = mean(as.numeric(coef),na.rm=T),
    sd = mean(as.numeric(`se.coef.`),na.rm=T),
    pval = mean(as.numeric(`Pr...z..`),na.rm=T),
    .groups = "drop"
  ) %>%
  mutate(
    lower = mean - 1.96*sd,
    upper = mean + 1.96*sd) %>%
  mutate(
    mean = exp(mean),
    lower = exp(lower),
    upper = exp(upper)
  )

strat_dat$endpt<-factor(strat_dat$endpt, levels=c("ACM","MACE","MI","HF","STROKE","REVASC"))
strat_dat %<>% arrange(endpt) %>%
  filter(grepl("adherence_yr1_qtqt",fit_var)) %>%
  mutate(
    fit_var_lbl = recode(
      fit_var,
      "adherence_yr1_qtqt_2" = "Q2 (8 - 12)",
      "adherence_yr1_qtqt_3" = "Q3 (13 - 15)",
      "adherence_yr1_qtqt_4" = "Q4 (> 15)"
    ),
  )

strat_dat_full<-strat_dat %>% 
  filter(stratum_var=="None") %>%
  mutate(stratum_var="")
  
fplt1<-forestplot.HR(
  df = strat_dat_full %>% filter(endpt %in% c("ACM","MACE")),  # long table
  x_idx1="stratum_var", # 1st layer index
  x_idx2="fit_var_lbl", # 2nd layer index
  y_idx="endpt", # 1st layer y index
  est="mean", # estimates
  lower="lower", # 95% CI lower bound
  upper="upper", # 95% CI upper bound
  pval="pval", # p value
  plt_par = list( # manual adjustment
    xlim = rep(list(c(0, 1.5)),2),
    vert_line = rep(list(c(0.3, 1.2)),2),
    ticks_at = rep(list(c(0.1, 0.5, 1, 1.2)),2)
  ), 
  ny = 2, # number of y groups (must be the same as groups of y_idx)
  idx_display = "Full\nModel",
)
# save figure
ggsave(
  file.path(path_to_outdir,"pap_adh_on_surv_mace_full.tiff"),
  plot = fplt1,
  dpi = 100,
  width = 12, 
  height = 6, 
  units = "in",
  device = 'tiff'
)

strat_dat2<-strat_dat %>% 
  filter(stratum_var!="None") %>%
  mutate(
    stratum_var_lbl = recode(
      stratum_var,
      "AGEGRP" = "01.Age",
      "SEX" = "02.Sex",
      "RACE_LABEL" = "03.Race",
      "LIS_DUAL_IND" = "04.Low-income-subsidy/Dual Eligibility",
      "HYSM" = "05.Hypersomnia",
      "INSM" = "06.Insomina",
      "MACE_HIST" = "07.MACE History",
      "OBES" = "08.Obesity",
      "COPD" = "09.COPD",
      "HTN" = "10.Hypertension",
      "T2DM" = "11.T2DM", 
      "ND" = "12.Anxiety Disorder",
      "AFIB" = "13.Atrial Fibrilation",
      "CKD" = "14.Chronic Kidney Disease",
      "CCI_CLASS" = "15.Charlson Comorbidity Index",
      "ACG" = "16.Anti-coagulant",
      "AHT" = "17.Anti-hypertensive",
      "ALP" = "18.Anti-lipemic",
      "BGR" = "19.Blood glucose regulator",
    ),
    stratum_val_lbl = recode(
      stratum_val,
      "agegrp1" = "65-69 years",
      "agegrp2" = "70-74 years",
      "agegrp3" = "75-79 years",
      "agegrp4" = "80+ years",
      "0" = "No",
      "1" = "Yes",
      " " = " ",
      "AA" = "Black",
      "AI" = "Native American",
      "Asian" = "Asian",
      "Other" = "Other",
      "Unknown" = "Unknown",
      "White" = "White",
      "F" = "Female",
      "M" = "Male"
    ),
    endpt_lbl = endpt
  ) %>%
  mutate(
    stratum_combine = paste0(stratum_var_lbl,":",stratum_val_lbl)
  )

fplt1<-forestplot.HR(
  df = strat_dat2 %>% 
    filter(endpt %in% c("ACM","MACE") & grepl("^(0[1-5])+",stratum_var_lbl)),  # long table
  x_idx1="stratum_combine", # 1st layer index
  x_idx2="fit_var_lbl", # 2nd layer index
  y_idx="endpt_lbl", # 1st layer y index
  est="mean", # estimates
  lower="lower", # 95% CI lower bound
  upper="upper", # 95% CI upper bound
  pval="pval", # p value
  plt_par = list( # manual adjustment
    xlim = rep(list(c(0, 1.5)),2),
    vert_line = rep(list(c(0.3, 1.2)),2),
    ticks_at = rep(list(c(0.1, 0.5, 1, 1.2)),2)
  ), 
  ny = 2, # number of y groups (must be the same as groups of y_idx)
  idx_display = "Stratification",
)
# save figure
ggsave(
  file.path(path_to_outdir,"pap_adh_on_surv_mace_str_p1.tiff"),
  plot = fplt1,
  dpi = 250,
  width = 12, 
  height = 18, 
  units = "in",
  device = 'tiff'
)

fplt2<-forestplot.HR(
  df = strat_dat2 %>% 
    filter(endpt %in% c("ACM","MACE") & 
             (grepl("^(0[6-9])+",stratum_var_lbl)|grepl("^(1[0-3])+",stratum_var_lbl))),  # long table
  x_idx1="stratum_combine", # 1st layer index
  x_idx2="fit_var_lbl", # 2nd layer index
  y_idx="endpt_lbl", # 1st layer y index
  est="mean", # estimates
  lower="lower", # 95% CI lower bound
  upper="upper", # 95% CI upper bound
  pval="pval", # p value
  plt_par = list(
    xlim = rep(list(c(0, 1.5)),2),
    vert_line = rep(list(c(0.3, 1.2)),2),
    ticks_at = rep(list(c(0.1, 0.5, 1, 1.2)),2)
  ), 
  ny = 2, # number of y groups (must be the same as groups of y_idx)
  idx_display = "Stratification",
)
# fplt2
# save figure
ggsave(
  file.path(path_to_outdir,"pap_adh_on_surv_mace_str_p2.tiff"),
  plot = fplt2,
  dpi = 250,
  width = 12, 
  height = 18,
  units = "in",
  device = "tiff"
)

fplt3<-forestplot.HR(
  df = strat_dat2 %>% 
    filter(endpt %in% c("ACM","MACE") & 
             (grepl("^(1[4-9])+",stratum_var_lbl)|grepl("^(2[0-3])+",stratum_var_lbl))),  # long table
  x_idx1="stratum_combine", # 1st layer index
  x_idx2="fit_var_lbl", # 2nd layer index
  y_idx="endpt_lbl", # 1st layer y index
  est="mean", # estimates
  lower="lower", # 95% CI lower bound
  upper="upper", # 95% CI upper bound
  pval="pval", # p value
  plt_par = list(
    xlim = rep(list(c(0, 1.5)),2),
    vert_line = rep(list(c(0.3, 1.2)),2),
    ticks_at = rep(list(c(0.1, 0.5, 1, 1.2)),2)
  ), 
  ny = 2, # number of y groups (must be the same as groups of y_idx)
  idx_display = "Stratification",
)
# fplt3
# save figure
ggsave(
  file.path(path_to_outdir,"pap_adh_on_surv_mace_str_p3.tiff"),
  plot = fplt3,
  dpi = 250,
  width = 12, 
  height = 18,
  units = "in",
  device = "tiff"
)

fplt1<-forestplot.HR(
  df = strat_dat_full %>% filter(!endpt %in% c("ACM","MACE")),  # long table
  x_idx1="stratum_var", # 1st layer index
  x_idx2="fit_var_lbl", # 2nd layer index
  y_idx="endpt", # 1st layer y index
  est="mean", # estimates
  lower="lower", # 95% CI lower bound
  upper="upper", # 95% CI upper bound
  pval="pval", # p value
  plt_par = list( # manual adjustment
    xlim = rep(list(c(0, 1.5)),4),
    vert_line = rep(list(c(0.3, 1.2)),4),
    ticks_at = rep(list(c(0.1, 0.5, 1, 1.2)),4)
  ), 
  ny = 4, # number of y groups (must be the same as groups of y_idx)
  idx_display = "Full\nModel",
)
# save figure
ggsave(
  file.path(path_to_outdir,"pap_adh_on_mace_events_full.tiff"),
  plot = fplt1,
  dpi = 100,
  width = 18, 
  height = 6, 
  units = "in",
  device = 'tiff'
)

fplt1<-forestplot.HR(
  df = strat_dat2 %>% 
    filter(stratum_var != "MACE_HIST") %>%
    filter(!endpt %in% c("ACM","MACE") & grepl("^(0[1-5])+",stratum_var_lbl)),  # long table
  x_idx1="stratum_combine", # 1st layer index
  x_idx2="fit_var_lbl", # 2nd layer index
  y_idx="endpt_lbl", # 1st layer y index
  est="mean", # estimates
  lower="lower", # 95% CI lower bound
  upper="upper", # 95% CI upper bound
  pval="pval", # p value
  plt_par = list( # manual adjustment
    xlim = rep(list(c(0, 1.5)),4),
    vert_line = rep(list(c(0.3, 1.2)),4),
    ticks_at = rep(list(c(0.1, 0.5, 1, 1.2)),4)
  ), 
  ny = 4, # number of y groups (must be the same as groups of y_idx)
  idx_display = "Stratification",
)
# fplt1
# save figure
ggsave(
  file.path(path_to_outdir,"pap_adh_on_mace_event_str_p1.tiff"),
  plot = fplt1,
  dpi = 250,
  width = 18, 
  height = 18, 
  units = "in",
  device = "tiff"
)

fplt2<-forestplot.HR(
  df = strat_dat2 %>% 
    filter(stratum_var != "MACE_HIST" & !endpt %in% c("ACM","MACE") &
            (grepl("^(0[6-9])+",stratum_var_lbl)|grepl("^(1[0-3])+",stratum_var_lbl))),  # long table
  x_idx1="stratum_combine", # 1st layer index
  x_idx2="fit_var_lbl", # 2nd layer index
  y_idx="endpt_lbl", # 1st layer y index
  est="mean", # estimates
  lower="lower", # 95% CI lower bound
  upper="upper", # 95% CI upper bound
  pval="pval", # p value
  plt_par = list(
    xlim = rep(list(c(0, 1.5)),4),
    vert_line = rep(list(c(0.3, 1.2)),4),
    ticks_at = rep(list(c(0.1, 0.5, 1, 1.2)),4)
  ), 
  ny = 4, # number of y groups (must be the same as groups of y_idx)
  idx_display = "Stratification",
)
# fplt2
# save figure
ggsave(
  file.path(path_to_outdir,"pap_adh_on_mace_event_str_p2.tiff"),
  plot = fplt2,
  dpi = 250,
  width = 18, 
  height = 18, 
  units = "in",
  device = "tiff"
)

fplt3<-forestplot.HR(
  df = strat_dat2 %>% 
    filter(stratum_var != "MACE_HIST" & !endpt %in% c("ACM","MACE") & 
            (grepl("^(1[4-9])+",stratum_var_lbl)|grepl("^(2[0-3])+",stratum_var_lbl))),  # long table
  x_idx1="stratum_combine", # 1st layer index
  x_idx2="fit_var_lbl", # 2nd layer index
  y_idx="endpt_lbl", # 1st layer y index
  est="mean", # estimates
  lower="lower", # 95% CI lower bound
  upper="upper", # 95% CI upper bound
  pval="pval", # p value
  plt_par = list(
    xlim = rep(list(c(0, 1.5)),4),
    vert_line = rep(list(c(0.3, 1.2)),4),
    ticks_at = rep(list(c(0.1, 0.5, 1, 1.2)),4)
  ), 
  ny = 4, # number of y groups (must be the same as groups of y_idx)
  idx_display = "Stratification",
)
# fplt2
# save figure
ggsave(
  file.path(path_to_outdir,"pap_adh_on_mace_event_str_p3.tiff"),
  plot = fplt3,
  dpi = 250,
  width = 18, 
  height = 18, 
  units = "in",
  device = "tiff"
)

