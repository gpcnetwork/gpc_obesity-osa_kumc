#################################################################
# Copyright (c) 2021-2025 University of Missouri                   
# Author: Xing Song, xsm7f@umsystem.edu                            
# File: extract.R
# Description: pre-analytic data from snowflake to R
#################################################################

rm(list=ls()); gc()
setwd("C:/repo/GPC-Analytics-Weight-Cohort")

# install.packages("pacman")
pacman::p_load(
  DBI,
  jsonlite,
  odbc,
  tidyverse,
  magrittr,
  dbplyr
)

path_to_data_folder<-file.path(
  getwd(),
  "osa-obesity-cpap-adherence",
  "data"
)

# make db connection
sf_conn <- DBI::dbConnect(drv = odbc::odbc(),
                         dsn = Sys.getenv("ODBC_DSN_NAME"),
                         uid = Sys.getenv("SNOWFLAKE_USER"),
                         pwd = Sys.getenv("SNOWFLAKE_PWD"))

# collect cohort for exposure analysis
dat<-tbl(sf_conn,in_schema("OBESITY_OSA","PAT_OSA_CPAP_DEID")) %>% collect()
saveRDS(dat,file=file.path(path_to_data_folder,"cpap_exposure_tbl1.rds"))

dat_add<-tbl(sf_conn,in_schema("OBESITY_OSA","PAT_OSA_CPAP_DEID_ADDCOV")) %>% collect()
saveRDS(dat_add,file=file.path(path_to_data_folder,"cpap_exposure_cov.rds"))


# collect cohort for dose-response analysis
dat2<-tbl(sf_conn,in_schema("OBESITY_OSA","PAT_OSA_CPAP_SUPPLY_DEID")) %>% collect()
saveRDS(dat2,file=file.path(path_to_data_folder,"cpap_dose_response_tbl1.rds"))

dat2_add<-tbl(sf_conn,in_schema("OBESITY_OSA","PAT_OSA_CPAP_SUPPLY_DEID_ADDCOV")) %>% collect()
saveRDS(dat2_add,file=file.path(path_to_data_folder,"cpap_dose_response_cov.rds"))


# close db connection
DBI::dbDisconnect(sf_conn)
