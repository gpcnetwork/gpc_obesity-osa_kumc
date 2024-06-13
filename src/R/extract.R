#################################################################
# Copyright (c) 2021-2025 University of Missouri                   
# Author: Xing Song, xsm7f@umsystem.edu                            
# File: extract.R
# Description: pre-analytic data from snowflake to R
#################################################################
rm(list=ls()); gc()

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
  "data"
)

# make db connection
sf_conn <- DBI::dbConnect(
  drv = odbc::odbc(),
  dsn = Sys.getenv("ODBC_DSN_NAME"),
  uid = Sys.getenv("SNOWFLAKE_USER"),
  pwd = Sys.getenv("SNOWFLAKE_PWD")
)

# collect cohort for exposure analysis
dat<-tbl(sf_conn,in_schema("SX_OSA_PAP","PAT_OSA_PAP_EXPOS")) %>% collect()
saveRDS(dat,file=file.path(path_to_data_folder,"pap_exposure_aset.rda"))

# collect cohort for dose-response analysis
dat<-tbl(sf_conn,in_schema("SX_OSA_PAP","PAT_OSA_PAP_ADHRN")) %>% collect()
saveRDS(dat,file=file.path(path_to_data_folder,"cpap_adherence_aset.rda"))

# close db connection
DBI::dbDisconnect(sf_conn)
