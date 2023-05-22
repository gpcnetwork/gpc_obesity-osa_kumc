#################################################################
# Copyright (c) 2021-2025 University of Missouri                   
# Author: Xing Song, xsm7f@umsystem.edu                            
# File: send_valueset.R
# Description: string-search published or curated valuesets for
#              targeted valueset of interest
# Note: this script can be run from any machine, but may need to 
#       modify where to obtain snowflake credentials
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
  dbplyr,
  devtools,
  stringdist,
  jsonlite
  )

source_url("https://raw.githubusercontent.com/sxinger/utils/master/extract_util.R")

tgt_schema <- "OBESITY_OSA"
tgt_tbl <- "OSA_CPAP_COV_VS"

# make db connection
sf_conn <- DBI::dbConnect(drv = odbc::odbc(),
                          dsn = Sys.getenv("ODBC_DSN_NAME"),
                          uid = Sys.getenv("SNOWFLAKE_USER"),
                          pwd = Sys.getenv("SNOWFLAKE_PWD"))

##==== curated valuesets
cov_vec<-c(
  "bariatric surgery",
  "chronic obstructive pulmonary disease",
  # "neurotic disorders",
  "hypersomnia",
  "insomina"
)

for (i in seq_along(cov_vec)){
  load_valueset(vs_template = "curated",
                vs_url = "https://raw.githubusercontent.com/RWD2E/phecdm/main/res/valueset_curated/vs-osa-comorb.json",
                vs_name_str = cov_vec[i],
                dry_run = FALSE,
                conn=sf_conn,
                write_to_schema = tgt_schema,
                write_to_tbl = tgt_tbl,
                overwrite = (i == 1))
}

##==== ecqm valuesets
load_valueset(vs_template = "ecqm",
              vs_url = "https://raw.githubusercontent.com/RWD2E/phecdm/main/res/valueset_autogen/ecqm-condition-diagnosis-problem.json",
              vs_name_str = "kidney failure",
              dry_run = FALSE,
              conn=sf_conn,
              write_to_schema = tgt_schema,
              write_to_tbl = tgt_tbl,
              overwrite = FALSE)

load_valueset(vs_template = "ecqm",
              vs_url = "https://raw.githubusercontent.com/RWD2E/phecdm/main/res/valueset_autogen/ecqm-medication.json",
              vs_name_str = "antipsychotic",
              dry_run = FALSE,
              conn=sf_conn,
              write_to_schema = tgt_schema,
              write_to_tbl = tgt_tbl,
              overwrite = FALSE)

load_valueset(vs_template = "ecqm",
              vs_url = "https://raw.githubusercontent.com/RWD2E/phecdm/main/res/valueset_autogen/ecqm-medication.json",
              vs_name_str = "Pharmacologic Therapy for Hypertension",
              dry_run = FALSE,
              conn=sf_conn,
              write_to_schema = tgt_schema,
              write_to_tbl = tgt_tbl,
              overwrite = FALSE)
