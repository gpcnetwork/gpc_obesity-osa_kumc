/*
# Copyright (c) 2021-2025 University of Missouri                   
# Author: Xing Song, xsm7f@umsystem.edu                            
# File: osa-cpap-add.sql
# Description: collect additional covariates
# Dependency: existence of table identifier($cov_vs)
#             run send_valueset.R to create
*/

use role GROUSE_ROLE_C_ANALYTICS; 
use warehouse GROUSE_WH;
use database GROUSE_DEID_ANALYTICS_DB;
use schema OBESITY_OSA;

set pat_cohort1 = 'PAT_OSA_ELIG';
set pat_cohort2 = 'PAT_OSA_CPAP_DEVICE';
set cdm_db_schema = 'GROUSE_DEID_DB.CMS_PCORNET_CDM';
set diagnosis = $cdm_db_schema || '.V_DEID_DIAGNOSIS';
set procedures = $cdm_db_schema || '.V_DEID_PROCEDURES';
set dispensing = $cdm_db_schema || '.V_DEID_DISPENSING';
set address_history = $cdm_db_schema || '.V_DEID_ADDRESS_HISTORY';
set address_geocode = $cdm_db_schema || '.V_DEID_ADDRESS_GEOCODE';
set obs_comm = $cdm_db_schema || '.V_DEID_OBS_COMM';
set cov_vs = 'OSA_CPAP_COV_VS';

-----------------------------------------------------------------------------------

/*check valuesets*/
select code_grp, code_type, code_subtype, count(distinct code) 
from identifier($cov_vs)
group by code_grp, code_type, code_subtype
order by code_grp, code_type, code_subtype;

/*map RXNORM to NDC codes for med*/
-- ontology.rxnorm.rxnsat
create or replace table VS_NDC as
select distinct
       rxn.CODE as RXCUI
      ,rxmap.ATV as NDC
      ,rxmap.SUPPRESS
from identifier($cov_vs) rxn
join ontology.rxnorm.rxnsat rxmap
on rxn.CODE = rxmap.RXCUI and rxn.CODE_TYPE = 'RXNORM' and
   rxmap.ATN = 'NDC'and rxmap.SAB = 'RXNORM' -- normalized 11-digit NDC codes
;

/*collect additional covariates for $PAT_COHORT*/
create or replace table PAT_OSA_CPAP_DEID_ADDCOV as
with add_dx as (
    select a.PATID, datediff(day,a.OSA_DX1_DATE, b.DX_DATE) as DAYS_OSA_TO_COV, vs.CODE_GRP as COV
    from identifier($pat_cohort1) a
    join identifier($DIAGNOSIS) b  on a.PATID = b.PATID
    join identifier($cov_vs) vs on vs.CODE_TYPE_CDM = b.DX_TYPE and vs.CODE = b.DX and VS.CODE_TYPE like 'icd%cm'
),   add_px as (
    select a.PATID, datediff(day,a.OSA_DX1_DATE, b.PX_DATE) as DAYS_OSA_TO_COV, vs.CODE_GRP as COV
    from identifier($pat_cohort1) a
    join identifier($PROCEDURES) b  on a.PATID = b.PATID
    join identifier($cov_vs) vs on vs.CODE_TYPE_CDM = b.PX_TYPE and vs.CODE = b.PX and vs.CODE_TYPE in ('icd9-cm','icd10-pcs','hcpcs')
-- ),   add_med as (
--     select a.PATID, datediff(day,a.OSA_DX1_DATE, b.DX_DATE) as DAYS_OSA_TO_COV, c.CODE_GRP as COV
--     from identifier($pat_cohort1) a
--     join identifier($DIAGNOSIS) b  on a.PATID = b.PATID
--     join VS_NDC vs on vs.CODE = b.NDC
),   add_all as (
    select * from add_dx
    union all
    select * from add_px
    -- union all 
    -- select * from add_med
),   add_all_reduce as (
    select PATID, COV, DAYS_OSA_TO_COV,
           row_number() over (partition by PATID, COV order by DAYS_OSA_TO_COV) rn
    from add_all
)
select PATID, COV, DAYS_OSA_TO_COV
from add_all_reduce 
where rn = 1
;


create or replace table PAT_OSA_CPAP_SUPPLY_DEID_ADDCOV as
with add_dx as (
    select a.PATID, datediff(day,a.OSA_DX1_DATE, b.DX_DATE) as DAYS_OSA_TO_COV, vs.CODE_GRP as COV
    from identifier($pat_cohort2) a
    join identifier($DIAGNOSIS) b on a.PATID = b.PATID
    join identifier($cov_vs) vs on vs.CODE_TYPE_CDM = b.DX_TYPE and vs.CODE = b.DX and VS.CODE_TYPE like 'icd%cm'
),   add_px as (
    select a.PATID, datediff(day,a.OSA_DX1_DATE, b.PX_DATE) as DAYS_OSA_TO_COV, vs.CODE_GRP as COV
    from identifier($pat_cohort2) a
    join identifier($PROCEDURES) b  on a.PATID = b.PATID
    join identifier($cov_vs) vs on vs.CODE_TYPE_CDM = b.PX_TYPE and vs.CODE = b.PX and vs.CODE_TYPE in ('icd9-cm','icd10-pcs','hcpcs')
-- ),   add_med as (
--     select a.PATID, datediff(day,a.OSA_DX1_DATE, b.DX_DATE) as DAYS_OSA_TO_COV, c.CODE_GRP as COV
--     from identifier($pat_cohort2) a
--     join identifier($DIAGNOSIS) b  on a.PATID = b.PATID
--     join VS_NDC vs on vs.CODE = b.NDC
),   add_all as (
    select * from add_dx
    union all
    select * from add_px
    -- union all 
    -- select * from add_med
),   add_all_reduce as (
    select PATID, COV, DAYS_OSA_TO_COV,
           row_number() over (partition by PATID, COV order by DAYS_OSA_TO_COV) as rn
    from add_all
)
select PATID, COV, DAYS_OSA_TO_COV 
from add_all_reduce 
where rn = 1
;



