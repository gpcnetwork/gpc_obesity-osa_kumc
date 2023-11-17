/*
# Copyright (c) 2021-2025 University of Missouri                   
# Author: Xing Song, xsm7f@umsystem.edu                            
# File: osa-outcome-wide.sql
# Description: collecting all phecodes for OSA population
*/

use role GROUSE_ROLE_C_ANALYTICS_DROC114; 
use warehouse GROUSE_WH;
use database GROUSE_DEID_ANALYTICS_DB_DROC114;
create schema if not exists SX_OSA_OUTCOME_WIDE;
use schema SX_OSA_OUTCOME_WIDE;

set cdm_db_schema = 'CMS_PCORNET_CDM';
set diagnosis = $cdm_db_schema || '.V_DEID_DIAGNOSIS';
set procedures = $cdm_db_schema || '.V_DEID_PROCEDURES';
set enrollment = $cdm_db_schema || '.V_DEID_ENROLLMENT';
set demographic = $cdm_db_schema || '.V_DEID_DEMOGRAPHIC';
set death = $cdm_db_schema || '.V_DEID_DEATH';
set dispensing = $cdm_db_schema || '.V_DEID_DISPENSING';

create or replace table PAT_OSA_ALL_PHECD as 
with phecd_map_cte as (
       select distinct dx.*,
              phe."phecode" as phecd_dxgrpcd, 
              ref."phecode_string" as phecd_dxgrp
       from SX_OSA_PAP.PAT_OSA_ALL_DX dx 
       join ONTOLOGY.GROUPER_VALUESETS.ICD10CM_PHECODEX phe on dx.DX = phe."icd10" and dx.DX_TYPE = '10'
              and phe."phecode" not like '%.%' -- only keep the integers
       join ONTOLOGY.GROUPER_VALUESETS.PHECODEX_REF ref on phe."phecode" = ref."phecode" 
       union
       select distinct dx.*,
              phe."phecode" as phecd_dxgrpcd, 
              ref."phecode_string" as phecd_dxgrp
       from SX_OSA_PAP.PAT_OSA_ALL_DX dx 
       join ONTOLOGY.GROUPER_VALUESETS.ICD9CM_PHECODEX phe on dx.DX = phe."icd9" and dx.DX_TYPE = '09'
              and phe."phecode" not like '%.%'  -- only keep the integers
       join ONTOLOGY.GROUPER_VALUESETS.PHECODEX_REF ref on phe."phecode" = ref."phecode"
), phecd_ord as (
    select  patid,
            phecd_dxgrpcd,
            phecd_dxgrp,
            dx_date as phecd_date,
            days_since_index,
            row_number() over (partition by patid order by days_since_index) rn_asc,
            row_number() over (partition by patid order by days_since_index desc) rn_desc
    from phecd_map_cte 
)
select patid,
       phecd_dxgrpcd,
       phecd_dxgrp,
       phecd_date,
       days_since_index
from phecd_ord
where rn_asc = 1 or rn_desc = 1       
;

select count(distinct patid) from PAT_OSA_ALL_PHECD;
-- 888,845

with cte_summ as (
    select phecd_dxgrp, count(distinct patid) as ptcnt, 'bef' as timepoint
    from PAT_OSA_ALL_PHECD
    where days_since_index <= 0
    group by phecd_dxgrp
    union
    select phecd_dxgrp, count(distinct patid) as ptcnt, 'aft' as timepoint
    from PAT_OSA_ALL_PHECD
    where days_since_index  > 0
    group by phecd_dxgrp
)
select phecd_int,
       bef as bef_cnt,
       round(bef/888845,3) as bef_prop,
       aft as aft_cnt,
       round(aft/888845,3) as aft_prop
from cte_summ
pivot (
    max(ptcnt) for timepoint in ('bef','aft')
)
as p(phecd_int,bef,aft)
where bef is not null and aft is not null
order by bef desc
;
