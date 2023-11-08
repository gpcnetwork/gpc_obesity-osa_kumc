/*
# Copyright (c) 2021-2025 University of Missouri                   
# Author: Xing Song, xsm7f@umsystem.edu                            
# File: osa-cpap-extract.sql
# Description: extracting osa population from Medicare claims and identify cpap exposure and study endpoints, 
#              as well as a few pre-defined covariates
*/

use role GROUSE_ROLE_C_ANALYTICS_DROC114; 
use warehouse GROUSE_WH;
use database GROUSE_DEID_ANALYTICS_DB_DROC114;
create schema if not exists SX_OSA_PAP;
use schema SX_OSA_PAP;

set cdm_db_schema = 'CMS_PCORNET_CDM';
set diagnosis = $cdm_db_schema || '.V_DEID_DIAGNOSIS';
set procedures = $cdm_db_schema || '.V_DEID_PROCEDURES';
set enrollment = $cdm_db_schema || '.V_DEID_ENROLLMENT';
set demographic = $cdm_db_schema || '.V_DEID_DEMOGRAPHIC';
set death = $cdm_db_schema || '.V_DEID_DEATH';
set address = $cdm_db_schema || '.V_DEID_ADDRESS_HISTORY';

-----------------------------------------------------------------------------------

/*denominator patient counts: */
-- all medicare
select count(distinct PATID) pat_cnt
from identifier($enrollment)
; -- 20,254,512

/*cohort initialization: medicare enrolled above 65 yo*/
create or replace table PAT_ENR_65UP as
with pat65up_cte as (
      select distinct d.PATID, d. BIRTH_DATE
      from identifier($demographic) d
      join identifier($enrollment) e
      on d.PATID = e.PATID
      where floor(datediff('day',d.BIRTH_DATE,e.ENR_START_DATE)/365.25) >= 65
)
select e.*, p.BIRTH_DATE
from identifier($enrollment) e
join pat65up_cte p on p.PATID = e.PATID
;
select count(distinct patid) from PAT_ENR_65UP;
-- 14,156,442

/*cohort initialization: at least 1 OSA*/
create or replace table PAT_OSA_INIT as
select distinct PATID, DX, DX_TYPE, DX_DATE
from identifier($diagnosis) 
where (DX_TYPE = '09' and DX in ('327.20','327.23','327.29','780.51','780.53','780.57')) or
      (DX_TYPE = '10' and DX in ('G47.30','G47.33','G47.39'))
;
select count(distinct PATID) from PAT_OSA_INIT;
-- 3,295,943

/*cohort inclusion: 2 or more OSA codes at different dates*/
create or replace table PAT_OSA_2DX as
select PATID 
      ,min(DX_DATE) as OSA_DX1_DATE
      ,count(distinct DX_DATE) as DX_DISTINCT_CNT
from PAT_OSA_INIT
group by PATID
having DX_DISTINCT_CNT >= 2
;
select count(distinct PATID) from PAT_OSA_2DX;
-- 2,613,147

/*cohort inclusion: 1 washup period priod to OSA_DX1_DATE*/
create or replace table PAT_OSA_ELIG as
with days_before_osa as (
      select osa.PATID
            ,osa.OSA_DX1_DATE
            ,osa.DX_DISTINCT_CNT
            ,enr.ENR_START_DATE::date as ENR_START_DATE
            ,enr.ENR_END_DATE::date as ENR_END_DATE
            ,floor(datediff('day',enr.BIRTH_DATE,osa.OSA_DX1_DATE)/365.25) as AGE_AT_OSA_DX1
            ,floor(datediff('day',enr.BIRTH_DATE,enr.ENR_START_DATE)/365.25) as AGE_AT_ENR_START
            ,datediff('day',enr.ENR_START_DATE,osa.OSA_DX1_DATE) as DAYS_ENR_START_TO_OSA
            ,datediff('day',enr.ENR_START_DATE,enr.ENR_END_DATE) as DAYS_ENR_START_TO_END
            ,datediff('day',osa.OSA_DX1_DATE::date,enr.ENR_END_DATE::date) as DAYS_OSA_TO_CENSOR
      from PAT_OSA_2DX osa
      join PAT_ENR_65UP enr on enr.patid = osa.patid
)
select * from days_before_osa 
where DAYS_ENR_START_TO_OSA >= 365 -- 1-full year observation window preceding first OSA diagnosis
;
select count(distinct patid) from PAT_OSA_ELIG;
-- 907,759

/*adverse outcome identification, pre-post: MACE, All-cause mortality*/
set pat_elig = 'PAT_OSA_ELIG';
create or replace table PAT_OSA_ENDPT as
with mace_event as (
    select dx.PATID
          ,'MI' ENDPOINT
          ,dx.DX_DATE ENDPOINT_DATE
    from identifier($DIAGNOSIS) dx
    where (dx.DX_TYPE = '09' and split_part(dx.DX,'.',1) in ( '410','412')) OR
          (dx.DX_TYPE = '10' and split_part(dx.DX,'.',1) in ( 'I21','I22')) OR 
          (dx.DX_TYPE = '10' and substr(dx.DX,1,5) in ( 'I25.2'))
    union
    select dx.PATID
          ,'STROKE' ENDPOINT
          ,dx.DX_DATE ENDPOINT_DATE
    from identifier($DIAGNOSIS) dx
    where (dx.DX_TYPE = '09' and split_part(dx.DX,'.',1) in ( '431','434')) OR
          (dx.DX_TYPE = '10' and split_part(dx.DX,'.',1) in ( 'I61','I62','I63','I64'))
    union
    select dx.PATID
          ,'HF' ENDPOINT
          ,dx.DX_DATE ENDPOINT_DATE
    from identifier($DIAGNOSIS) dx
    where (dx.DX_TYPE = '09' and split_part(dx.DX,'.',1) in ( '428')) OR
          (dx.DX_TYPE = '10' and split_part(dx.DX,'.',1) in ( 'I50'))
    union
    select px.PATID
          ,'REVASC' ENDPOINT
          ,px.PX_DATE ENDPOINT_DATE
    from identifier($PROCEDURES) px
    where (px.PX_TYPE = '09' and substr(px.PX,1,4) in ( '36.0','36.1')) OR
          (px.PX_TYPE = '10' and substr(px.PX,1,3) in ( '021','027')) OR
          (px.PX_TYPE = 'CH' and
           px.PX in ( '92920','92921','92924','92925','92928','92929'
                     ,'92933','92934','92937','92938','92941','92943'
                     ,'92944','92980','92981','92982','92984','92995'
                     ,'92996','92973','92974'))
), mace_occur as (
select p.PATID
      ,datediff('day',p.OSA_DX1_DATE,m.ENDPOINT_DATE) as DAYS_SINCE_OSA
      ,case when datediff('day',p.OSA_DX1_DATE,m.ENDPOINT_DATE)<=0 then (m.ENDPOINT || '_BEF')
            else (m.ENDPOINT || '_AFT') end as MACE_BEF_AFT
from identifier($pat_elig) p
left join mace_event m on p.PATID = m.PATID 
), mace_pivot as (
select * from mace_occur
pivot (min(DAYS_SINCE_OSA) for MACE_BEF_AFT 
       in ('MI_BEF','MI_AFT','STROKE_BEF','STROKE_AFT','HF_BEF','HF_AFT','REVASC_BEF','REVASC_AFT')) 
       as p(PATID,MI_BEF,MI_AFT,STROKE_BEF,STROKE_AFT,HF_BEF,HF_AFT,REVASC_BEF,REVASC_AFT)
) 
select p.PATID
      ,p.OSA_DX1_DATE
      ,p.AGE_AT_OSA_DX1
      ,m.MI_BEF
      ,m.MI_AFT
      ,m.STROKE_BEF
      ,m.STROKE_AFT
      ,m.HF_BEF
      ,m.HF_AFT
      ,m.REVASC_BEF
      ,m.REVASC_AFT
      ,NULLIF(least(NVL(m.MI_BEF, 999999),NVL(m.STROKE_BEF, 999999),NVL(m.HF_BEF, 999999),NVL(m.REVASC_BEF, 999999)),999999) as MACE_BEF
      ,NULLIF(least(NVL(m.MI_AFT, 999999),NVL(m.STROKE_AFT, 999999),NVL(m.HF_AFT, 999999),NVL(m.REVASC_AFT, 999999)),999999) as MACE_AFT 
      ,d.DEATH_DATE
      ,datediff('day',p.OSA_DX1_DATE,d.DEATH_DATE) as DAYS_OSA_TO_DEATH
from identifier($pat_elig) p
left join mace_pivot m on p.PATID = m.PATID
left join identifier($DEATH) d on p.PATID = d.PATID
;

select count(distinct patid) from PAT_OSA_ENDPT
where MACE_BEF is null
;
-- 594,984

/*exposure: OSA diagnostics*/
create or replace table PAT_OSA_TEST as
select distinct
       p.PATID
      ,p.OSA_DX1_DATE
      ,px.PX as TEST_CODE
      ,px.PX_DATE as TEST_DATE
      ,datediff('day',p.OSA_DX1_DATE,px.PX_DATE) as DAYS_DIAGNOSTIC_TO_OSA
from identifier($pat_elig) p
join identifier($PROCEDURES) px
on p.PATID = px.PATID
where px.PX in ('95782','95783','95800','95801',
                '95803','95805','95806','95807',
                '95808','95810','95811') and 
      px.PX_TYPE = 'CH'
;
select count(distinct patid) from PAT_OSA_TEST;
-- 404,099

/*exposure: Oral device/appliance*/
create or replace table PAT_OSA_MAD as
select distinct
       p.PATID
      ,p.OSA_DX1_DATE
      ,px.PX as MAD_CODE
      ,px.PX_DATE as MAD_DATE
      ,datediff('day',px.PX_DATE,p.OSA_DX1_DATE) as DAYS_OSA_TO_MAD_INIT
from identifier($pat_elig) p
join identifier($procedures) px
on p.PATID = px.PATID
where px.PX in ('E0485','E0486')
      and px.PX_TYPE = 'CH' 
      -- and px.PX_DATE > p.OSA_DX1_DATE
;
select count(distinct patid) from PAT_OSA_MAD;
-- 7,941

/*exposure: CPAP initialization*/
create or replace table PAT_OSA_PAP_DEVICE as
      with cpap_first_date as (
      select p.PATID
            ,p.OSA_DX1_DATE
            ,min(px.PX_DATE) as CPAP_INIT_DATE
            ,datediff('day',p.OSA_DX1_DATE,min(px.PX_DATE)) as DAYS_OSA_TO_CPAP_INIT
      from identifier($pat_elig) p
      join identifier($procedures) px
      on p.PATID = px.PATID
      where px.PX in ('E0601', 'E0470', 'E0471', '94660') and px.PX_TYPE = 'CH'
      group by p.PATID,p.OSA_DX1_DATE
)
select * from cpap_first_date
where DAYS_OSA_TO_CPAP_INIT > 0
;
select count(distinct patid) from PAT_OSA_PAP_DEVICE;
-- 296,203

/*exposure: CPAP adherence*/
create or replace table PAT_OSA_PAP_SUPPLY as
select d.PATID
      ,e.ENR_START_DATE
      ,e.DAYS_ENR_START_TO_OSA
      ,d.OSA_DX1_DATE
      ,d.CPAP_INIT_DATE
      ,d.DAYS_OSA_TO_CPAP_INIT
      ,listagg(px.PX,'|') within group (order by px.PX) as HCPCS_LIST
      ,max(case when px.PX in ('E0601','E0470','E0471') then 1 else 0 end) as DEVICE_IND
      ,px.PX_DATE as CPAP_SUPPLY_DATE
      ,datediff('day',d.CPAP_INIT_DATE::date,px.PX_DATE::date) as DAYS_CPAP_INIT_TO_SUPPLY
      ,floor(datediff('day',d.CPAP_INIT_DATE::date,px.PX_DATE::date)/365.25) as YEAR_CPAP_INIT_TO_SUPPLY
      ,e.ENR_END_DATE
      ,datediff('day',px.PX_DATE::date,e.ENR_END_DATE::date) as DAYS_SUPPLY_TO_ENR_END
from PAT_OSA_PAP_DEVICE d
join identifier($pat_elig) e on d.PATID = e.PATID
left join identifier($procedures) px on d.PATID = px.PATID
where px.PX_TYPE = 'CH' and px.PX_DATE >= d.CPAP_INIT_DATE and
      px.PX in ('A4604','A7027','A7028','A7029',
                'A7030','A7031','A7032','A7033',
                'A7034','A7035','A7036','A7037',
                'A7038','A7039','A7044','A7046',
                'E0561','E0562',
                'E0601','E0470','E0471')
group by d.PATID,e.ENR_START_DATE,e.DAYS_ENR_START_TO_OSA,d.OSA_DX1_DATE,d.CPAP_INIT_DATE,
         d.DAYS_OSA_TO_CPAP_INIT,px.PX_DATE,e.ENR_END_DATE
order by d.PATID, px.PX_DATE
;

/*covariates: demographics */
create or replace table PAT_OSA_DEMO as
select distinct
       p.PATID
      ,d.BIRTH_DATE
      ,round(datediff('day',d.BIRTH_DATE,p.OSA_DX1_DATE)/365.25) AGE_AT_OSA_DX1
      ,d.SEX
      ,d.RACE
      ,d.HISPANIC
from identifier($pat_elig) p
join identifier($demographic) d 
on p.patid = d.patid
;

/*covariates: diagnoses*/
create or replace table PAT_OSA_ALL_DX as
select distinct
       p.PATID
      ,dx.DX
      ,dx.DX_TYPE
      ,dx.DX_DATE
      ,datediff('day',p.OSA_DX1_DATE,dx.DX_DATE) as DAYS_SINCE_INDEX
from identifier($pat_elig) p
join identifier($diagnosis) dx 
on p.patid = dx.patid
;
select count(distinct patid) from PAT_OSA_ALL_DX;
-- 907,759

/*covariates: selective comorb*/
select * from Z_REF_COMORB;
create or replace table PAT_OSA_SEL_DX as 
with cte_comorb_map as (
      select dx.patid,
             dx.dx_date,
             dx.days_since_index,
             cmb.code_grp,
             cmb.full as code_grp_lbl,
             row_number() over (partition by dx.patid,cmb.code_grp, case when dx.days_since_index <=0 then 1 else 0 end order by dx.days_since_index) as rn
      from PAT_OSA_ALL_DX dx
      join Z_REF_COMORB cmb 
      on dx.dx like cmb.code || '%' and 
         dx.dx_type = cmb.code_type
)
select patid
      ,dx_date as comorb_date
      ,days_since_index
      ,code_grp
      ,code_grp_lbl
from cte_comorb_map 
where rn = 1 
; 
select code_grp_lbl, count(distinct patid) from PAT_OSA_SEL_DX
group by code_grp_lbl;
 

/*covariates: cci*/
select * from Z_REF_CCI;
create or replace table PAT_OSA_CCI as
with cte_cci_map as (
      select dx.patid,
             dx.dx_date,
             dx.days_since_index,
             cci.code_grp,
             cci.full as code_grp_lbl,
             cci.score as cci_score,
             row_number() over (partition by dx.patid,cci.code_grp, case when dx.days_since_index <=0 then 1 else 0 end order by dx.days_since_index) as rn
      from PAT_OSA_ALL_DX dx
      join Z_REF_CCI cci 
      on dx.dx like cci.code || '%' and 
         dx.dx_type = cci.code_type
)
select patid
      ,dx_date as cci_date
      ,days_since_index
      ,code_grp
      ,code_grp_lbl
      ,cci_score
from cte_cci_map 
where rn = 1
; 

/*covariates: medications*/
select * from z_ref_rx_ndc;
create or replace table PAT_OSA_SEL_RX as
with cte_cls as (
      select p.PATID
            ,case when z.VACLS in ('CV350') then 'ALP'
                  when z.VACLS in ('BL110') then 'ACG'
                  when z.VACLS in ('HS500') then 'BGR' 
                  else 'AHT' 
             end as RXCLS
            ,rx.DISPENSE_DATE
            ,datediff('day',p.OSA_DX1_DATE,dx.DISPENSE_DATE) as DAYS_SINCE_OSA
      from identifier($pat_elig) p
      join identifier($dispensing) rx 
      on p.patid = dx.patid
      join z_ref_rx_ndc z 
      on rx.ndc = z.ndc
), cte_ord as (
      select patid
            ,rxcls
            ,dispense_date
            ,days_since_osa
            ,row_number() over (partition by patid,rxcls,case_when days_since_osa<=0 then 1 else 0 end order by abs(days_since_osa)) as rn
      from cte_cls
)
select patid
      ,rxcls
      ,dispense_date
      ,days_since_osa
from cte_ord
where rn = 1
;
