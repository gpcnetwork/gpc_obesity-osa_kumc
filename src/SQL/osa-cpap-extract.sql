/*
# Copyright (c) 2021-2025 University of Missouri                   
# Author: Xing Song, xsm7f@umsystem.edu                            
# File: osa-cpap-extract.sql
# Description: extracting osa population from Medicare claims and identify cpap exposure and study endpoints, 
#              as well as a few pre-defined covariates
*/

use role GROUSE_ROLE_C_ANALYTICS; 
use warehouse GROUSE_WH;
use database GROUSE_DEID_ANALYTICS_DB;

create schema if not exists OBESITY_OSA;
use schema OBESITY_OSA;

set cdm_db_schema = 'GROUSE_DEID_DB.CMS_PCORNET_CDM';
set diagnosis = $cdm_db_schema || '.V_DEID_DIAGNOSIS';
set procedures = $cdm_db_schema || '.V_DEID_PROCEDURES';
set enrollment = $cdm_db_schema || '.V_DEID_ENROLLMENT';
set demographic = $cdm_db_schema || '.V_DEID_DEMOGRAPHIC';
set death = $cdm_db_schema || '.V_DEID_DEATH';

-----------------------------------------------------------------------------------

/*denominator patient counts: */
-- all medicare
select count(distinct PATID) pat_cnt
from identifier($enrollment)
; -- 13,706,965

-- all medicare above 65 yo
select count(distinct d.PATID) as pat_cnt
from identifier($demographic) d
join identifier($enrollment) e
on d.PATID = e.PATID
where floor(datediff('day',d.BIRTH_DATE,e.ENR_START_DATE)/365.25) >= 65
; -- 9,273,239

/*cohort inclusion: 5-consecutive year of enrollment*/
create or replace table PAT_ENR5YR as
with enr_cte as (
    select enr.PATID
          ,demo.BIRTH_DATE
          ,min(enr.ENR_START_DATE::date) as ENR_START_DATE
          ,max(enr.ENR_END_DATE::date) as ENR_END_DATE
          ,floor(datediff('day',demo.BIRTH_DATE,min(enr.ENR_START_DATE))/365.25) as AGE_AT_ENR_START
          ,datediff('day',min(enr.ENR_START_DATE),max(enr.ENR_END_DATE)) as DAYS_ENR_START_TO_END
    from identifier($enrollment) enr
    join identifier($demographic) demo on enr.PATID = demo.PATID
    where enr.RAW_BASIS like 'AB%' 
      -- or enr.RAW_BASIS like 'MC%' 
    group by enr.PATID,demo.BIRTH_DATE
)
select * from enr_cte where DAYS_ENR_START_TO_END >= 365.25*5
;
select count(distinct patid) from PAT_ENR5YR;
-- 6,772,162

/*cohort initialization: at least 1 OSA*/
create or replace table PAT_OSA_INIT as
select distinct PATID, DX, DX_TYPE, DX_DATE
from identifier($diagnosis) 
where (DX_TYPE = '09' and DX in ('327.20','327.23','327.29','780.51','780.53','780.57')) or
      (DX_TYPE = '10' and DX in ('G47.30','G47.33','G47.39'))
;
select count(distinct PATID) from PAT_OSA_INIT;
-- 1,795,214

/*cohort inclusion: 2 or more OSA codes at different dates*/
create or replace table PAT_OSA_2DX as
select PATID 
      ,min(DX_DATE) as OSA_DX1_DATE
      ,count(distinct DX_DATE) as DX_DISTINCT_CNT
from PAT_OSA_INIT
group by PATID
having DX_DISTINCT_CNT >= 2
;
select count(distinct PATID), count(*) from PAT_OSA_2DX;
-- 1,368,901

/*cohort inclusion: 5-consecutive year of enrollment*/
create or replace table PAT_OSA_INCLD as
select osa.PATID
      ,osa.PATID_DEID
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
join PAT_ENR5YR enr on enr.patid = osa.patid
;
select count(distinct patid) from PAT_OSA_INCLD;
-- 658,642

/*eligible cohort*/
create or replace table PAT_OSA_ELIG as
select PATID
      ,OSA_DX1_DATE
      ,AGE_AT_OSA_DX1
      ,AGE_AT_ENR_START
      ,ENR_START_DATE
      ,ENR_END_DATE
      ,DAYS_ENR_START_TO_OSA
      ,DAYS_OSA_TO_CENSOR
from PAT_OSA_INCLD
where DAYS_ENR_START_TO_OSA >= 365 -- 1-full year observation window preceding first OSA diagnosis
      and AGE_AT_ENR_START >= 65   -- most likely to be Medicare-eligible due to age
;
select count(distinct patid) from PAT_OSA_ELIG;
-- 225,132

/*summarize consort counts*/ 
select '>= 1 OSA diagnosis code' as EVENT, count(distinct patid) as PAT_CNT from PAT_OSA_INIT
union
select '>= 2 OSA diagnosis code' as EVENT, count(distinct patid) as PAT_CNT from PAT_OSA_2DX
union
select '5-year conscutive' as EVENT, count(distinct patid) as PAT_CNT from PAT_OSA_INCLD
union
select '5-year conscutive, age at enrollment>= 65' as EVENT, count(distinct patid) as PAT_CNT 
from PAT_OSA_INCLD where AGE_AT_ENR_START >= 65 
union
select '5-year conscutive, age at enrollment>= 65, 1-full-year enrollment before OSA' as EVENT, count(distinct patid) as PAT_CNT from PAT_OSA_ELIG
order by PAT_CNT desc
;

/*exposure: OSA diagnostics*/
-- create or replace table PAT_OSA_TEST as
-- select distinct
--        p.PATID
--       ,p.OSA_DX1_DATE
--       ,px.PX as TEST_CODE
--       ,px.PX_DATE as TEST_DATE
--       ,datediff('day',p.OSA_DX1_DATE,px.PX_DATE) as DAYS_DIAGNOSTIC_TO_OSA
-- from PAT_OSA_ELIG p
-- join identifier($PROCEDURES) px
-- on p.PATID = px.PATID
-- where px.PX in ('95782','95783','95800','95801',
--                 '95803','95805','95806','95807',
--                 '95808','95810','95811') and 
--       px.PX_TYPE = 'CH'
-- ;

/*exposure: Oral device/appliance*/
create or replace table PAT_MAD as
select distinct
       p.PATID
      ,p.OSA_DX1_DATE
      ,px.PX as MAD_CODE
      ,px.PX_DATE as MAD_DATE
      ,datediff('day',px.PX_DATE,p.OSA_DX1_DATE) as DAYS_OSA_TO_MAD_INIT
from PAT_OSA_ELIG p
join identifier($procedures) px
on p.PATID = px.PATID
where px.PX in ('E0485','E0486')
      and px.PX_TYPE = 'CH' 
      and px.PX_DATE > p.OSA_DX1_DATE
;

select count(distinct patid) from PAT_MAD;
-- 2,132

/*exposure: CPAP initialization*/
create or replace table PAT_OSA_CPAP_DEVICE as
with cpap_first_date as (
select p.PATID
      ,p.OSA_DX1_DATE
      ,min(px.PX_DATE) as CPAP_INIT_DATE
      ,datediff('day',p.OSA_DX1_DATE,min(px.PX_DATE)) as DAYS_OSA_TO_CPAP_INIT
from PAT_OSA_ELIG p
join identifier($procedures) px
on p.PATID = px.PATID
where px.PX in ('E0601', 'E0470', 'E0471', '94660') and px.PX_TYPE = 'CH'
group by p.PATID,p.OSA_DX1_DATE
)
select * from cpap_first_date
where DAYS_OSA_TO_CPAP_INIT > 0
;

select count(distinct patid) from PAT_OSA_CPAP_DEVICE;
-- 112,429

/*exposure: CPAP adherence*/
create or replace table PAT_OSA_CPAP_SUPPLY as
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
from PAT_OSA_CPAP_DEVICE d
join PAT_OSA_ELIG e on d.PATID = e.PATID
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

-- create or replace table PAT_OSA_CPAP_SUPPLY_YEARLY as
-- select distinct 
--        PATID, YEAR_CPAP_INIT_TO_SUPPLY,
--        max(YEAR_CPAP_INIT_TO_SUPPLY) over (partition by PATID) as MAX_SUPPLY_YEAR,
--        count(distinct CPAP_SUPPLY_DATE) over (partition by PATID,YEAR_CPAP_INIT_TO_SUPPLY) as YEARLY_SUPPLY_COUNT
-- from PAT_OSA_CPAP_SUPPLY
-- order by PATID, YEAR_CPAP_INIT_TO_SUPPLY
-- ;

/*outcome identification: MACE, All-cause mortality*/
create or replace table PAT_OSA_ENDPT as
with mace_event as (
    select dx.PATID
          ,'MI' ENDPOINT
          ,dx.DX_DATE ENDPOINT_DATE
    from identifier($DIAGNOSIS) dx
    where (dx.DX_TYPE = '09' and split_part(dx.DX,'.',1) in ( '410','412')) OR
          (dx.DX_TYPE = '10' and split_part(dx.DX,'.',1) in ( 'I21','I22','I23'))
    union
    select dx.PATID
          ,'STROKE' ENDPOINT
          ,dx.DX_DATE ENDPOINT_DATE
    from identifier($DIAGNOSIS) dx
    where (dx.DX_TYPE = '09' and split_part(dx.DX,'.',1) in ( '431','434')) OR
          (dx.DX_TYPE = '10' and split_part(dx.DX,'.',1) in ( 'I61','I62','I63'))
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
from PAT_OSA_ELIG p
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
from PAT_OSA_ELIG p
left join mace_pivot m on p.PATID = m.PATID
left join identifier($DEATH) d on p.PATID = d.PATID
;

-- with denom as (
-- select count(distinct patid) as N from PAT_OSA_ELIG
-- )
-- select 'Eligible' as EVENT, count(distinct patid) as PAT_CNT, round(count(distinct patid)/denom.N,3) as PAT_PROP 
-- from PAT_OSA_ELIG cross join denom group by denom.N
-- union
-- select 'MACE_recurrent', count(distinct patid), round(count(distinct patid)/denom.N,3) 
-- from PAT_OSA_ENDPT cross join denom where MACE_BEF is not null and MACE_AFT is not null group by denom.N
-- union
-- select 'MACE_new', count(distinct patid), round(count(distinct patid)/denom.N,3) 
-- from PAT_OSA_ENDPT cross join denom where MACE_BEF is null and MACE_AFT is not null group by denom.N
-- union
-- select 'Death', count(distinct patid), round(count(distinct patid)/denom.N,3) 
-- from PAT_OSA_ENDPT cross join denom where DEATH_DATE is not null group by denom.N
-- union 
-- select 'CPAP initialization', count(distinct patid), round(count(distinct patid)/denom.N,3) as PAT_CNT 
-- from PAT_OSA_CPAP_DEVICE cross join denom group by denom.N
-- order by PAT_CNT desc
-- ;

/*
Eligible	212,445	1
CPAP initialization	107,994	0.508
MACE_recurrent	66,111	0.311
MACE_new	40,236	0.189
Death	23,815	0.112
*/


/*covariates: demographics, hypertension, obesity, T2DM, ...*/
-- Obesity: ICD9 (278.00, 278.01, 278.03) or ICD10 (E66.0,E66.01, E66.09, E66.1, E66.2, E66.8, E66.9, Z68.30-Z68.45)
-- T2DM: ICD9 (250,357.2,362.0[1-7]) or ICD10 (E10, E11, E08.42, E13.42)
-- hypertension: ICD9 (401,402,403,404,405) or ICD10 (I10,I11,I12,I13,R03)
-- Smoking status

create or replace table PAT_OSA_COVARIATES as
with cov_event as (
    select dx.PATID
          ,'Obesity' ENDPOINT
          ,dx.DX_DATE ENDPOINT_DATE
    from identifier($DIAGNOSIS) dx
    where (dx.DX_TYPE = '09' and dx.DX in ('278.00', '278.01', '278.03')) OR
          (dx.DX_TYPE = '10' and (split_part(dx.DX,'.',1) in ('E66') or substr(dx.DX,1,5) in ('Z68.3','Z68.4')))
    union
    select dx.PATID
          ,'T2DM' ENDPOINT
          ,dx.DX_DATE ENDPOINT_DATE
    from identifier($DIAGNOSIS) dx
    where (dx.DX_TYPE = '09' and (split_part(dx.DX,'.',1) in ('250') or substr(dx.DX,1,5) in ('357.2','362.0'))) OR
          (dx.DX_TYPE = '10' and (split_part(dx.DX,'.',1) in ('E10','E11') or substr(dx.DX,1,6) in ('E08.42', 'E13.42')))
    union
    select dx.PATID
          ,'HTN' ENDPOINT
          ,dx.DX_DATE ENDPOINT_DATE
    from identifier($DIAGNOSIS) dx
    where (dx.DX_TYPE = '09' and split_part(dx.DX,'.',1) in ('401','402','403','404','405')) OR
          (dx.DX_TYPE = '10' and split_part(dx.DX,'.',1) in ('I10','I11','I12','I13','R03'))
), cov_occur as (
select p.PATID
      ,datediff('day',p.OSA_DX1_DATE,m.ENDPOINT_DATE) as DAYS_SINCE_OSA
      ,case when datediff('day',p.OSA_DX1_DATE,m.ENDPOINT_DATE)<=0 then (m.ENDPOINT || '_BEF')
            else (m.ENDPOINT || '_AFT') end as MACE_BEF_AFT
from PAT_OSA_ELIG p
left join cov_event m on p.PATID = m.PATID 
), cov_pivot as (
select * from cov_occur
pivot (min(DAYS_SINCE_OSA) for MACE_BEF_AFT 
       in ('Obesity_BEF','Obesity_AFT','T2DM_BEF','T2DM_AFT','HTN_BEF','HTN_AFT')) 
       as p(PATID,Obesity_BEF,Obesity_AFT,T2DM_BEF,T2DM_AFT,HTN_BEF,HTN_AFT)
), cov_pivot2 as (
select * from cov_occur
pivot (max(DAYS_SINCE_OSA) for MACE_BEF_AFT 
       in ('Obesity_BEF','Obesity_AFT','T2DM_BEF','T2DM_AFT','HTN_BEF','HTN_AFT')) 
       as p(PATID,Obesity_BEF,Obesity_AFT,T2DM_BEF,T2DM_AFT,HTN_BEF,HTN_AFT)
), cov_ses as (
select distinct PATID, 1 as LIS_DUAL_IND
from identifier($enrollment)
where raw_basis like 'LIS%' or raw_basis like 'DUAL%'
)
select p.PATID
      ,p.OSA_DX1_DATE
      ,p.AGE_AT_OSA_DX1
      ,d.BIRTH_DATE
      ,d.SEX
      ,d.RACE
      ,d.HISPANIC
      ,coalesce(cov_ses.LIS_DUAL_IND,0) as LIS_DUAL_IND
      ,cov.Obesity_BEF as Obesity_BEF_min
      ,cov2.Obesity_BEF as Obesity_BEF_max
      ,cov.Obesity_AFT as Obesity_AFT_min
      ,cov2.Obesity_AFT as Obesity_AFT_max
      ,cov.T2DM_BEF as T2DM_BEF_min
      ,cov2.T2DM_BEF as T2DM_BEF_max
      ,cov.T2DM_AFT as T2DM_AFT_min
      ,cov2.T2DM_AFT as T2DM_AFT_max
      ,cov.HTN_BEF as HTN_BEF_min
      ,cov2.HTN_BEF as HTN_BEF_max
      ,cov.HTN_AFT as HTN_AFT_min
      ,cov2.HTN_AFT as HTN_AFT_max
from PAT_OSA_ELIG p
left join cov_ses on p.PATID = cov_ses.PATID
left join cov_pivot cov on p.PATID = cov.PATID
left join cov_pivot2 cov2 on p.PATID = cov2.PATID
left join identifier($demographic) d on p.PATID = d.PATID
;

/*de-identified dataset for Aim 1: association between CPAP usage and health outcomes*/
create or replace table PAT_OSA_CPAP_DEID as
select p.PATID
      ,p.DAYS_OSA_TO_CENSOR
      ,d.DAYS_OSA_TO_CPAP_INIT
      ,ep.MI_AFT AS DAYS_OSA_TO_MI
      ,ep.STROKE_AFT AS DAYS_OSA_TO_STROKE
      ,ep.HF_AFT AS DAYS_OSA_TO_HF
      ,ep.REVASC_AFT AS DAYS_OSA_TO_REVASC
      ,ep.MACE_AFT AS DAYS_OSA_TO_MACE
      ,ep.DAYS_OSA_TO_DEATH
      ,cov.AGE_AT_OSA_DX1
      ,cov.SEX
      ,cov.RACE
      ,cov.HISPANIC
      ,cov.LIS_DUAL_IND
      ,case when ep.MI_BEF < 0 then 1 else 0 end as MI_history
      ,case when ep.STROKE_BEF < 0 then 1 else 0 end as STROKE_history
      ,case when ep.HF_BEF < 0 then 1 else 0 end as HF_history
      ,case when ep.REVASC_BEF < 0 then 1 else 0 end as REVASC_history
      ,case when ep.MACE_BEF < 0 then 1 else 0 end as MACE_history
      ,case when cov.Obesity_BEF_min < 0 then 1 else 0 end as Obesity_history --obesity dx before OSA dx1
      ,case when cov.T2DM_BEF_min < 0 then 1 else 0 end as T2DM_history --T2DM dx before OSA dx1
      ,case when cov.HTN_BEF_min < 0 then 1 else 0 end as HTN_history --HTN dx before OSA dx1
from PAT_OSA_ELIG p
left join PAT_OSA_CPAP_DEVICE d on p.PATID = d.PATID
left join PAT_OSA_ENDPT ep on p.PATID = ep.PATID
left join PAT_OSA_COVARIATES cov on p.PATID = cov.PATID
;

/*de-id dataset for Aim 2: association between CPAP adherence and health outcomes*/
create or replace table PAT_OSA_CPAP_SUPPLY_DEID as 
select distinct
       p.PATID
      ,s.DAYS_ENR_START_TO_OSA
      ,s.DAYS_OSA_TO_CPAP_INIT
      ,s.DAYS_CPAP_INIT_TO_SUPPLY
      ,s.DAYS_SUPPLY_TO_ENR_END + s.DAYS_CPAP_INIT_TO_SUPPLY as DAYS_CPAP_INIT_TO_CENSOR
      ,cov.AGE_AT_OSA_DX1
      ,cov.SEX
      ,cov.RACE
      ,cov.HISPANIC
      ,cov.LIS_DUAL_IND
      ,ep.MI_AFT  - s.DAYS_OSA_TO_CPAP_INIT as DAYS_CPAP_INIT_TO_MI
      ,ep.STROKE_AFT - s.DAYS_OSA_TO_CPAP_INIT as DAYS_CPAP_INIT_TO_STROKE
      ,ep.HF_AFT - s.DAYS_OSA_TO_CPAP_INIT as DAYS_CPAP_INIT_TO_HF
      ,ep.REVASC_AFT - s.DAYS_OSA_TO_CPAP_INIT as DAYS_CPAP_INIT_TO_REVASC
      ,ep.MACE_AFT - s.DAYS_OSA_TO_CPAP_INIT as DAYS_CPAP_INIT_TO_MACE
      ,ep.DAYS_OSA_TO_DEATH - s.DAYS_OSA_TO_CPAP_INIT as DAYS_CPAP_INIT_TO_DEATH
      ,case when ep.MI_BEF < 0 then 1 else 0 end as MI_history
      ,case when ep.STROKE_BEF < 0 then 1 else 0 end as STROKE_history
      ,case when ep.HF_BEF < 0 then 1 else 0 end as HF_history
      ,case when ep.REVASC_BEF < 0 then 1 else 0 end as REVASC_history
      ,case when ep.MACE_BEF < 0 then 1 else 0 end as MACE_history
      ,case when cov.Obesity_BEF_min < 0 then 1 else 0 end as Obesity_history --obesity dx before OSA dx1
      ,case when cov.T2DM_BEF_min < 0 then 1 else 0 end as T2DM_history --T2DM dx before OSA dx1
      ,case when cov.HTN_BEF_min < 0 then 1 else 0 end as HTN_history --HTN dx before OSA dx1
      ,cov.Obesity_AFT_min
      ,cov.Obesity_AFT_max
      ,cov.T2DM_AFT_min
      ,cov.T2DM_AFT_max
      ,cov.HTN_AFT_min
      ,cov.HTN_AFT_max
from PAT_OSA_ELIG p
join PAT_OSA_CPAP_SUPPLY s on p.PATID = s.PATID -- only CPAP user group
left join PAT_OSA_ENDPT ep on p.PATID = ep.PATID
left join PAT_OSA_COVARIATES cov on p.PATID = cov.PATID
order by PATID,DAYS_CPAP_INIT_TO_SUPPLY
;
