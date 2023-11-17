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
set dispensing = $cdm_db_schema || '.V_DEID_DISPENSING';
set pat_elig = 'PAT_OSA_ELIG';


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
            ,max(enr.ENR_END_DATE::date) over (partition by osa.PATID) as CENSOR_DATE
            ,max(datediff('day',osa.OSA_DX1_DATE::date,enr.ENR_END_DATE::date)) over (partition by osa.PATID) as DAYS_OSA_TO_CENSOR
      from PAT_OSA_2DX osa
      join PAT_ENR_65UP enr on enr.patid = osa.patid
)
select distinct
       PATID
      ,OSA_DX1_DATE
      ,DX_DISTINCT_CNT
      ,AGE_AT_OSA_DX1
      ,CENSOR_DATE
      ,DAYS_OSA_TO_CENSOR
from days_before_osa 
where DAYS_ENR_START_TO_OSA >= 365 and -- 1-full year observation window preceding first OSA diagnosis
      AGE_AT_OSA_DX1 >= 65 and -- above 65yo at initial OSA
      DAYS_OSA_TO_CENSOR > 0 -- pre-mature censor
;

select count(distinct patid), count(*) from PAT_OSA_ELIG;
-- 888,845

/*adverse outcome identification, pre-post: MACE, All-cause mortality*/
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
            px.PX in (
                   '92920','92921','92924','92925','92928','92929'
                  ,'92933','92934','92937','92938','92941','92943'
                  ,'92944','92980','92981','92982','92984','92995'
                  ,'92996','92973','92974'
            ))    
)
select p.PATID
      ,e.ENDPOINT
      ,e.ENDPOINT_DATE
      ,datediff('day',p.OSA_DX1_DATE,e.ENDPOINT_DATE) as DAYS_OSA_TO_ENDPOINT
      ,d.DEATH_DATE
      ,datediff('day',p.OSA_DX1_DATE,d.DEATH_DATE) as DAYS_OSA_TO_DEATH
      ,c.CENSOR_DATE
      ,c.DAYS_OSA_TO_CENSOR
from identifier($pat_elig) p
inner join PAT_OSA_ELIG c on p.PATID = c.PATID
left join mace_event e on p.PATID = e.PATID 
left join identifier($death) d on p.PATID = d.PATID
where DAYS_OSA_TO_CENSOR >= 0
;

select count(distinct patid), count(*) from PAT_OSA_ENDPT;

/*exposure: OSA diagnostics*/
-- create or replace table PAT_OSA_TEST as
-- select distinct
--        p.PATID
--       ,p.OSA_DX1_DATE
--       ,px.PX as TEST_CODE
--       ,px.PX_DATE as TEST_DATE
--       ,datediff('day',p.OSA_DX1_DATE,px.PX_DATE) as DAYS_DIAGNOSTIC_TO_OSA
-- from identifier($pat_elig) p
-- join identifier($PROCEDURES) px
-- on p.PATID = px.PATID
-- where px.PX in (
--       '95782','95783','95800','95801',
--       '95803','95805','95806','95807',
--       '95808','95810','95811'
--       ) and 
--       px.PX_TYPE = 'CH'
-- ;
-- select count(distinct patid) from PAT_OSA_TEST;
-- 396,019

/*exposure: Oral device/appliance*/
-- create or replace table PAT_OSA_MAD as
-- select distinct
--        p.PATID
--       ,p.OSA_DX1_DATE
--       ,px.PX as MAD_CODE
--       ,px.PX_DATE as MAD_DATE
--       ,datediff('day',px.PX_DATE,p.OSA_DX1_DATE) as DAYS_OSA_TO_MAD_INIT
-- from identifier($pat_elig) p
-- join identifier($procedures) px
-- on p.PATID = px.PATID
-- where px.PX in ('E0485','E0486')
--       and px.PX_TYPE = 'CH' 
--       -- and px.PX_DATE > p.OSA_DX1_DATE
-- ;
-- select count(distinct patid) from PAT_OSA_MAD;
-- 7,853

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
select count(distinct patid), count(*)
 from PAT_OSA_PAP_DEVICE;
-- 290,015

/*exposure: CPAP adherence*/
create or replace table PAT_OSA_PAP_SUPPLY as
select d.PATID
      ,d.OSA_DX1_DATE
      ,d.CPAP_INIT_DATE
      ,d.DAYS_OSA_TO_CPAP_INIT
      ,listagg(px.PX,'|') within group (order by px.PX) as HCPCS_LIST
      ,max(case when px.PX in ('E0601','E0470','E0471','94660') then 1 else 0 end) as DEVICE_IND
      ,px.PX_DATE as CPAP_SUPPLY_DATE
      ,datediff('day',d.CPAP_INIT_DATE::date,px.PX_DATE::date) as DAYS_CPAP_INIT_TO_SUPPLY
      ,floor(datediff('day',d.CPAP_INIT_DATE::date,px.PX_DATE::date)/365.25) as YEAR_CPAP_INIT_TO_SUPPLY
      ,e.CENSOR_DATE
      ,datediff('day',px.PX_DATE::date,e.CENSOR_DATE) as DAYS_SUPPLY_TO_CENSOR
from PAT_OSA_PAP_DEVICE d
join identifier($pat_elig) e on d.PATID = e.PATID
left join identifier($procedures) px on d.PATID = px.PATID
where px.PX_TYPE = 'CH' and px.PX_DATE >= d.CPAP_INIT_DATE and
      px.PX in 
      (
            'A4604','A7027','A7028','A7029',
            'A7030','A7031','A7032','A7033',
            'A7034','A7035','A7036','A7037',
            'A7038','A7039','A7044','A7046',
            'E0561','E0562',
            'E0601','E0470','E0471',
            '94660'
      )
group by d.PATID,d.OSA_DX1_DATE,d.CPAP_INIT_DATE,
         d.DAYS_OSA_TO_CPAP_INIT,px.PX_DATE,e.CENSOR_DATE
order by d.PATID, px.PX_DATE
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
-- 888,845

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
select code_grp,code_grp_lbl, count(distinct patid) 
from PAT_OSA_SEL_DX
group by code_grp,code_grp_lbl;
 

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
select code_grp,code_grp_lbl, count(distinct patid) 
from PAT_OSA_CCI
group by code_grp,code_grp_lbl;

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
            ,datediff('day',p.OSA_DX1_DATE,rx.DISPENSE_DATE) as DAYS_SINCE_INDEX
      from identifier($pat_elig) p
      join identifier($dispensing) rx on p.patid = rx.patid
      join z_ref_rx_ndc z on rx.ndc = z.ndc
), cte_ord as (
      select patid
            ,rxcls
            ,dispense_date
            ,days_since_index
            ,row_number() over (partition by patid,rxcls,case when days_since_index<=0 then 1 else 0 end order by abs(days_since_index)) as rn
      from cte_cls
)
select patid
      ,rxcls
      ,dispense_date
      ,days_since_index
from cte_ord
where rn = 1
;

select rxcls, count(distinct patid) 
from PAT_OSA_SEL_RX
group by rxcls
;


-- final analytical dataset
create or replace table PAT_OSA_PAP_EXPOS as 
with cte_cci as (
      select patid, sum(cci_score) as CCI_SCORE
      from PAT_OSA_CCI
      where days_since_index <=0 
      group by patid
), cte_sel_dx as (
      select * from
      (
            select patid, code_grp, 1 as ind from PAT_OSA_SEL_DX
            where days_since_index <= 0 and 
                  code_grp in (
                        'obes',
                        't2dm',
                        'ckd',
                        'copd',
                        'afib',
                        'htn',
                        'nd',
                        'hysm',
                        'insm'
                        )
      )
      pivot 
      (
            max(ind) for code_grp 
            in  ('obes','t2dm','ckd','copd','afib','htn','nd','hysm','insm') 
      )
      as p(patid,obes,t2dm,ckd,copd,afib,htn,nd,hysm,insm)
), cte_sel_rx as (
      select * from
      (
            select patid, rxcls, 1 as ind from PAT_OSA_SEL_RX
            where days_since_index <= 0 and days_since_index >= -365
      )
      pivot 
      (
            max(ind) for rxcls 
            in  ('ACG','AHT','ALP','BGR') 
      )
      as p(patid,ACG,AHT,ALP,BGR)
), cte_sdoh as (
      select distinct patid, 1 as lis_dual_ind 
      from identifier($enrollment)
      where raw_basis in ('LIS','DUAL')
), cte_mace_hist as (
      select * from
      (
            select patid, endpoint, days_osa_to_endpoint 
            from PAT_OSA_ENDPT
            where days_osa_to_endpoint <= 0 
      )
      pivot (min(days_osa_to_endpoint) for ENDPOINT in ('MI','HF','STROKE','REVASC'))
      as p(PATID, MI,HF,STROKE,REVASC)
), cte_mace_endpt as (
      select * from
      (
            select patid, death_date, days_osa_to_death, endpoint, days_osa_to_endpoint 
            from PAT_OSA_ENDPT
            where days_osa_to_endpoint > 0 
      )
      pivot (min(days_osa_to_endpoint) for ENDPOINT in ('MI','HF','STROKE','REVASC'))
      as p(PATID, DEATH_DATE, DAYS_OSA_TO_DEATH, MI,HF,STROKE,REVASC)
)
select distinct
       p.PATID
      ,p.OSA_DX1_DATE
      ,d.BIRTH_DATE
      ,round(datediff('day',d.BIRTH_DATE,p.OSA_DX1_DATE)/365.25) AGE_AT_OSA_DX1
      ,d.SEX
      ,d.RACE
      ,d.HISPANIC
      ,mh.MI as MI_HIST
      ,ep.MI as DAYS_OSA_TO_MI
      ,mh.STROKE as STROKE_HIST 
      ,ep.STROKE as DAYS_OSA_TO_STROKE
      ,mh.HF as HF_HIST
      ,ep.HF as DAYS_OSA_TO_HF
      ,mh.REVASC as REVASC_HIST
      ,ep.REVASC as DAYS_OSA_TO_REVASC
      ,NULLIF(least(NVL(mh.MI, 999999),NVL(mh.STROKE, 999999),NVL(mh.HF, 999999),NVL(mh.REVASC, 999999)),999999) as MACE_HIST
      ,NULLIF(least(NVL(ep.MI, 999999),NVL(ep.STROKE, 999999),NVL(ep.HF, 999999),NVL(ep.REVASC, 999999)),999999) as DAYS_OSA_TO_MACE 
      ,ep.DEATH_DATE
      ,ep.DAYS_OSA_TO_DEATH
      ,case when ep.DEATH_DATE is not null then 1 else 0 end as DEATH_IND
      ,cs.CENSOR_DATE
      ,cs.DAYS_OSA_TO_CENSOR
      ,dv.CPAP_INIT_DATE
      ,dv.DAYS_OSA_TO_CPAP_INIT
      ,case when dv.DAYS_OSA_TO_CPAP_INIT is not null then 1 else 0 end as CPAP_IND
      ,cci.cci_score
      ,dx.obes
      ,dx.t2dm
      ,dx.ckd
      ,dx.copd
      ,dx.afib
      ,dx.htn
      ,dx.nd
      ,dx.hysm
      ,dx.insm
      ,rx.ACG
      ,rx.AHT
      ,rx.ALP
      ,rx.BGR
      ,sdh.LIS_DUAL_IND
from identifier($pat_elig) p
inner join identifier($demographic) d on p.patid = d.patid
inner join PAT_OSA_ENDPT cs  on p.patid = cs.patid
left join cte_mace_endpt ep on p.patid = ep.patid
left join cte_mace_hist mh on p.patid = mh.patid 
left join PAT_OSA_PAP_DEVICE dv on p.patid = dv.patid
left join cte_cci cci on cci.patid = p.patid
left join cte_sel_dx dx on dx.patid = p.patid
left join cte_sel_rx rx on rx.patid = p.patid
left join cte_sdoh sdh on sdh.patid = p.patid
where ep.DAYS_OSA_TO_DEATH > 0 or ep.DAYS_OSA_TO_DEATH is null
;

select count(distinct patid), count(*) from PAT_OSA_PAP_EXPOS;
-- 888,835


create or replace table PAT_OSA_PAP_ADHRN as 
with cte_cci as (
      select patid, sum(cci_score) as CCI_SCORE
      from PAT_OSA_CCI
      where days_since_index <= 365 
      group by patid
), cte_sel_dx as (
      select * from
      (
            select patid, code_grp, 1 as ind from PAT_OSA_SEL_DX
            where days_since_index <= 365 and 
                  code_grp in (
                        'obes',
                        't2dm',
                        'ckd',
                        'copd',
                        'afib',
                        'htn',
                        'nd',
                        'hysm',
                        'insm'
                        )
      )
      pivot 
      (
            max(ind) for code_grp 
            in  ('obes','t2dm','ckd','copd','afib','htn','nd','hysm','insm') 
      )
      as p(patid,obes,t2dm,ckd,copd,afib,htn,nd,hysm,insm)
), cte_sel_rx as (
      select * from
      (
            select patid, rxcls, 1 as ind from PAT_OSA_SEL_RX
            where days_since_index <= 365 and days_since_index >= 0
      )
      pivot 
      (
            max(ind) for rxcls 
            in  ('ACG','AHT','ALP','BGR') 
      )
      as p(patid,ACG,AHT,ALP,BGR)
), cte_sdoh as (
      select distinct patid, 1 as lis_dual_ind 
      from identifier($enrollment)
      where raw_basis in ('LIS','DUAL')
), cte_supply_yr1 as (
      select patid, count(distinct CPAP_SUPPLY_DATE) as cpap_yr1 
      from PAT_OSA_PAP_SUPPLY
      where DAYS_CPAP_INIT_TO_SUPPLY <= 365
      group by patid
), cte_supply_t as (
      select patid, max(DAYS_CPAP_INIT_TO_SUPPLY+DAYS_SUPPLY_TO_CENSOR) as DAYS_CPAP_INIT_TO_SUPPLY_CENSOR 
      from PAT_OSA_PAP_SUPPLY
      group by patid
), cte_mace_hist as (
      select * from
      (
            select a.patid, a.endpoint, a.days_osa_to_endpoint 
            from PAT_OSA_ENDPT a
            where exists (
                  select 1 from PAT_OSA_PAP_DEVICE dv 
                  where a.patid = dv.patid and 
                        a.days_osa_to_endpoint <= dv.days_osa_to_cpap_init + 365
            ) 
      )
      pivot (min(days_osa_to_endpoint) for ENDPOINT in ('MI','HF','STROKE','REVASC'))
      as p(PATID, MI,HF,STROKE,REVASC)
), cte_mace_endpt as (
      select * from
      (
            select a.patid, a.death_date, a.days_osa_to_death, a.censor_date,days_osa_to_censor, a.endpoint, a.days_osa_to_endpoint 
            from PAT_OSA_ENDPT a
            where exists (
                  select 1 from PAT_OSA_PAP_DEVICE dv 
                  where a.patid = dv.patid and 
                        a.days_osa_to_endpoint > dv.days_osa_to_cpap_init + 365
            )  
      )
      pivot (min(days_osa_to_endpoint) for ENDPOINT in ('MI','HF','STROKE','REVASC'))
      as p(PATID, DEATH_DATE, DAYS_OSA_TO_DEATH, CENSOR_DATE, DAYS_OSA_TO_CENSOR, MI,HF,STROKE,REVASC)
)
select distinct
       p.PATID
      ,p.OSA_DX1_DATE
      ,d.BIRTH_DATE
      ,round(datediff('day',d.BIRTH_DATE,p.OSA_DX1_DATE)/365.25) AGE_AT_OSA_DX1
      ,d.SEX
      ,d.RACE
      ,d.HISPANIC
      ,mh.MI as MI_HIST
      ,ep.MI as DAYS_OSA_TO_MI
      ,mh.STROKE as STROKE_HIST 
      ,ep.STROKE as DAYS_OSA_TO_STROKE
      ,mh.HF as HF_HIST
      ,ep.HF as DAYS_OSA_TO_HF
      ,mh.REVASC as REVASC_HIST
      ,ep.REVASC as DAYS_OSA_TO_REVASC
      ,NULLIF(least(NVL(mh.MI, 999999),NVL(mh.STROKE, 999999),NVL(mh.HF, 999999),NVL(mh.REVASC, 999999)),999999) as MACE_HIST
      ,NULLIF(least(NVL(ep.MI, 999999),NVL(ep.STROKE, 999999),NVL(ep.HF, 999999),NVL(ep.REVASC, 999999)),999999) as DAYS_OSA_TO_MACE 
      ,ep.DEATH_DATE
      ,ep.DAYS_OSA_TO_DEATH
      ,cs.CENSOR_DATE
      ,cs.DAYS_OSA_TO_CENSOR
      ,dv.CPAP_INIT_DATE
      ,dv.DAYS_OSA_TO_CPAP_INIT
      ,syr1.cpap_yr1
      ,send.DAYS_CPAP_INIT_TO_SUPPLY_CENSOR
      ,cci.cci_score
      ,dx.obes
      ,dx.t2dm
      ,dx.ckd
      ,dx.copd
      ,dx.afib
      ,dx.htn
      ,dx.nd
      ,dx.hysm
      ,dx.insm
      ,rx.ACG
      ,rx.AHT
      ,rx.ALP
      ,rx.BGR
      ,sdh.LIS_DUAL_IND
from identifier($pat_elig) p
inner join identifier($demographic) d on p.patid = d.patid
inner join PAT_OSA_PAP_DEVICE dv on p.patid = dv.patid
inner join PAT_OSA_ENDPT cs  on p.patid = cs.patid
left join cte_mace_endpt ep on p.patid = ep.patid
left join cte_mace_hist mh on p.patid = mh.patid 
left join cte_supply_yr1 syr1 on p.patid = syr1.patid
left join cte_supply_t send on p.patid = send.patid
left join cte_cci cci on cci.patid = p.patid
left join cte_sel_dx dx on dx.patid = p.patid
left join cte_sel_rx rx on rx.patid = p.patid
left join cte_sdoh sdh on sdh.patid = p.patid
where DAYS_CPAP_INIT_TO_SUPPLY_CENSOR >= 365
;

select count(distinct patid), count(*) from PAT_OSA_PAP_ADHRN;
-- 238,541

