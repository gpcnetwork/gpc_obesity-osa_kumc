
/*
review the SAVE trial inclusion criteria: 
- patients age between 45 and 75 years old (might be hard to define in observational cohort, we will not restrict age for now)
- with a diagnosis of coronary artery disease or cerebrovascular disease
- and a diagnosis of obstructive sleep apnea
*/

-- identify coronary artery disease or cerebrovascular disease
create or replace table CAD_CVD_DX_ALL as
select *
from deidentified_pcornet_cdm.cdm_c015r031.deid_diagnosis
where
    -- coronary artery disease
    (
     dx like '410%' or dx like '411%' or dx like '412%' or dx like '413%' or dx like '414%' or 
     dx like 'I20%' or dx like 'I21%' or dx like 'I22%' or dx like 'I23%' or dx like 'I24%' or dx like 'I25%'
    ) or 
    -- cerebrovascular disease
    (
     dx like '430%' or dx like '431%' or dx like '432%' or dx like '433%' or dx like '434%' or dx like '435%' or dx like '436%' or dx like '437%' or dx like '438%' or
     dx like 'I60%' or dx like 'I61%' or dx like 'I62%' or dx like 'I63%' or dx like 'I64%' or dx like 'I65%' or dx like 'I66%' or dx like 'I67%' or dx like 'I68%' or dx like 'I69%' 
    ) and 
    dx_date <= '2024-01-31' -- data release date
;
select count(distinct patid) from CAD_CVD_DX_ALL;
-- 125,439

-- identify OSA 
create or replace table OSA_DX_ALL as
select *
from deidentified_pcornet_cdm.cdm_c015r031.deid_diagnosis
where (DX_TYPE = '09' and DX in ('327.20','327.23','327.29','780.51','780.53','780.57')) or
      (DX_TYPE = '10' and DX in ('G47.30','G47.33','G47.39')) and 
      dx_date <= '2024-01-31' -- data release date
;
select count(distinct patid) from OSA_DX_ALL;
-- 77,857

-- collect initial OSA diagnosis date
create or replace table OSA_1DX as 
select patid, 
       min(dx_date) as OSA_1DX_DATE
from OSA_DX_ALL
group by patid
;

-- identify cohort with 1st OSA succeed a CAD or CVD event
create or replace table CAD_CVD_OSA as 
select distinct
       a.patid, 
       min(b.dx_date::date) as CAD_CVD_1DX_DATE, 
       a.OSA_1DX_DATE::date as OSA_1DX_DATE -- remove the time chunk
from OSA_1DX a 
join CAD_CVD_DX_ALL b 
on a.patid = b.patid 
where b.dx_date <= a.OSA_1DX_DATE
group by a.patid,a.OSA_1DX_DATE::date
;
select count(distinct patid), count(*) from CAD_CVD_OSA;
-- 16,768

select * from CAD_CVD_OSA limit 5;

-- identify cohort with 1st OSA succeed a CAD or CVD event, who are also undergoing usual care
select * from MACE_OSA_CPAP.MEDICATION limit 5;

create or replace table CAD_CVD_OSA_USUAL_ALL as 
select distinct
       a.*, 
       m.cls,
       m.subcls,
       m.ing,
       p.rx_order_date
from CAD_CVD_OSA a 
join deidentified_pcornet_cdm.cdm_c015r031.deid_prescribing p 
on a.patid = p.patid 
join MACE_OSA_CPAP.MEDICATION m 
on p.rxnorm_cui = m.rxcui
-- this time window is to mimic the intervention assignment around the time of OSA diagnosis
where datediff('day',p.rx_order_date,a.OSA_1DX_DATE) between -90 and 90
;
select count(distinct patid) from CAD_CVD_OSA_USUAL_ALL;
-- 11,735

-- summarize meds dates and counts of distinct drug classes used for usual care
-- 1 patient per row
create or replace table CAD_CVD_OSA_USUAL_PER_CLS as 
select patid, 
       case when MED_ANTIHTN is not null then 1 else 0 end as MED_ANTIHTN_IND,
       case when MED_ANTIDIABETIC is not null then 1 else 0 end as MED_ANTIDIABETIC_IND,
       case when MED_ASPIRIN_ANTITHROM is not null then 1 else 0 end as MED_ASPIRIN_ANTITHROM_IND,
       case when MED_STATI_LIPID is not null then 1 else 0 end as MED_STATI_LIPID_IND
from (
    select distinct patid, cls, rx_order_date
    from CAD_CVD_OSA_USUAL_ALL
)
pivot (
    min(rx_order_date) for 
    cls in ('antihypertensive','antidiabetic','antirhtombotic','antilipidemic')
) 
as p(patid,MED_ANTIHTN,MED_ANTIDIABETIC,MED_ASPIRIN_ANTITHROM,MED_STATI_LIPID)
;
select count(distinct patid), count(*) from CAD_CVD_OSA_USUAL_PER_CLS;
select * from CAD_CVD_OSA_USUAL_PER_CLS limit 5;


create or replace table CAD_CVD_OSA_USUAL as 
select a.PATID, 
       d.BIRTH_DATE,
       d.SEX,
       d.RACE,
       d.HISPANIC,
       a.CAD_CVD_1DX_DATE,
       datediff(year,d.birth_date,a.CAD_CVD_1DX_DATE) as AGE_AT_CAD_CVD,
       a.OSA_1DX_DATE,
       datediff(year,d.birth_date,a.OSA_1DX_DATE) as AGE_AT_OSA,
       per.MED_ANTIHTN_IND,
       per.MED_ANTIDIABETIC_IND,
       per.MED_ASPIRIN_ANTITHROM_IND,
       per.MED_STATI_LIPID_IND,
       count(distinct b.cls) as MED_CLS_COUNT,
       count(distinct b.subcls) as MED_SUBCLS_COUNT,
       min(b.rx_order_date::date) as MED_DATE1,
       max(b.rx_order_date::date) as MED_DATE2
from CAD_CVD_OSA_USUAL_ALL a 
join CAD_CVD_OSA_USUAL_PER_CLS per
on a.patid = per.patid
join CAD_CVD_OSA_USUAL_ALL b
on a.patid = b.patid
join deidentified_pcornet_cdm.cdm_c015r031.deid_demographic d
on a.patid = d.patid
where datediff(year,d.birth_date,a.OSA_1DX_DATE) >= 18
group by a.PATID, a.CAD_CVD_1DX_DATE,a.OSA_1DX_DATE,
         per.MED_ANTIHTN_IND,per.MED_ANTIDIABETIC_IND,per.MED_ASPIRIN_ANTITHROM_IND,per.MED_STATI_LIPID_IND,
         d.birth_date,d.sex,d.race,d.hispanic
;
select count(distinct patid), count(*) from CAD_CVD_OSA_USUAL;
-- 11,735 -> 11,700

select * from CAD_CVD_OSA_USUAL limit 5;

-- CAD_CVD_OSA_USUAL is our eligible cohort
-- identify PAP initialization
create or replace table CAD_CVD_OSA_PAP_DATE1 as 
select a.patid, 
       min(px.admit_date::date) as PAP_DATE1,
       datediff('day',a.OSA_1DX_DATE,min(px.admit_date::date)) as DAYS_OSA_PAP
from CAD_CVD_OSA_USUAL a 
join deidentified_pcornet_cdm.cdm_c015r031.deid_procedures px
where px.px in ('E0601','E0470','E0471','94660'/*,'95811'*/) and 
      px.admit_date between a.CAD_CVD_1DX_DATE and '2024-01-31'
group by a.patid,a.OSA_1DX_DATE
;
select count(distinct patid), count(*) from CAD_CVD_OSA_PAP_DATE1;
-- 11,660 -> 11625

select ceil(DAYS_OSA_PAP/365) as YEAR_APART,
       count(distinct patid) as pat_cnt, 
       round(count(distinct patid)/11660, 3) as pat_prop
from CAD_CVD_OSA_PAP_DATE1
group by ceil(DAYS_OSA_PAP/365)
order by YEAR_APART;
-- it seems that almost everyone in our eligible cohort has being billed for PAP at some point
-- > 44% had ther first PAP date beyond 1 year before first OSA 


-- identify the treated group: usual-care + PAP
-- modify the definition of treated group as:
--     those who have been using PAP between -1 and +1 year since first OSA diagnosis
create or replace table CAD_CVD_OSA_USUAL_PAP as 
select a.*, pap.PAP_DATE1
from CAD_CVD_OSA_USUAL a 
join CAD_CVD_OSA_PAP_DATE1 pap
on a.patid = pap.patid
where datediff('day',a.OSA_1DX_DATE,pap.PAP_DATE1) between -90 and 90
;
select count(distinct patid), count(*) from CAD_CVD_OSA_USUAL_PAP;
-- 3,597 -> 3,589

-- identify the control group: usual-care only
-- modify the definition of control group as:
--     those who haven't used PAP between -1 and +1 year since first OSA diagnosis
create or replace table CAD_CVD_OSA_USUAL_NOPAP as 
select a.*, pap.PAP_DATE1
from CAD_CVD_OSA_USUAL a 
left join CAD_CVD_OSA_PAP_DATE1 pap
on a.patid = pap.patid
where datediff('day',a.OSA_1DX_DATE,pap.PAP_DATE1) < -90 or 
      datediff('day',a.OSA_1DX_DATE,pap.PAP_DATE1) > 90 or 
      pap.PAP_DATE1 is null
;
select count(distinct patid), count(*) from CAD_CVD_OSA_USUAL_NOPAP;
-- 8,138 -> 8,111

/*
identify MACE event sequence, which will be used for: 
1. reassign the treated into control if MACE occurred before PAP
2. identify MACE outcomes
Ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8571870/
*/

create or replace table MACE_DX_SEQ as 
-- myocardial infarction
select distinct patid, dx_date::date as MACE_DATE, 'MI' as MACE_TYPE
from deidentified_pcornet_cdm.cdm_c015r031.deid_diagnosis
where dx like '410%' or dx like '411%' or 
      dx like 'I21%' or dx like 'I22%' or dx like 'I23%'   
union
-- hospitalization for heart failure
select distinct patid, dx_date::date as MACE_DATE, 'HF_HOSP'
from deidentified_pcornet_cdm.cdm_c015r031.deid_diagnosis
where dx like '428%' or dx like 'I50%' and enc_type in ('EI','IP')
union
-- stroke
select distinct patid, dx_date::date as MACE_DATE, 'STROKE'
from deidentified_pcornet_cdm.cdm_c015r031.deid_diagnosis
where dx like '433%' or dx like '434%' or dx like '436%' or 
      dx like 'I63%' or dx like 'I65%' or dx like 'I66%'
union
-- revascularization procedure
select distinct patid, px_date::date as MACE_DATE,'REVASC'
from deidentified_pcornet_cdm.cdm_c015r031.deid_procedures
where (PX_TYPE = '09' and substr(PX,1,4) in ( '36.0','36.1')) OR
      (PX_TYPE = '10' and substr(PX,1,3) in ( '021','027')) OR
      (PX_TYPE = 'CH' and
      (PX in ( '92920','92921','92924','92925','92928','92929'
              ,'92933','92934','92937','92938','92941','92943'
              ,'92944','92980','92981','92982','92984','92995'
              ,'92996','92973','92974'))) 
union
-- all-cause death
select patid, min(death_date::date) as MACE_DATE, 'DEATH'
from deidentified_pcornet_cdm.cdm_c015r031.deid_death
where death_date is not null
group by patid
;
select mace_type, count(distinct patid)
from MACE_DX_SEQ
group by mace_type
;

select * from MACE_DX_SEQ order by patid limit 5;

-- from the treated group, there might be someone who developed MACE before PAP was initiated
-- they could be either excluded or reassigned to control group
create or replace table CAD_CVD_OSA_USUAL_PAP_REASSIGN as 
select distinct a.*
from CAD_CVD_OSA_USUAL_PAP a
join MACE_DX_SEQ b 
on a.patid = b.patid 
where a.cad_cvd_1dx_date < b.mace_date and
      b.mace_date <= a.pap_date1
;
select count(distinct patid), count(*) from CAD_CVD_OSA_USUAL_PAP_REASSIGN;
-- 1,020

/*
update the treated and control group and finalize the eligible cohort
also, create the exposure indicator PAP_IND
*/

create or replace table PAT_ELIG as 
select patid,
       birth_date,
       sex,
       race,
       hispanic,
       cad_cvd_1dx_date,
       age_at_cad_cvd,
       osa_1dx_date,
       age_at_osa,
       med_cls_count,
       med_subcls_count,
       med_antidiabetic_ind,
       med_antihtn_ind,
       med_aspirin_antithrom_ind,
       med_stati_lipid_ind,
       pap_date1,
       1 as PAP_IND
from CAD_CVD_OSA_USUAL_PAP
where patid not in (select patid from CAD_CVD_OSA_USUAL_PAP_REASSIGN) -- exclude reassigned cases
union 
select patid,
       birth_date,
       sex,
       race,
       hispanic,
       cad_cvd_1dx_date,
       age_at_cad_cvd,
       osa_1dx_date,
       age_at_osa,
       med_cls_count,
       med_subcls_count,
       med_antidiabetic_ind,
       med_antihtn_ind,
       med_aspirin_antithrom_ind,
       med_stati_lipid_ind,
       pap_date1,
       0 as PAP_IND
from CAD_CVD_OSA_USUAL_NOPAP
union 
select patid,
       birth_date,
       sex,
       race,
       hispanic,
       cad_cvd_1dx_date,
       age_at_cad_cvd,
       osa_1dx_date,
       age_at_osa,
       med_cls_count,
       med_subcls_count,
       med_antidiabetic_ind,
       med_antihtn_ind,
       med_aspirin_antithrom_ind,
       med_stati_lipid_ind,
       pap_date1,
       0 as PAP_IND
from CAD_CVD_OSA_USUAL_PAP_REASSIGN
;
select PAP_IND, count(distinct patid), count(*) from PAT_ELIG
group by PAP_IND;
-- 1	2569
-- 0	9131

create or replace table MACE_ANY_POS as 
select a.patid, a.pap_ind, 
       min(MACE_DATE) as MACE_ANY_DATE
from PAT_ELIG a 
join MACE_DX_SEQ m 
on a.patid = m.patid 
where m.mace_date > a.cad_cvd_1dx_date and 
      (m.mace_date > a.pap_date1 or a.pap_date1 is null)
group by a.patid, a.pap_ind
;
select pap_ind, count(distinct patid), count(*) from MACE_ANY_POS
group by pap_ind;
-- 1	1246
-- 0	5315

create or replace table PAT_CENSOR as 
select a.patid, a.pap_ind, a.osa_1dx_date,
       max(e.admit_date::date) as censor_date
from PAT_ELIG a 
join deidentified_pcornet_cdm.cdm_c015r031.deid_encounter e
on a.patid = e.patid
group by a.patid, a.pap_ind, a.osa_1dx_date
;
select PAP_IND, count(distinct patid), count(*) from PAT_CENSOR
group by PAP_IND;
-- 1	2569
-- 0	9131

create or replace table MACE_OUTCOME as 
select cs.patid, cs.pap_ind,
       case when m.MACE_ANY_DATE is not null then 1 
            else 0 
       end as mace_recur_ind,
       coalesce(m.MACE_ANY_DATE,cs.censor_date) as endpoint_date,
       datediff(day,cs.osa_1dx_date,coalesce(m.MACE_ANY_DATE,cs.censor_date)) as days_osa_to_endpoint_date
from PAT_CENSOR cs 
left join MACE_ANY_POS m 
on cs.patid = m.patid
;

select * from MACE_OUTCOME limit 5;

select PAP_IND, mace_recur_ind, count(distinct patid), count(*) 
from MACE_OUTCOME
group by PAP_IND,mace_recur_ind;

/* 
collect other baseline information
baseline is defined as before first OSA diagnosis
*/

-- baseline BMI
create or replace table BASE_BMI as 
with ht as (
    select a.patid, median(v.ht) as ht, 
    from pat_elig a 
    join deidentified_pcornet_cdm.cdm_c015r031.deid_vital v 
    on a.patid = v.patid
    where datediff(year,a.birth_date,v.measure_date) >= 18
    group by a.patid
),  wt as (
    select patid, wt
    from (
        select a.patid, v.wt, 
               row_number() over (partition by a.patid order by abs(datediff(day,a.osa_1dx_date,v.measure_date::date))) as rn 
        from pat_elig a 
        join deidentified_pcornet_cdm.cdm_c015r031.deid_vital v 
        on a.patid = v.patid 
    )
    where rn = 1
), bmi as (
    select patid, original_bmi as bmi
    from (
        select a.patid, v.original_bmi, 
               row_number() over (partition by a.patid order by abs(datediff(day,a.osa_1dx_date,v.measure_date::date))) as rn 
        from pat_elig a 
        join deidentified_pcornet_cdm.cdm_c015r031.deid_vital v 
        on a.patid = v.patid 
        where v.original_bmi is not null and v.original_bmi between 10 and 200
    )
    where rn = 1
), bmi_calc as (
    select a.patid, 
           round(ht.ht) as ht, 
           round(wt.wt) as wt, 
           round(wt.wt/(ht.ht*ht.ht)*703) as bmi_calc
    from pat_elig a
    join ht on a.patid = ht.patid
    join wt on a.patid = wt.patid
    where ht.ht > 0 and wt.wt is not null and round(wt.wt/(ht.ht*ht.ht)*703) between 10 and 200
)
select a.patid, 
       bmi_calc.ht,
       bmi_calc.wt,
       coalesce(round(bmi.bmi),bmi_calc.bmi_calc) as bmi
from pat_elig a 
left join bmi_calc on a.patid = bmi_calc.patid
left join bmi on a.patid = bmi.patid
where coalesce(round(bmi.bmi),bmi_calc.bmi_calc) is not null
;

select count(distinct patid), count(*) from BASE_BMI;
-- 11,322

select * from base_bmi limit 5;


-- baseline SMOKING STATUS
CREATE OR REPLACE TABLE SMOKING_HIS AS
with smoking_valid as (
    SELECT a.patid,
           a.pap_ind,
           b.smoking,
           row_number() over (partition by a.patid order by abs(datediff('day',b.measure_date,a.osa_1dx_date))) as rn
    FROM PAT_ELIG a
    JOIN deidentified_pcornet_cdm.cdm_c015r031.DEID_VITAL b
    ON a.patid = b.patid
    where b.smoking not in ('UN','OT','NI') and b.smoking is not null and 
          datediff('day',a.osa_1dx_date,b.measure_date) <= 90
)
select patid, pap_ind,
       case when smoking in ('01','02','03','05','07','08') then 'SM' 
            when smoking in ('04') then 'NSM'
       else 'NI' end as smoking_status 
from smoking_valid
where rn = 1
; 

select smoking_status, count(distinct patid)
from SMOKING_HIS
group by smoking_status
;

-- baseline indicator of CAD 
create or replace table CAD_HIS as 
select pe.patid, pe.pap_ind, 
       min(dx_date::date) as CAD_HIS_DATE, 1 as CAD_HIS
from pat_elig pe
join CAD_CVD_DX_ALL a
on pe.patid = a.patid
where (
    a.dx like '410%' or a.dx like '411%' or a.dx like '412%' or a.dx like '413%' or a.dx like '414%' or 
    a.dx like 'I20%' or a.dx like 'I21%' or a.dx like 'I22%' or a.dx like 'I23%' or a.dx like 'I24%' or a.dx like 'I25%'
) and a.dx_date <= pe.osa_1dx_date
group by pe.patid, pe.pap_ind
;
select pap_ind, count(distinct patid), count(*)
from CAD_HIS
group by pap_ind
;
-- 1	2013
-- 0	7014

-- baseline indicator of CVD 
create or replace table CVD_HIS as 
select pe.patid, pe.pap_ind,
       min(dx_date::date) as CVD_HIS_DATE, 1 as CVD_HIS
from pat_elig pe  
join CAD_CVD_DX_ALL a
on pe.patid = a.patid
where  (
    a.dx like '430%' or a.dx like '431%' or a.dx like '432%' or a.dx like '433%' or a.dx like '434%' or a.dx like '435%' or a.dx like '436%' or a.dx like '437%' or a.dx like '438%' or
    a.dx like 'I60%' or a.dx like 'I61%' or a.dx like 'I62%' or a.dx like 'I63%' or a.dx like 'I64%' or a.dx like 'I65%' or a.dx like 'I66%' or a.dx like 'I67%' or a.dx like 'I68%' or a.dx like 'I69%' 
) and a.dx_date <= pe.osa_1dx_date
group by pe.patid, pe.pap_ind
;
select pap_ind, count(distinct patid), count(*)
from CVD_HIS
group by pap_ind
;
-- 1	696
-- 0	4242

-- baseline indicator of IHD
create or replace table IHD_HIS as 
select pe.patid, pe.pap_ind, 
       min(dx_date::date) as IHD_HIS_DATE, 1 as IHD_HIS
from pat_elig pe  
join CAD_CVD_DX_ALL a
on pe.patid = a.patid
where  (
    a.dx like '411%' or a.dx like '412%' or a.dx like '413%' or a.dx like '414%' or
    a.dx like 'I21%' or a.dx like 'I22%' or a.dx like 'I23%' or a.dx like 'I24%' 
) and a.dx_date <= pe.osa_1dx_date
group by pe.patid, pe.pap_ind
;
select pap_ind, count(distinct patid), count(*)
from IHD_HIS
group by pap_ind
;
-- 1	936
-- 0	5212

-- baseline indicator of any stroke
create or replace table STROKE_HIS as 
select pe.patid, pe.pap_ind, 
       min(mace_date) as STROKE_HIS_DATE, 1 as STROKE_HIS
from pat_elig pe  
join MACE_DX_SEQ a
on pe.patid = a.patid
where mace_type = 'STROKE' and 
      a.mace_date <= pe.osa_1dx_date
group by pe.patid, pe.pap_ind
;
select pap_ind, count(distinct patid), count(*)
from STROKE_HIS
group by pap_ind
;
-- 1	430
-- 0	3266

-- baseline indicator of myocardial infarction
create or replace table MI_HIS as 
select pe.patid, pe.pap_ind, 
       min(mace_date) as MI_HIS_DATE, 1 as MI_HIS
from pat_elig pe  
join MACE_DX_SEQ a
on pe.patid = a.patid
where mace_type = 'MI' and 
      a.mace_date <= pe.osa_1dx_date
group by pe.patid, pe.pap_ind
;
select pap_ind, count(distinct patid), count(*)
from MI_HIS
group by pap_ind
;
-- 1	251
-- 0	2416

-- baseline indicator of any heart failure
create or replace table HF_HIS as 
select pe.patid, pe.pap_ind, 
       min(mace_date) as HF_HIS_DATE, 1 as HF_HIS
from pat_elig pe  
join MACE_DX_SEQ a
on pe.patid = a.patid
where mace_type = 'HF_HOSP' and 
      a.mace_date <= pe.osa_1dx_date
group by pe.patid, pe.pap_ind
;
select pap_ind, count(distinct patid), count(*)
from HF_HIS
group by pap_ind
;
-- 1	581
-- 0	3174

-- baseline indicator of revascularization survery
create or replace table REVASC_HIS as 
select pe.patid, pe.pap_ind, 
       min(mace_date) as REVASC_HIS_DATE, 1 as REVASC_HIS
from pat_elig pe  
join MACE_DX_SEQ a
on pe.patid = a.patid
where mace_type = 'REVASC' and 
      a.mace_date <= pe.osa_1dx_date
group by pe.patid, pe.pap_ind
;
select pap_ind, count(distinct patid), count(*)
from REVASC_HIS
group by pap_ind
;
-- 1	80
-- 0	1322

-- baseline indicator of hypertention
CREATE OR REPLACE TABLE HTN AS
SELECT pe.PATID, pe.pap_ind,
       MIN(b.DX_DATE::date) AS HTN_HIS_DATE, 1 as HTN_HIS
FROM pat_elig pe 
JOIN deidentified_pcornet_cdm.cdm_c015r031.DEID_DIAGNOSIS b
ON pe.patid = b.patid
WHERE b.dx LIKE 'I10%' or b.dx LIKE 'I11%' or b.dx LIKE 'I12%' or b.dx LIKE 'I13%' or 
      b.dx LIKE '401%' or b.dx LIKE '402%' or b.dx LIKE '403%' or b.dx LIKE '404%' or b.dx LIKE '405%'
group by pe.PATID, pe.pap_ind
;  
select pap_ind, count(distinct patid), count(*)
from HTN
group by pap_ind
;
-- 1	2368
-- 0	8636

-- baseline indicator for diabetes  
CREATE OR REPLACE TABLE DM2 AS
SELECT pe.PATID, pe.pap_ind,
       MIN(b.DX_DATE::date) AS DM_HIS_DATE, 1 as DM_HIS
FROM pat_elig pe 
JOIN deidentified_pcornet_cdm.cdm_c015r031.DEID_DIAGNOSIS b
ON pe.patid = b.patid
WHERE b.dx LIKE 'E08%' or b.dx LIKE 'E09%' or b.dx LIKE 'E10%' or b.dx LIKE 'E11%' or b.dx LIKE 'E12%' or b.dx LIKE 'E13%' or
      b.dx LIKE '250%' or b.dx LIKE '357.2%' or b.dx LIKE '362.0%'
group by pe.PATID, pe.pap_ind
;
select pap_ind, count(distinct patid), count(*)
from DM2
group by pap_ind
;
-- 1	1567
-- 0	5956

-- baseline CCI 
create or replace table CCI_HIS as 
with dx_all as (
    select a.patid, 
           c.code_grp,
           c.score,
           row_number() over (partition by a.patid,c.code_grp order by datediff(day,coalesce(b.dx_date,b.admit_date),a.osa_1dx_date)) as rn
    from pat_elig a
    join deidentified_pcornet_cdm.cdm_c015r031.deid_diagnosis b
    on a.patid = b.patid and coalesce(b.dx_date,b.admit_date)<= a.osa_1dx_date
    join shared_db.depression.cci_ref c 
    on b.dx = c.code and b.dx_type = c.code_type
), dx_cci_tot as (
    select patid, sum(score) as cci_tot
    from dx_all 
    where rn = 1
    group by patid
)
select a.patid, 
       tot.cci_tot,
       case when tot.cci_tot between 1 and 2 then 'CCI_GRP1'
            when tot.cci_tot between 3 and 4 then 'CCI_GRP2'
            when tot.cci_tot >= 5 then 'CCI_GRP3'
       end as cci_class
from pat_elig a 
join dx_cci_tot tot on a.patid = tot.patid
;

select cci_class, count(distinct patid), count(*) 
from CCI_HIS
group by cci_class;
-- CCI_GRP1	5975
-- CCI_GRP2	1841
-- CCI_GRP3	578

/*
final analytic dataset
*/
create or replace table MACE_OSA_CPAP as 
select pe.patid,
       pe.pap_ind,
       oc.mace_recur_ind,
       oc.endpoint_date,
       oc.days_osa_to_endpoint_date,
       pe.sex,
       case when pe.race = '05' then 'WH'
            when pe.race = '03' then 'AA'
            when pe.race = '02' then 'AS'
            when pe.race in ('01','04','05','07','OT') then 'OT'
            else 'NI'
       end as race, 
       case when pe.hispanic in ('Y','N') then pe.hispanic
            else 'NI'
       end as hispanic,
       pe.AGE_AT_CAD_CVD,
       pe.AGE_AT_OSA,
       pe.MED_CLS_COUNT,
       pe.MED_SUBCLS_COUNT,
       pe.med_antihtn_ind,
       pe.med_antidiabetic_ind,
       pe.med_aspirin_antithrom_ind,
       pe.med_stati_lipid_ind,
       -- pe.cad_cvd_1dx_date,
       -- pe.osa_1dx_date.
       -- pe.pap_date1,
       -- oc.mace_any_date,
       bmi.bmi,         
       case when bmi.bmi < 18.5 then 'underweight'
            when bmi.bmi between 18.5 and 24.9 then 'normal'
            when bmi.bmi between 25.0 and 29.9 then 'overweight'
            when bmi.bmi between 30.0 and 39.9 then 'obesity'
            when bmi.bmi >= 40 then 'severe_obesity'
       end as bmi_class, 
       coalesce(sm.smoking_status,'NI') as smoking_status,
       coalesce(cad.cad_his,0) as cad_his,
       coalesce(cvd.cvd_his,0) as cvd_his,
       coalesce(hf.hf_his,0) as hf_his,
       coalesce(ihd.ihd_his,0) as ihd_his,
       coalesce(mi.mi_his,0) as mi_his,
       coalesce(re.revasc_his,0) as revasc_his,
       coalesce(st.stroke_his,0) as stroke_his,
       coalesce(htn.htn_his,0) as htn_his,
       coalesce(dm2.dm_his,0) as dm_his,
       coalesce(cci.cci_tot,0) as cci_tot,
       coalesce(cci.cci_class,'CCI_GRP0') as cci_class
from PAT_ELIG pe 
join MACE_OUTCOME oc on pe.patid = oc.patid
join base_bmi bmi on pe.patid = bmi.patid
left join SMOKING_HIS sm on pe.patid = sm.patid
left join CAD_HIS cad on pe.patid = cad.patid
left join CVD_HIS cvd on pe.patid = cvd.patid
left join HF_HIS hf on pe.patid = hf.patid
left join IHD_HIS ihd on pe.patid = ihd.patid
left join MI_HIS mi on pe.patid = mi.patid
left join REVASC_HIS re on pe.patid = re.patid
left join STROKE_HIS st on pe.patid = st.patid 
left join HTN on pe.patid = htn.patid
left join DM2 on pe.patid = dm2.patid
left join CCI_HIS cci on pe.patid = cci.patid
;

select count(distinct patid), count(*) from MACE_OSA_CPAP;
-- 11322


