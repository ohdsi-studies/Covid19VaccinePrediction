CREATE TABLE #flu_vaccines (
  concept_id bigint NOT NULL
)
;

INSERT INTO #flu_vaccines (concept_id)
select c.descendant_concept_id as concept_id
  from @cdm_database_schema.CONCEPT_ancestor c
  where ancestor_concept_id in (40213145,42903442,40213150,40213159,40225028,40213156,40213151,40213327,40213148,40213158,36878713,42873961,40225038,40213146,40213143,36879025,40213157,45776076,40213149,40213147,40213152,42903441,40213141,40213153,40213144,40213142,40213155,40164828)
;


DELETE FROM @target_database_schema.@target_cohort_table where cohort_definition_id = @target_cohort_id;

insert into @target_database_schema.@target_cohort_table (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
select @target_cohort_id as cohort_definition_id, subject_id, cohort_start_date, cohort_end_date
from
	(
  select t0.subject_id, t1.vaccine_date as cohort_start_date, t1.vaccine_date as cohort_end_date 
  from
  (
    select op1.person_id as subject_id, datefromparts(2017,1,1) as cohort_start_date, datefromparts(2017,1,1) as cohort_end_date
    from @cdm_database_schema.observation_period op1
    where dateadd(dd,365,op1.observation_period_start_date) <= datefromparts(2017,1,1)
    and op1.observation_period_end_date >= datefromparts(2017,1,1)
  ) t0
  inner join
  (
  	   select de1.person_id, min(drug_exposure_start_date) as vaccine_date
  		  from @cdm_database_schema.drug_exposure de1
  		  inner join #flu_vaccines wp1
  		  on de1.drug_concept_id = wp1.concept_id 
  		  where drug_exposure_start_date >= '1/1/2017'
  		  and drug_exposure_start_date < '1/1/2018'
  		  group by de1.person_id
  ) t1
  on t0.subject_id = t1.person_id  

	union all


  select t0.subject_id, t1.vaccine_date as cohort_start_date, t1.vaccine_date as cohort_end_date 
  from
  (
    select op1.person_id as subject_id, datefromparts(2018,1,1) as cohort_start_date, datefromparts(2018,1,1) as cohort_end_date
    from @cdm_database_schema.observation_period op1
    where dateadd(dd,365,op1.observation_period_start_date) <= datefromparts(2018,1,1)
    and op1.observation_period_end_date >= datefromparts(2018,1,1)
  ) t0
  inner join
  (
  	   select de1.person_id, min(drug_exposure_start_date) as vaccine_date
  		  from @cdm_database_schema.drug_exposure de1
  		  inner join #flu_vaccines wp1
  		  on de1.drug_concept_id = wp1.concept_id 
  		  where drug_exposure_start_date >= '1/1/2018'
  		  and drug_exposure_start_date < '1/1/2019'
  		  group by de1.person_id
  ) t1
  on t0.subject_id = t1.person_id   


	 union all


  select t0.subject_id, t1.vaccine_date as cohort_start_date, t1.vaccine_date as cohort_end_date 
  from
  (
    select op1.person_id as subject_id, datefromparts(2019,1,1) as cohort_start_date, datefromparts(2019,1,1) as cohort_end_date
    from @cdm_database_schema.observation_period op1
    where dateadd(dd,365,op1.observation_period_start_date) <= datefromparts(2019,1,1)
    and op1.observation_period_end_date >= datefromparts(2019,1,1)
  ) t0
  inner join
  (
  	   select de1.person_id, min(drug_exposure_start_date) as vaccine_date
  		  from @cdm_database_schema.drug_exposure de1
  		  inner join #flu_vaccines wp1
  		  on de1.drug_concept_id = wp1.concept_id 
  		  where drug_exposure_start_date >= '1/1/2019'
  		  and drug_exposure_start_date < '1/1/2020'
  		  group by de1.person_id
  ) t1
  on t0.subject_id = t1.person_id   
) t3
;

TRUNCATE TABLE #flu_vaccines;
DROP TABLE #flu_vaccines;
