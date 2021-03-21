DELETE FROM @target_database_schema.@target_cohort_table where cohort_definition_id = @target_cohort_id;

insert into @target_database_schema.@target_cohort_table (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
select @target_cohort_id as cohort_definition_id, subject_id, cohort_start_date, cohort_end_date
from
	(
  select t0.subject_id, t2.visit_start_date as cohort_start_date, t2.visit_start_date as cohort_end_date 
  from
  (
    select op1.person_id as subject_id, datefromparts(2017,1,1) as cohort_start_date, datefromparts(2017,1,1) as cohort_end_date
    from @cdm_database_schema.observation_period op1
    where dateadd(dd,365,op1.observation_period_start_date) <= datefromparts(2017,1,1)
    and op1.observation_period_end_date >= datefromparts(2017,1,1)
  ) t0
  inner join
  (
    select person_id, visit_start_date
  	from
  	(
  	   select person_id, row_number() over (partition by person_id order by newid()) as rn1, visit_start_date
  		  from @cdm_database_schema.visit_occurrence
  		  where visit_start_date >= '1/1/2017'
  		  and visit_start_date < '1/1/2018'
  	) t1
  	where rn1 = 1
  ) t2
  on t0.subject_id = t2.person_id  

	union all


  select t0.subject_id, t2.visit_start_date as cohort_start_date, t2.visit_start_date as cohort_end_date 
  from
  (
    select op1.person_id as subject_id, datefromparts(2018,1,1) as cohort_start_date, datefromparts(2018,1,1) as cohort_end_date
    from @cdm_database_schema.observation_period op1
    where dateadd(dd,365,op1.observation_period_start_date) <= datefromparts(2018,1,1)
    and op1.observation_period_end_date >= datefromparts(2018,1,1)
  ) t0
  inner join
  (
    select person_id, visit_start_date
  	from
  	(
  	   select person_id, row_number() over (partition by person_id order by newid()) as rn1, visit_start_date
  		  from @cdm_database_schema.visit_occurrence
  		  where visit_start_date >= '1/1/2018'
  		  and visit_start_date < '1/1/2019'
  	) t1
  	where rn1 = 1
  ) t2
  on t0.subject_id = t2.person_id  


	 union all


  select t0.subject_id, t2.visit_start_date as cohort_start_date, t2.visit_start_date as cohort_end_date 
  from
  (
    select op1.person_id as subject_id, datefromparts(2019,1,1) as cohort_start_date, datefromparts(2019,1,1) as cohort_end_date
    from @cdm_database_schema.observation_period op1
    where dateadd(dd,365,op1.observation_period_start_date) <= datefromparts(2019,1,1)
    and op1.observation_period_end_date >= datefromparts(2019,1,1)
  ) t0
  inner join
  (
    select person_id, visit_start_date
  	from
  	(
  	   select person_id, row_number() over (partition by person_id order by newid()) as rn1, visit_start_date
  		  from @cdm_database_schema.visit_occurrence
  		  where visit_start_date >= '1/1/2019'
  		  and visit_start_date < '1/1/2020'
  	) t1
  	where rn1 = 1
  ) t2
  on t0.subject_id = t2.person_id  
) t3
;
