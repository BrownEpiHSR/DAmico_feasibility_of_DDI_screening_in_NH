/*
OVERALL PROJECT: Feasibility of applying pharmacoepidemiologic drug-drug interaction screening 
				 methods to a population of nursing home residents: An application to clopidogrel

				DESCRIPTION: 
				Screen potential precipitant drugs for DDIs with clopidogrel resulting in major bleed and fall related injury (FRI)
				in a sample of NH residents using the self controlled case series (SCCS) design. 


PROGRAM: conditional_poisson_models_NHDDI_screening_github

				DESCRIPTION: 
				This code uses a person-day level dataset that has already been restricted to nursing home (NH) residents who are
				taking the object drug (clopiodgrel). Other eligibility criteria have also already been applied. Person-day level dataset
				contain a person ID, date, indicator for precipitant exposure (0=no exposure, 1=exposure) and indicator for outcome event
				on that day (0=no outcome event, 1=outcome event). 

				Program generates basic counts statistics, fits conditional Poisson models, and combines results for all precipitant + outcome
				pairs into a single summary document. It also preps results to be imported into R for semi-Bayes adjusted, which helps reduce
				bias from multiple estimation. 

				Conditional Poisson models are run using the %POISREG macro (and %ELEMENT helper macro) from https://sccs-studies.info/sas.html

Programmer: Adam DAmico

Date: 19Dec2024

Version History:
*/

/*SET OPTIONS SO LOG CONTAINS USEFUL INFORMATION ON MACRO RUNS*/
option mprint mlogic;

/*SET DIRECTORIES*/
libname duck "P:\nhddi\a5d\Aim1_Screening\Antiplatelets\Data\DerivedData\day_level_ow_fw_cov";
libname cardinal "P:\nhddi\a5d\Aim1_Screening\Antiplatelets\Data\DerivedData\SCCSmodels_SAS";

/*READ IN MACROS*/
%INCLUDE "P:\nhddi\a5d\Aim1_Screening\Antiplatelets\Code\macros\poisreg.sas";
%INCLUDE "P:\nhddi\a5d\Aim1_Screening\Antiplatelets\Code\macros\element.sas";

/*SET LOCATION TO SAVE LOG*/
proc printto new log="P:\nhddi\a5d\Aim1_Screening\Antiplatelets\Output\7_SCCS_models\sumstats_SAS_simplest_models_&sysdate..txt";
run;


/*CREATE EMPTY DATATSET TO PUT RESULTS INTO*/
data cardinal.sumstats_SAS_simplest_models;
	length 
		model $255
		outcome $5
		precip $255
		PARAMS $255
		Estimate 8.
		StErr 8.
		z  8.
		chisq 8.
		pvalue 8.
		ll 8.
		ul 8.
		expest 8.
		expll 8.
		expul 8.
		N_total_days 8.
		N_0exp_0out 8.
		N_1exp_0out 8.
		N_0exp_1out 8.
		N_1exp_1out 8.;
	if PARAMS='' then delete;
run;

/*OPTIONAL: GENERATE BASIC COUNTS FOR DAYS IN THE COHORT AS A WHOLE*/
%macro chicago (XYZ, drug);
proc sql;

	/*all days*/
	create table temp0 as
		select count(*) as N_total_days from duck.anal_&XYZ._all_precip2;

	/*days without outcome and without exposure*/
	create table temp1 as
		select count(*) as N_0exp_0out from duck.anal_&XYZ._all_precip2
		where day_outcome_&XYZ.=0 & day_&drug._a=0;

	/*days without outcome and with exposure*/
	create table temp2 as
		select count(*) as N_1exp_0out from duck.anal_&XYZ._all_precip2
		where day_outcome_&XYZ.=0 & day_&drug._a=1;

	/*days with outcome and without exposure*/
	create table temp3 as
		select count(*) as N_0exp_1out from duck.anal_&XYZ._all_precip2
		where day_outcome_&XYZ.=1 & day_&drug._a=0;

	/*days with outcome and with exposure*/
	create table temp4 as
		select count(*) as N_1exp_1out from duck.anal_&XYZ._all_precip2
		where day_outcome_&XYZ.=1 & day_&drug._a=1;
quit;

data basic_counts; 
	merge temp0 temp1 temp2 temp3 temp4;
run;


/*RUN CONDITIONAL POISSON MODEL*/
%poisreg(data=duck.anal_&XYZ._all_precip2,
	y=day_outcome_&XYZ.,
	covar= day_&drug._a ,
	class= ,
	offset=,
	elim=bene_id_18900,
	prntyn=N,
	outdata=SCCS_&drug.);

/*ADD SOME METADATA ABOUT MODEL*/
data SCCS_&drug.2; set SCCS_&drug.;
	if PARAMS='INT' then delete; 
	model='Simplest model (event~precip)';
	outcome="&XYZ.";
	precip="&drug.";
run;

data SCCS_&drug.3; 
	merge SCCS_&drug.2 basic_counts;
run;

/*APPEND MODEL SPECIFIC STATS TO SUMMARY FILE*/
proc append 
	base=cardinal.sumstats_SAS_simplest_models 
	data=SCCS_&drug.3;
run;

/*CLEAR WORK DIRECTORY*/
proc datasets nolist nodetails lib=work kill; 
quit;

%mend;

/*RUN MACROS FOR EACH OUTCOME+PRECIPITANT PAIR*/

*major bleed analysis;
%chicago(mb, alprazolam);
%chicago(mb, amlodipine);
%chicago(mb, amox_clav);
%chicago(mb, amox_no_clav);
%chicago(mb, atorvastatin);
%chicago(mb, azithromycin);
%chicago(mb, carvedilol);
%chicago(mb, cephalexin);
%chicago(mb, ciprofloxacin);
%chicago(mb, clonidine);
%chicago(mb, donepezil);
%chicago(mb, doxycycline);
%chicago(mb, escitalopram);
%chicago(mb, furosemide);
%chicago(mb, gabapentin);
%chicago(mb, hydralazine);
%chicago(mb, hydrocodone);
%chicago(mb, isosorbide_m);
%chicago(mb, levetiracetam);
%chicago(mb, levofloxacin);
%chicago(mb, lisinopril);
%chicago(mb, lorazepam);
%chicago(mb, losartan);
%chicago(mb, metoprolol_s);
%chicago(mb, metoprolol_t);
%chicago(mb, mirtazapine);
%chicago(mb, nitrofurantoin);
%chicago(mb, omeprazole);
%chicago(mb, ondansetron);
%chicago(mb, oxycodone);
%chicago(mb, pantoprazole);
%chicago(mb, potassium_chl);
%chicago(mb, prednisone);
%chicago(mb, sertraline);
%chicago(mb, simvastatin);
%chicago(mb, sucralfate);
%chicago(mb, sulfameth_tri);
%chicago(mb, tamsulosin);
%chicago(mb, tramadol);
%chicago(mb, trazodone);

*FRI analysis;
%chicago(fri, alprazolam);
%chicago(fri, amiodarone);
%chicago(fri, amlodipine);
%chicago(fri, amox_clav);
%chicago(fri, amox_no_clav);
%chicago(fri, atorvastatin);
%chicago(fri, azithromycin);
%chicago(fri, buspirone);
%chicago(fri, carvedilol);
%chicago(fri, cephalexin);
%chicago(fri, ciprofloxacin);
%chicago(fri, citalopram);
%chicago(fri, clonidine);
%chicago(fri, divalproex);
%chicago(fri, donepezil);
%chicago(fri, doxycycline);
%chicago(fri, duloxetine);
%chicago(fri, escitalopram);
%chicago(fri, furosemide);
%chicago(fri, gabapentin);
%chicago(fri, hydralazine);
%chicago(fri, hydrochlorothia);
%chicago(fri, hydrocodone);
%chicago(fri, isosorbide_m);
%chicago(fri, levetiracetam);
%chicago(fri, levofloxacin);
%chicago(fri, levothyroxine);
%chicago(fri, lisinopril);
%chicago(fri, lorazepam);
%chicago(fri, losartan);
%chicago(fri, meloxicam);
%chicago(fri, memantine);
%chicago(fri, metformin);
%chicago(fri, metoprolol_s);
%chicago(fri, metoprolol_t);
%chicago(fri, mirtazapine);
%chicago(fri, morphine);
%chicago(fri, nitrofurantoin);
%chicago(fri, omeprazole);
%chicago(fri, ondansetron);
%chicago(fri, oxycodone);
%chicago(fri, pantoprazole);
%chicago(fri, potassium_chl);
%chicago(fri, pravastatin);
%chicago(fri, prednisone);
%chicago(fri, quetiapine);
%chicago(fri, risperidone);
%chicago(fri, sertraline);
%chicago(fri, simvastatin);
%chicago(fri, spironolactone);
%chicago(fri, sulfameth_tri);
%chicago(fri, tamsulosin);
%chicago(fri, tramadol);
%chicago(fri, trazodone);
%chicago(fri, zolpidem);

/*RESET THE DESTINATION FOR SAS LOG TO DEFAULT*/
proc printto;
run;

/*EXPORT SUMMARY STATS TO EXCEL*/
proc export 
  data=cardinal.sumstats_SAS_simplest_models 
  dbms=xlsx 
  outfile="P:\nhddi\a5d\Aim1_Screening\Antiplatelets\Output\7_SCCS_models\sumstats_SAS_simplest_models_&sysdate..xlsx" 
  replace;
run;


/*DATA PREP TO APPLY SEMI-BAYES ADJUSTMENT IN R						
	- crete separate version of summary stats for each outcome		
	- calculate variance from standard error						
	- limit to precipitant, estimate and variance
	- export as csv 											*/

*major bleed;
data sumstats_csv_mb; set cardinal.sumstats_SAS_simplest_models;
	keep precip Est Var;
	Est=Estimate;
	Var=StErr*StErr;
	where outcome='mb';
run;

proc export 
	data=sumstats_csv_mb
	dbms=csv 
	outfile="P:\nhddi\a5d\Aim1_Screening\Antiplatelets\Output\7_SCCS_models\sumstats_SAS_simplest_models_mb_&sysdate..csv" 
	replace;
run;

*FRI;
data sumstats_csv_fri; set cardinal.sumstats_SAS_simplest_models;
	keep precip Est Var;
	Est=Estimate;
	Var=StErr*StErr;
	where outcome='fri';
run;

proc export 
	data=sumstats_csv_fri 
	dbms=csv 
	outfile="P:\nhddi\a5d\Aim1_Screening\Antiplatelets\Output\7_SCCS_models\sumstats_SAS_simplest_models_fri_&sysdate..csv" 
	replace;
run;



