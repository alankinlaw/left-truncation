dm 'log; clear; output; clear;';


%let    PATH = ----SPECIFY YOUR OWN WORKING DIRECTORY HERE----;
libname results "&PATH\documents";


options mprint;

/*********************************************************************/
/*                                                                   */
/* PROJECT : Placental Abruption, Early Pregnancy Loss, and          */
/*           Perinatal Mortality: Impact of Left Truncation Bias     */
/*                                                                   */
/* AUTHORS : Cande Ananth (ananthcv@rwjms.rutgers.edu>
/*			             											 */
/*           Hillary Graham <hlg55@rwjms.rutgers.edu> 				 */
/*			 Kinlaw, Alan <akinlaw@email.unc.edu>					 */
/*                                                                   */
/* PURPOSE : To carry out simulations for assessing the impact of    */
/*			 left truncation on studies of the association between   */
/*			 placental abruption and mortality                       */
/*                                                                   */
/* PROGRAM : ABRUPTION-MORTALITY-SIMULATIONS_GIT.sas                 */
/*                                                                   */
/* LANGUAGE: SAS on PC (9.4)                                         */
/*                                                                   */
/* CREATED : 16 July 2022                                            */
/* UPDATED : 25 August 2022                                          */ 
/*																	 */	
/* License and warranty information: 								 */
/*			 The published material is shared under a GNU General    */
/*			 Public License v3.0. It is being distributed without    */
/*           warranty of any kind, either expressed or implied.      */
/* 			 The responsibility for the interpretation and use of    */
/*		     the material lies with the reader. In no event shall    */
/*			 the Author be liable for damages arising from its use.  */
/*                                                                   */
/*********************************************************************/

* PARAMETERS:  	
	sim  		simulation scenario  
	seed  		random number seed  
	cohorts		number of cohorts to simulate
	n  			number of individuals in each cohort  
	p_a  		prevalence of placental abruption 
	lp_ap_na 		lower bound on uniform probability distribution for abnormal placentation among non abruption  
	rp_ap_na 		range of uniform probability distribution for abnormal placentation among no abruption
	lp_ap_pa 		lower bound on uniform probability distribution for abnormal placentation abruption 
	rp_ap_pa 		range of uniform probability distribution for abnormal placentation among abruption
	a_pnm  		log-odds of perinatal mortality among those without abnormal placentation   
	b_pnm_ap 		log-odds ratio of perinatal mortality for those with abnormal placentation  
	a_epl  		log-odds of early pregnancy loss among non abruption without abnormal placentation  
	b_epl_ap 		log-odds ratio of early pregnancy loss for those with abnormal placentation  
	b_epl_a 		log-odds ratio of early pregnancy loss for abruption  
	b_epl_as 		departure from multiplicativity of log-odds ratio for effects of abruption and abnormal placentation
  ;

* First, submit %LT macro ;
%macro lt(sim,seed,cohorts,n,p_z,lp_x_z0,rp_x_z0,lp_x_z1,rp_x_z1,a_c,b_c_x,b_c_z,b_c_xz,a_y,b_y_m,b_y_z,b_y_mz);

DATA work.unobservable;
	call streaminit(&seed);
	do sample = 1 to &cohorts;
		do id = 1 to &n;

		a = rand("uniform");
		b = rand("uniform");

		* set p_z and generate z=1 and z=0, 
		  based on random uniform # a ;
			if a ge &p_z then z = 0;
			else z = 1;

		* generate individual probabilities of x=1, 
		  based on random uniform # b, and z-specific risks of x ;
			if z = 0 then p = b*&rp_x_z0 + &lp_x_z0;
			else if z = 1 then p = b*&rp_x_z1 + &lp_x_z1;
		* do weighted coin flip for x=1 for each individual ;
			x = rand("binomial",p,1);

		* use logistic regression to simulate predicted probabilities of c=1 ;
			prob_c = 1/(1+exp(-(&a_c + &b_c_x*x + &b_c_z*z + &b_c_xz*x*z)));
		* do weighted coin flip for early pregnancy loss (epl) for each individual ;
			c = rand("binomial",prob_c,1);

		* deterministically code a post-competing-event variable M and outcome Y ;
			if c = 1 then do; 
				m = .;
				y = 1; *deterministic because early loss is a type of Y=1;
				end;
			else if c = 0 then do; 
				m = x; *we assume all non-early losses had same value for abruption throughout pregnancy ;
				* use logistic regression to simulate individual probabilities of y=1 ;
					prob_y = 1/(1+exp(-(&a_y + &b_y_m*m + &b_y_z*z + &b_y_mz*m*z)));
				* do weighted coin flip for preeclampsia (pe) for each individual ;
					y = rand("binomial",prob_y,1);
				end;

		* code a mutually exclusive categorical variable for x and z associated with y ;
			if x = 1 and z = 1 then do; xz01 = 0; xz10 = 0; xz11 = 1; end;
			else if x = 1 and z = 0 then do; xz01 = 0; xz10 = 1; xz11 = 0; end;
			else if x = 0 and z = 1 then do; xz01 = 1; xz10 = 0; xz11 = 0; end;
			else if x = 0 and z = 0 then do; xz01 = 0; xz10 = 0; xz11 = 0; end;

		output; 
		end;
	end; run;

* remove all competing events from the 'observed' dataset 
	(i.e., make left truncation occur) ;

DATA work.observed;
	set work.unobservable;
	if c = 1 then delete;
run;

* calculate risks and contrasts for total effect of x on y 
	stratified on z in unobservable and observed populations ;
	ods output Estimates=est_unobservable ;
		proc genmod data = work.unobservable descending; 
		by sample;
		model y = xz10 xz01 xz11 / link=log dist=bin ;
		estimate 'risk xz00' int 1 / exp; 
		estimate 'risk xz10' int 1 xz10 1 / exp; 
		estimate 'risk xz01' int 1 xz01 1 / exp;
		estimate 'risk xz11' int 1 xz11 1 / exp;
		run; ods rtf close;
		data select_unobservable_1;
			set est_unobservable;
			where label in('risk xz00','risk xz10','risk xz01','risk xz11');
			keep sample label meanestimate;
			run;
			proc transpose data = select_unobservable_1 out = select_unobservable_2;
				by sample;
				run;
				proc sql;
					create table select_unobservable_3_&sim as
					select 	sample, 
							col1 as risk_xz00, 
							col2 as risk_xz10, 
							col3 as risk_xz01, 
							col4 as risk_xz11,
							col2/col1 as rr_x_z0,
							col4/col3 as rr_x_z1,
							col2-col1 as rd_x_z0,
							col4-col3 as rd_x_z1
					from select_unobservable_2;
					quit;
					proc means data=select_unobservable_3_&sim mean std noprint;
						var risk_xz00
							risk_xz10
							risk_xz01
							risk_xz11
							rr_x_z0
							rr_x_z1
							rd_x_z0
							rd_x_z1;
						output out=results.xy_unobservable_&n._&sim; run;

ods output Estimates=est_observed ;
		proc genmod data = observed descending; 
		by sample;
		model y = xz10 xz01 xz11 / link=log dist=bin ;
		estimate 'risk xz00' int 1 / exp; 
		estimate 'risk xz10' int 1 xz10 1 / exp; 
		estimate 'risk xz01' int 1 xz01 1 / exp;
		estimate 'risk xz11' int 1 xz11 1 / exp;
		run; ods rtf close;
		data select_observed_1;
			set est_observed;
			where label in('risk xz00','risk xz10','risk xz01','risk xz11');
			keep sample label meanestimate;
			run;
			proc transpose data = select_observed_1 out = select_observed_2;
				by sample;
				run;
				proc sql;
					create table select_observed_3_&sim as
					select 	sample, 
							col1 as risk_xz00, 
							col2 as risk_xz10, 
							col3 as risk_xz01, 
							col4 as risk_xz11,
							col2/col1 as rr_x_z0,
							col4/col3 as rr_x_z1,
							col2-col1 as rd_x_z0,
							col4-col3 as rd_x_z1
					from select_observed_2;
					quit;
					proc means data=select_observed_3_&sim mean std noprint;
						var risk_xz00
							risk_xz10
							risk_xz01
							risk_xz11
							rr_x_z0
							rr_x_z1
							rd_x_z0
							rd_x_z1;
						output out=results.xy_observed_&n._&sim; run;



* calculate risks and contrasts for total effect of z on y ;
	ods output Estimates=est_unobservable ;
		proc genmod data = unobservable descending; 
		by sample;
		model y = z / link=log dist=bin ;
		estimate 'risk z0' int 1 / exp; 
		estimate 'risk z1' int 1 z 1 / exp; 
		run; ods rtf close;
		data select_unobservable_1;
			set est_unobservable;
			where label in('risk z0','risk z1');
			keep sample label meanestimate;
			run;
			proc transpose data = select_unobservable_1 out = select_unobservable_2;
				by sample;
				run;
				proc sql;
					create table select_unobservable_3_&sim as
					select 	sample, 
							col1 as risk_z0, 
							col2 as risk_z1, 
							col2/col1 as rr_z,
							col2-col1 as rd_z
					from select_unobservable_2;
					quit;
					proc means data=select_unobservable_3_&sim mean std noprint;
						var risk_z0
							risk_z1
							rr_z
							rd_z;
						output out=results.zy_unobservable_&n._&sim; run;

	ods output Estimates=est_observed ;
		proc genmod data = observed descending; 
		by sample;
		model y = z / link=log dist=bin ;
		estimate 'risk z0' int 1 / exp; 
		estimate 'risk z1' int 1 z 1 / exp; 
		run; ods rtf close;
		data select_observed_1;
			set est_observed;
			where label in('risk z0','risk z1');
			keep sample label meanestimate;
			run;
			proc transpose data = select_observed_1 out = select_observed_2;
				by sample;
				run;
				proc sql;
					create table select_observed_3_&sim as
					select 	sample, 
							col1 as risk_z0, 
							col2 as risk_z1, 
							col2/col1 as rr_z,
							col2-col1 as rd_z
					from select_observed_2;
					quit;
					proc means data=select_observed_3_&sim mean std noprint;
						var risk_z0
							risk_z1
							rr_z
							rd_z;
						output out=results.zy_observed_&n._&sim;
run;

%mend;


*Second, specify parameter values and execute macro -- examples are shown below from manuscript assumptions;
* sim setup  1 in manuscript; %lt(sim=1001, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=0.747214401830221, b_c_xz=0.0637158143861078, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=0.701449344272941, b_y_mz=0.137112006695)
* sim setup  2 in manuscript; %lt(sim=1013, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=0.747214401830221, b_c_xz=0.0637158143861078, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=1.41141053587434, b_y_mz=0.495772744598063)
* sim setup  3 in manuscript; %lt(sim=1002, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=1.70767635201751, b_c_z=0.747214401830221, b_c_xz=0.288877529856554, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=0.701449344272941, b_y_mz=0.137112006695)
* sim setup  4 in manuscript; %lt(sim=1014, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=1.70767635201751, b_c_z=0.747214401830221, b_c_xz=0.288877529856554, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=1.41141053587434, b_y_mz=0.495772744598063)
* sim setup  5 in manuscript; %lt(sim=1003, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=1.55814461804655, b_c_xz=0.233614851181505, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=0.701449344272941, b_y_mz=0.137112006695)
* sim setup  6 in manuscript; %lt(sim=1015, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=1.55814461804655, b_c_xz=0.233614851181505, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=1.41141053587434, b_y_mz=0.495772744598063)
* sim setup  7 in manuscript; %lt(sim=1004, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=1.70767635201751, b_c_z=1.55814461804655, b_c_xz=1.8758425864386, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=0.701449344272941, b_y_mz=0.137112006695)
* sim setup  8 in manuscript; %lt(sim=1016, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=1.70767635201751, b_c_z=1.55814461804655, b_c_xz=1.8758425864386, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=1.41141053587434, b_y_mz=0.495772744598063)
* sim setup  9 in manuscript; %lt(sim=1005, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=2.09714111877924, b_c_xz=0.505548566665147, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=0.701449344272941, b_y_mz=0.137112006695)
* sim setup 10 in manuscript; %lt(sim=1017, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=2.09714111877924, b_c_xz=0.505548566665147, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=1.41141053587434, b_y_mz=0.495772744598063)
* sim setup 11 in manuscript; %lt(sim=1007, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=0.747214401830221, b_c_xz=0.0637158143861078, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=0.705367291894721, b_y_mz=0.133183027765782)
* sim setup 12 in manuscript; %lt(sim=1019, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=0.747214401830221, b_c_xz=0.0637158143861078, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=1.42341202407639, b_y_mz=0.483723083052569)
* sim setup 13 in manuscript; %lt(sim=1008, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=1.70767635201751, b_c_z=0.747214401830221, b_c_xz=0.288877529856554, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=0.705367291894721, b_y_mz=0.133183027765782)
* sim setup 14 in manuscript; %lt(sim=1020, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=1.70767635201751, b_c_z=0.747214401830221, b_c_xz=0.288877529856554, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=1.42341202407639, b_y_mz=0.483723083052569)
* sim setup 15 in manuscript; %lt(sim=1009, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=1.55814461804655, b_c_xz=0.233614851181505, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=0.705367291894721, b_y_mz=0.133183027765782)
* sim setup 16 in manuscript; %lt(sim=1021, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=1.55814461804655, b_c_xz=0.233614851181505, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=1.42341202407639, b_y_mz=0.483723083052569)
* sim setup 17 in manuscript; %lt(sim=1010, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=1.70767635201751, b_c_z=1.55814461804655, b_c_xz=1.8758425864386, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=0.705367291894721, b_y_mz=0.133183027765782)
* sim setup 18 in manuscript; %lt(sim=1022, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=1.70767635201751, b_c_z=1.55814461804655, b_c_xz=1.8758425864386, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=1.42341202407639, b_y_mz=0.483723083052569)
* sim setup 19 in manuscript; %lt(sim=1011, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=2.09714111877924, b_c_xz=0.505548566665147, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=0.705367291894721, b_y_mz=0.133183027765782)
* sim setup 20 in manuscript; %lt(sim=1023, seed=123, cohorts=200, n=100000, p_z=0.1, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=2.09714111877924, b_c_xz=0.505548566665147, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=1.42341202407639, b_y_mz=0.483723083052569)
* sim setup 21 in manuscript; %lt(sim=1025, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=0.747214401830221, b_c_xz=0.0637158143861078, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=0.701449344272941, b_y_mz=0.137112006695)
* sim setup 22 in manuscript; %lt(sim=1037, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=0.747214401830221, b_c_xz=0.0637158143861078, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=1.41141053587434, b_y_mz=0.495772744598063)
* sim setup 23 in manuscript; %lt(sim=1026, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=1.70767635201751, b_c_z=0.747214401830221, b_c_xz=0.288877529856554, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=0.701449344272941, b_y_mz=0.137112006695)
* sim setup 24 in manuscript; %lt(sim=1038, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=1.70767635201751, b_c_z=0.747214401830221, b_c_xz=0.288877529856554, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=1.41141053587434, b_y_mz=0.495772744598063)
* sim setup 25 in manuscript; %lt(sim=1027, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=1.55814461804655, b_c_xz=0.233614851181505, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=0.701449344272941, b_y_mz=0.137112006695)
* sim setup 26 in manuscript; %lt(sim=1039, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=1.55814461804655, b_c_xz=0.233614851181505, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=1.41141053587434, b_y_mz=0.495772744598063)
* sim setup 27 in manuscript; %lt(sim=1028, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=1.70767635201751, b_c_z=1.55814461804655, b_c_xz=1.8758425864386, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=0.701449344272941, b_y_mz=0.137112006695)
* sim setup 28 in manuscript; %lt(sim=1040, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=1.70767635201751, b_c_z=1.55814461804655, b_c_xz=1.8758425864386, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=1.41141053587434, b_y_mz=0.495772744598063)
* sim setup 29 in manuscript; %lt(sim=1029, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=2.09714111877924, b_c_xz=0.505548566665147, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=0.701449344272941, b_y_mz=0.137112006695)
* sim setup 30 in manuscript; %lt(sim=1041, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=2.09714111877924, b_c_xz=0.505548566665147, a_y=-4.79538731978483, b_y_m=2.79539187222536, b_y_z=1.41141053587434, b_y_mz=0.495772744598063)
* sim setup 31 in manuscript; %lt(sim=1031, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=0.747214401830221, b_c_xz=0.0637158143861078, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=0.705367291894721, b_y_mz=0.133183027765782)
* sim setup 32 in manuscript; %lt(sim=1043, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=0.747214401830221, b_c_xz=0.0637158143861078, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=1.42341202407639, b_y_mz=0.483723083052569)
* sim setup 33 in manuscript; %lt(sim=1032, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=1.70767635201751, b_c_z=0.747214401830221, b_c_xz=0.288877529856554, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=0.705367291894721, b_y_mz=0.133183027765782)
* sim setup 34 in manuscript; %lt(sim=1044, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=1.70767635201751, b_c_z=0.747214401830221, b_c_xz=0.288877529856554, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=1.42341202407639, b_y_mz=0.483723083052569)
* sim setup 35 in manuscript; %lt(sim=1033, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=1.55814461804655, b_c_xz=0.233614851181505, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=0.705367291894721, b_y_mz=0.133183027765782)
* sim setup 36 in manuscript; %lt(sim=1045, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=1.55814461804655, b_c_xz=0.233614851181505, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=1.42341202407639, b_y_mz=0.483723083052569)
* sim setup 37 in manuscript; %lt(sim=1034, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=1.70767635201751, b_c_z=1.55814461804655, b_c_xz=1.8758425864386, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=0.705367291894721, b_y_mz=0.133183027765782)
* sim setup 38 in manuscript; %lt(sim=1046, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=1.70767635201751, b_c_z=1.55814461804655, b_c_xz=1.8758425864386, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=1.42341202407639, b_y_mz=0.483723083052569)
* sim setup 39 in manuscript; %lt(sim=1035, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=2.09714111877924, b_c_xz=0.505548566665147, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=0.705367291894721, b_y_mz=0.133183027765782)
* sim setup 40 in manuscript; %lt(sim=1047, seed=123, cohorts=200, n=100000, p_z=0.2, lp_x_z0=0.001, rp_x_z0=0.019, lp_x_z1=0.005, rp_x_z1=0.045, a_c=-2.94443897916644, b_c_x=0.747214401830221, b_c_z=2.09714111877924, b_c_xz=0.505548566665147, a_y=-4.41077604795987, b_y_m=2.41071011825686, b_y_z=1.42341202407639, b_y_mz=0.483723083052569)




*Third, compile results from each simulation setup and output into excel;
%macro extractmean(association,version,n,sim);
data &association._&version._&n._&sim._2;
	set results.&association._&version._&n._&sim;
	where _stat_ = 'MEAN';
	n_per_iteration = &n;
	num_iterations = _freq_;
	sim = &sim;
	drop _type_ _freq_ _stat_;
	run;
%mend;

%extractmean(xy,unobservable,100000,1001);
%extractmean(xy,unobservable,100000,1002);
%extractmean(xy,unobservable,100000,1003);
%extractmean(xy,unobservable,100000,1004);
%extractmean(xy,unobservable,100000,1005);
%extractmean(xy,unobservable,100000,1007);
%extractmean(xy,unobservable,100000,1008);
%extractmean(xy,unobservable,100000,1009);
%extractmean(xy,unobservable,100000,1010);
%extractmean(xy,unobservable,100000,1011);
%extractmean(xy,unobservable,100000,1013);
%extractmean(xy,unobservable,100000,1014);
%extractmean(xy,unobservable,100000,1015);
%extractmean(xy,unobservable,100000,1016);
%extractmean(xy,unobservable,100000,1017);
%extractmean(xy,unobservable,100000,1019);
%extractmean(xy,unobservable,100000,1020);
%extractmean(xy,unobservable,100000,1021);
%extractmean(xy,unobservable,100000,1022);
%extractmean(xy,unobservable,100000,1023);
%extractmean(xy,unobservable,100000,1025);
%extractmean(xy,unobservable,100000,1026);
%extractmean(xy,unobservable,100000,1027);
%extractmean(xy,unobservable,100000,1028);
%extractmean(xy,unobservable,100000,1029);
%extractmean(xy,unobservable,100000,1031);
%extractmean(xy,unobservable,100000,1032);
%extractmean(xy,unobservable,100000,1033);
%extractmean(xy,unobservable,100000,1034);
%extractmean(xy,unobservable,100000,1035);
%extractmean(xy,unobservable,100000,1037);
%extractmean(xy,unobservable,100000,1038);
%extractmean(xy,unobservable,100000,1039);
%extractmean(xy,unobservable,100000,1040);
%extractmean(xy,unobservable,100000,1041);
%extractmean(xy,unobservable,100000,1043);
%extractmean(xy,unobservable,100000,1044);
%extractmean(xy,unobservable,100000,1045);
%extractmean(xy,unobservable,100000,1046);
%extractmean(xy,unobservable,100000,1047);

DATA results.xy_mean_unobservable_100000;
	set xy_unobservable_100000_1001_2
		xy_unobservable_100000_1002_2
		xy_unobservable_100000_1003_2
		xy_unobservable_100000_1004_2
		xy_unobservable_100000_1005_2
		xy_unobservable_100000_1007_2
		xy_unobservable_100000_1008_2
		xy_unobservable_100000_1009_2
		xy_unobservable_100000_1010_2
		xy_unobservable_100000_1011_2
		xy_unobservable_100000_1013_2
		xy_unobservable_100000_1014_2
		xy_unobservable_100000_1015_2
		xy_unobservable_100000_1016_2
		xy_unobservable_100000_1017_2
		xy_unobservable_100000_1019_2
		xy_unobservable_100000_1020_2
		xy_unobservable_100000_1021_2
		xy_unobservable_100000_1022_2
		xy_unobservable_100000_1023_2
		xy_unobservable_100000_1025_2
		xy_unobservable_100000_1026_2
		xy_unobservable_100000_1027_2
		xy_unobservable_100000_1028_2
		xy_unobservable_100000_1029_2
		xy_unobservable_100000_1031_2
		xy_unobservable_100000_1032_2
		xy_unobservable_100000_1033_2
		xy_unobservable_100000_1034_2
		xy_unobservable_100000_1035_2
		xy_unobservable_100000_1037_2
		xy_unobservable_100000_1038_2
		xy_unobservable_100000_1039_2
		xy_unobservable_100000_1040_2
		xy_unobservable_100000_1041_2
		xy_unobservable_100000_1043_2
		xy_unobservable_100000_1044_2
		xy_unobservable_100000_1045_2
		xy_unobservable_100000_1046_2
		xy_unobservable_100000_1047_2
		;
run;


%extractmean(xy,observed,100000,1001);
%extractmean(xy,observed,100000,1002);
%extractmean(xy,observed,100000,1003);
%extractmean(xy,observed,100000,1004);
%extractmean(xy,observed,100000,1005);
%extractmean(xy,observed,100000,1007);
%extractmean(xy,observed,100000,1008);
%extractmean(xy,observed,100000,1009);
%extractmean(xy,observed,100000,1010);
%extractmean(xy,observed,100000,1011);
%extractmean(xy,observed,100000,1013);
%extractmean(xy,observed,100000,1014);
%extractmean(xy,observed,100000,1015);
%extractmean(xy,observed,100000,1016);
%extractmean(xy,observed,100000,1017);
%extractmean(xy,observed,100000,1019);
%extractmean(xy,observed,100000,1020);
%extractmean(xy,observed,100000,1021);
%extractmean(xy,observed,100000,1022);
%extractmean(xy,observed,100000,1023);
%extractmean(xy,observed,100000,1025);
%extractmean(xy,observed,100000,1026);
%extractmean(xy,observed,100000,1027);
%extractmean(xy,observed,100000,1028);
%extractmean(xy,observed,100000,1029);
%extractmean(xy,observed,100000,1031);
%extractmean(xy,observed,100000,1032);
%extractmean(xy,observed,100000,1033);
%extractmean(xy,observed,100000,1034);
%extractmean(xy,observed,100000,1035);
%extractmean(xy,observed,100000,1037);
%extractmean(xy,observed,100000,1038);
%extractmean(xy,observed,100000,1039);
%extractmean(xy,observed,100000,1040);
%extractmean(xy,observed,100000,1041);
%extractmean(xy,observed,100000,1043);
%extractmean(xy,observed,100000,1044);
%extractmean(xy,observed,100000,1045);
%extractmean(xy,observed,100000,1046);
%extractmean(xy,observed,100000,1047);

data results.xy_mean_observed_100000;
	set xy_observed_100000_1001_2
		xy_observed_100000_1002_2
		xy_observed_100000_1003_2
		xy_observed_100000_1004_2
		xy_observed_100000_1005_2
		xy_observed_100000_1007_2
		xy_observed_100000_1008_2
		xy_observed_100000_1009_2
		xy_observed_100000_1010_2
		xy_observed_100000_1011_2
		xy_observed_100000_1013_2
		xy_observed_100000_1014_2
		xy_observed_100000_1015_2
		xy_observed_100000_1016_2
		xy_observed_100000_1017_2
		xy_observed_100000_1019_2
		xy_observed_100000_1020_2
		xy_observed_100000_1021_2
		xy_observed_100000_1022_2
		xy_observed_100000_1023_2
		xy_observed_100000_1025_2
		xy_observed_100000_1026_2
		xy_observed_100000_1027_2
		xy_observed_100000_1028_2
		xy_observed_100000_1029_2
		xy_observed_100000_1031_2
		xy_observed_100000_1032_2
		xy_observed_100000_1033_2
		xy_observed_100000_1034_2
		xy_observed_100000_1035_2
		xy_observed_100000_1037_2
		xy_observed_100000_1038_2
		xy_observed_100000_1039_2
		xy_observed_100000_1040_2
		xy_observed_100000_1041_2
		xy_observed_100000_1043_2
		xy_observed_100000_1044_2
		xy_observed_100000_1045_2
		xy_observed_100000_1046_2
		xy_observed_100000_1047_2
		;
run;

PROC EXPORT DATA= results.xy_mean_unobservable_100000
            OUTFILE= "&PATH\documents\xy_mean_unobservable_100000.xlsx" 
            DBMS=xlsx REPLACE;
			SHEET="sheet"; 
			NEWFILE=YES;
run;

PROC EXPORT DATA= results.xy_mean_observed_100000
            OUTFILE= "&PATH\documents\xy_mean_observed_100000.xlsx" 
            DBMS=xlsx REPLACE;
			SHEET="sheet"; 
			NEWFILE=YES;
run;


%extractmean(zy,unobservable,100000,1001);
%extractmean(zy,unobservable,100000,1002);
%extractmean(zy,unobservable,100000,1003);
%extractmean(zy,unobservable,100000,1004);
%extractmean(zy,unobservable,100000,1005);
%extractmean(zy,unobservable,100000,1007);
%extractmean(zy,unobservable,100000,1008);
%extractmean(zy,unobservable,100000,1009);
%extractmean(zy,unobservable,100000,1010);
%extractmean(zy,unobservable,100000,1011);
%extractmean(zy,unobservable,100000,1013);
%extractmean(zy,unobservable,100000,1014);
%extractmean(zy,unobservable,100000,1015);
%extractmean(zy,unobservable,100000,1016);
%extractmean(zy,unobservable,100000,1017);
%extractmean(zy,unobservable,100000,1019);
%extractmean(zy,unobservable,100000,1020);
%extractmean(zy,unobservable,100000,1021);
%extractmean(zy,unobservable,100000,1022);
%extractmean(zy,unobservable,100000,1023);
%extractmean(zy,unobservable,100000,1025);
%extractmean(zy,unobservable,100000,1026);
%extractmean(zy,unobservable,100000,1027);
%extractmean(zy,unobservable,100000,1028);
%extractmean(zy,unobservable,100000,1029);
%extractmean(zy,unobservable,100000,1031);
%extractmean(zy,unobservable,100000,1032);
%extractmean(zy,unobservable,100000,1033);
%extractmean(zy,unobservable,100000,1034);
%extractmean(zy,unobservable,100000,1035);
%extractmean(zy,unobservable,100000,1037);
%extractmean(zy,unobservable,100000,1038);
%extractmean(zy,unobservable,100000,1039);
%extractmean(zy,unobservable,100000,1040);
%extractmean(zy,unobservable,100000,1041);
%extractmean(zy,unobservable,100000,1043);
%extractmean(zy,unobservable,100000,1044);
%extractmean(zy,unobservable,100000,1045);
%extractmean(zy,unobservable,100000,1046);
%extractmean(zy,unobservable,100000,1047);

data results.zy_mean_unobservable_100000;
	set zy_unobservable_100000_1001_2
		zy_unobservable_100000_1002_2
		zy_unobservable_100000_1003_2
		zy_unobservable_100000_1004_2
		zy_unobservable_100000_1005_2
		zy_unobservable_100000_1007_2
		zy_unobservable_100000_1008_2
		zy_unobservable_100000_1009_2
		zy_unobservable_100000_1010_2
		zy_unobservable_100000_1011_2
		zy_unobservable_100000_1013_2
		zy_unobservable_100000_1014_2
		zy_unobservable_100000_1015_2
		zy_unobservable_100000_1016_2
		zy_unobservable_100000_1017_2
		zy_unobservable_100000_1019_2
		zy_unobservable_100000_1020_2
		zy_unobservable_100000_1021_2
		zy_unobservable_100000_1022_2
		zy_unobservable_100000_1023_2
		zy_unobservable_100000_1025_2
		zy_unobservable_100000_1026_2
		zy_unobservable_100000_1027_2
		zy_unobservable_100000_1028_2
		zy_unobservable_100000_1029_2
		zy_unobservable_100000_1031_2
		zy_unobservable_100000_1032_2
		zy_unobservable_100000_1033_2
		zy_unobservable_100000_1034_2
		zy_unobservable_100000_1035_2
		zy_unobservable_100000_1037_2
		zy_unobservable_100000_1038_2
		zy_unobservable_100000_1039_2
		zy_unobservable_100000_1040_2
		zy_unobservable_100000_1041_2
		zy_unobservable_100000_1043_2
		zy_unobservable_100000_1044_2
		zy_unobservable_100000_1045_2
		zy_unobservable_100000_1046_2
		zy_unobservable_100000_1047_2
		;
run;


%extractmean(zy,observed,100000,1001);
%extractmean(zy,observed,100000,1002);
%extractmean(zy,observed,100000,1003);
%extractmean(zy,observed,100000,1004);
%extractmean(zy,observed,100000,1005);
%extractmean(zy,observed,100000,1007);
%extractmean(zy,observed,100000,1008);
%extractmean(zy,observed,100000,1009);
%extractmean(zy,observed,100000,1010);
%extractmean(zy,observed,100000,1011);
%extractmean(zy,observed,100000,1013);
%extractmean(zy,observed,100000,1014);
%extractmean(zy,observed,100000,1015);
%extractmean(zy,observed,100000,1016);
%extractmean(zy,observed,100000,1017);
%extractmean(zy,observed,100000,1019);
%extractmean(zy,observed,100000,1020);
%extractmean(zy,observed,100000,1021);
%extractmean(zy,observed,100000,1022);
%extractmean(zy,observed,100000,1023);
%extractmean(zy,observed,100000,1025);
%extractmean(zy,observed,100000,1026);
%extractmean(zy,observed,100000,1027);
%extractmean(zy,observed,100000,1028);
%extractmean(zy,observed,100000,1029);
%extractmean(zy,observed,100000,1031);
%extractmean(zy,observed,100000,1032);
%extractmean(zy,observed,100000,1033);
%extractmean(zy,observed,100000,1034);
%extractmean(zy,observed,100000,1035);
%extractmean(zy,observed,100000,1037);
%extractmean(zy,observed,100000,1038);
%extractmean(zy,observed,100000,1039);
%extractmean(zy,observed,100000,1040);
%extractmean(zy,observed,100000,1041);
%extractmean(zy,observed,100000,1043);
%extractmean(zy,observed,100000,1044);
%extractmean(zy,observed,100000,1045);
%extractmean(zy,observed,100000,1046);
%extractmean(zy,observed,100000,1047);

data results.zy_mean_observed_100000;
	set zy_observed_100000_1001_2
		zy_observed_100000_1002_2
		zy_observed_100000_1003_2
		zy_observed_100000_1004_2
		zy_observed_100000_1005_2
		zy_observed_100000_1007_2
		zy_observed_100000_1008_2
		zy_observed_100000_1009_2
		zy_observed_100000_1010_2
		zy_observed_100000_1011_2
		zy_observed_100000_1013_2
		zy_observed_100000_1014_2
		zy_observed_100000_1015_2
		zy_observed_100000_1016_2
		zy_observed_100000_1017_2
		zy_observed_100000_1019_2
		zy_observed_100000_1020_2
		zy_observed_100000_1021_2
		zy_observed_100000_1022_2
		zy_observed_100000_1023_2
		zy_observed_100000_1025_2
		zy_observed_100000_1026_2
		zy_observed_100000_1027_2
		zy_observed_100000_1028_2
		zy_observed_100000_1029_2
		zy_observed_100000_1031_2
		zy_observed_100000_1032_2
		zy_observed_100000_1033_2
		zy_observed_100000_1034_2
		zy_observed_100000_1035_2
		zy_observed_100000_1037_2
		zy_observed_100000_1038_2
		zy_observed_100000_1039_2
		zy_observed_100000_1040_2
		zy_observed_100000_1041_2
		zy_observed_100000_1043_2
		zy_observed_100000_1044_2
		zy_observed_100000_1045_2
		zy_observed_100000_1046_2
		zy_observed_100000_1047_2
		;

run;

PROC EXPORT DATA= results.zy_mean_unobservable_100000
            OUTFILE= "&PATH\documents\zy_mean_unobservable_100000.xlsx" 
            DBMS=xlsx REPLACE;
			SHEET="sheet"; 
			NEWFILE=YES;
run;

PROC EXPORT DATA= results.zy_mean_observed_100000
            OUTFILE= "&PATH\documents\zy_mean_observed_100000.xlsx" 
            DBMS=xlsx REPLACE;
			SHEET="sheet"; 
			NEWFILE=YES;
run;

quit;
