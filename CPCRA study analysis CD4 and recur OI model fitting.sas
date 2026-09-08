proc datasets lib=work kill; run; quit;
DM 'log;clear;out;clear;odsresults;clear;';
***********************************************************************************
*************************Model I***************************************************;
*%LET InD=...;
*%LET out=...;

libname ddiddc "...\DDIDDC";

* Result for correlated random effects;


data ddiddc;
set ddiddc.ddiddc;
stratum=stratum-1;
gender=gender-1;
hemobl=hemobl-12;
id=seq;
trt=randgrp-1;
run;

data one;
set ddiddc;
array cd4_all {11} cd4bl cd402 cd404 cd406 cd408 cd410 cd412 cd414 cd416 cd418 cd420;
aa=1;
retain CD4;
do i=1 to 11;
*if cd4_all{i}^=. then cd4=cd4_all{i};
*else if cd4_all{i}=. then cd4=cd4;
cd4=cd4_all{i};
	stoptime=2*(i-1);
	event=3;
	output;
end;
run;

Data one;
set one;
if cd4^=.;
cd4=cd4/100;
run;

Data CD4; set one;
Keep ID CD4 Stoptime Event;Run;



*****CD4 dataset***;
***************Recurrent event*******;
data two;
set ddiddc;
aa=1;
array event_all{20} MAC PCP1 PCP2 PCP3 CANE1 CANE2 CANE3 CMV WAST 
KSV ADC CRYC TB LYMP PML TOXO CRYS OMYC HIST HZ1;
array time_all{20} T2MAC T2PCP1 T2PCP2 T2PCP3 T2CANE1 T2CANE2 T2CANE3 T2CMV T2WAST 
T2KSV T2ADC T2CRYC T2TB T2LYMP T2PML T2TOXO T2CRYS T2OMYC T2HIST T2HZ1;

do i=1 to 20;
	event=event_all{i};
	stoptime=time_all{i}/30;
	output;
end;
run;

data two;
set two;
if event=1;
run;

proc sort data=two;
by id stoptime;
run;
* Get the quantiles for recurrent events;
proc univariate data=two noprint;
var stoptime; 
output out=quant_r pctlpts=0 25 50 75 100 pctlpre=qr; 
run;
data quant_r;
set quant_r;
aa=1;
run;

* The dataset for death event;
data three;
set ddiddc;
by seq;
if first.seq;
stoptime=t2death/30;
event=death*2;		* Set event=2 for death;
aa=1;
run;
* Get the quantiels for death time;
proc univariate data=three noprint;
var stoptime; 
output out=quant_d pctlpts=0 25 50 75 100 pctlpre=qd; 
where event=2;
run;
data quant_d;
set quant_d;
aa=1;
run;
* Merge the recurrent and death event times;
data four;
set one two three;
run;

proc sort data=four;
by id stoptime;
run;
* Merge data with the quantiles;

data four2;
Retain ID Age Gender Race Trt Stratum Prevoi Hemobl CD4 stoptime event;
merge four quant_r quant_d;
by aa;
if event=3 then event=4;
if event=0 then event=3;
Keep ID  Age Gender Race Trt Stratum Prevoi Hemobl CD4 stoptime event qr0--qd100 aa;
run;


* Calculate W for each stoptime;

Data four3; set four2; by id;
retain CD4_new;
if first.id then CD4_new=CD4;
else if CD4^=. then CD4_new=CD4;
drop CD4;
rename CD4_new=CD4;
run;

/***************Need create a new variable W_dur for the last W value in the duration************/

proc sort data=four3; by event ID;run;

data four3_2; set four3;
by event id;
retain last_CD4;
if first.id then do;
	   CD4_dur=0; 
	   last_CD4=CD4;	
end;
else do;
	CD4_dur=last_CD4;
	last_CD4=CD4;
end;
if event^=4 then CD4_dur=CD4;
drop last_CD4;
run;


* Calculate the number of recurrent events as a mediator for death;
proc sort data=four3_2; by ID stoptime;run;

/*calculate the nevent only for recurrent and terminal event*/
data four4;set four3_2;
by id;
retain nevent_tmp;
if first.id then do;
	   nevent_tmp=0;   	
end;
else if event=1 then do;
	nevent_tmp=nevent_tmp+1;
end;
else if event^=1 then do;
nevent_tmp=nevent_tmp;
end;

if event=1 then nevent=max(0,nevent_tmp-1);
else nevent=nevent_tmp;
drop nevent_tmp;
run;

proc sort data=four4; by ID stoptime;run;
/*calculate the start and stop for all the events*/
data four5;set four4;
by id;
retain last_stop;
if first.id then do;
	   start=0;
	   stop=stoptime;
	   last_stop=stoptime;	
end;
else do;
	start=last_stop;
	stop=stoptime;
	last_stop=stoptime;
end;
run;

* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data five;
set four5;
array quant_r {5} qr0 qr25 qr50 qr75 qr100;
array quant_d {5} qd0 qd25 qd50 qd75 qd100;

array dur_r {4} dur_r1-dur_r4;

array event_r {4} event_r1-event_r4;
array event_d {4} event_d1-event_d4;

do i=1 to 4;
	dur_r{i}=0;
	
	event_r{i}=0;
	event_d{i}=0;
end;

last_start=start;

do i=2 to 5;
	if stop>quant_r{i} then do;
		if last_start<= quant_r{i} then do;
			dur_r{i-1}=quant_r{i}-last_start;
			last_start=quant_r{i};
		end;
	end;
	else do;
		dur_r{i-1}=stop - last_start;
		i=5;
	end;
end;
* For recurrent event;
if event=1 then do;
	do i=2 to 5;
		if stoptime<=quant_r{i} then do;
			event_r{i-1}=1;
			i=5;
		end;
	end;

end;

/* If death or censored observation */
else if event in (2,3) then do; 
	do i=2 to 5;
		if stoptime<=quant_d{i} then do;
			event_d{i-1}=(event=2);
			i=5;
		end;
	end;
end;
drop i;
run;

* For each of mediator "nevent", we need to calculate the duration in each quantile interval of death time;
data five2;
set five;
array quant {5} qd0 qd25 qd50 qd75 qd100;
array dur {4} dur1-dur4;

last_start=start;
do i=1 to 4;
	dur{i}=0;
end;

do i=2 to 5;
	if stop>quant{i} then do;
		if last_start<= quant{i} then do;
			dur{i-1}=quant{i}-last_start;
			last_start=quant{i};
		end;
	end;
	else do;
		dur{i-1}=stop - last_start;
		i=5;
	end;
end;
run;
/*
*****Estimate the initial value for the CD4;
Data CD4; set five2; if event=4;run;
proc mixed data=CD4;
class ID;
model CD4=prevoi gender trt stratum hemobl stoptime / s; 
random ID ;
repeated / subject=ID Type=AR(1);

run;
*/

proc nlmixed data=five2 qpoints=5 MAXITER=5000 corr cov;

parms log_lam_m1=-3 log_lam_m2=-3 log_lam_m3=-3 log_lam_m4=-3
		log_lam_y1=-5.5 log_lam_y2=-5 log_lam_y3=-4.6 log_lam_y4=-3.9
		betax1=0.6 betax2=-0.6 betaz=0 betax3=-0.1 betax4=-0.2
        betaw=-0.5
		etax1=1.4 etax2=0 etaz=-0.4 etax3=-0.1 etax4=-0.5 etam=0 etaw=-0.7
		alpha0=4 alphax1=-1 alphax2=0.3 alphaz=-0.1 alphax3=0.1 alphax4=0.3 alphat=0
        delta1=1 delta2=-0.5 deltam=0
		log_varv=-1.6 log_vare=0.2 log_varu=0.2;
covab=0;
* dur{k} is the duration of nevent in each of the death quantiles, where nevent is a constant in each dur{k};
array dur {4} dur1-dur4;
array baseh {4} log_lam_y1-log_lam_y4;
base_haz_r=exp(log_lam_m1) * event_r1 + exp(log_lam_m2)* event_r2 
+ exp(log_lam_m3) * event_r3 + exp(log_lam_m4) * event_r4;

****need modify this one since here assume the jumppoint for W and M are same;
cum_base_haz_r=exp(log_lam_m1) * dur_r1 * exp(betaw*CD4_dur) + exp(log_lam_m2) * dur_r2 * exp(betaw*CD4_dur)
+ exp(log_lam_m3) *dur_r3*exp(betaw*CD4_dur) + exp(log_lam_m4) * dur_r4 * exp(betaw*CD4_dur);

mu1_1= betax1 * prevoi + betax2 * gender + betaz * trt + betax3 * stratum + betax4 * hemobl + betaw * CD4 + vi + deltam*ui;	/* for recurrent event */
mu1_2= betax1 * prevoi + betax2 * gender + betaz * trt + betax3 * stratum + betax4 * hemobl + vi + deltam*ui;	/* for recurrent event */

loglik1=-exp(mu1_2) * cum_base_haz_r;

sum2=0;
do k=1 to 4 ;
	/* cumulative baseline hazard for time dependent measure */
	sum2=sum2 + exp(baseh{k}) * dur{k} * exp(etam * nevent + etaw * CD4_dur);
	
end;

mu2= etax1 * prevoi + etax2 * gender + etaz * trt + etax3 * stratum + etax4 * hemobl + delta1 * vi + delta2 * ui;	/* for death event */

loglik0=-exp(mu2) * sum2;


if event=2 then do;  /* for death event */
    base_haz_d=exp(log_lam_y1) * event_d1 + exp(log_lam_y2) * event_d2 + exp(log_lam_y3) * event_d3 + exp(log_lam_y4) * event_d4;
	
	mu4= etax1 * prevoi + etax2 * gender + etaz * trt + etax3 * stratum + etax4 * hemobl + etam * nevent  + etaw * CD4 + delta1 * vi + delta2 * ui;	
   
end;
if event=4 then do;
elpson = CD4 - alpha0 - alphax1 * prevoi - alphax2 * gender - alphaz * trt - alphax3 * stratum - alphax4 * hemobl - alphat * stoptime - ui;
loglik2=-0.5*log_vare-0.5*(elpson*elpson/exp(log_vare));
end;

if event=1 then loglik= log(base_haz_r) + mu1_1 + loglik0 + loglik1; 			/*log likelihood for recurrent event */
if event=2 then loglik= loglik0 + log(base_haz_d) + mu4 + loglik1;	/*log likelihood for death */
if event=3 then loglik= loglik0 + loglik1;							/*log likelihood for censoring */
if event=4 then loglik= loglik0 + loglik1 + loglik2;

model id ~ general(loglik);

random vi ui ~ normal([0,0],  [exp(log_varv),covab,exp(log_varu)]) subject=id;

ods output ParameterEstimates=est1 FitStatistics=fit1 CorrMatParmEst=corr1 CovMatParmEst=cov1; 
run;





data _null_;
	set cov1;
	file "&out.\partIII_covI_real_cd4100_quantile.txt";
	put parameter 
		log_lam_m1 log_lam_m2 log_lam_m3 log_lam_m4
		log_lam_y1 log_lam_y2 log_lam_y3 log_lam_y4
		betax1 betax2 betaz betax3 betax4
        betaw
		etax1 etax2 etaz etax3 etax4 etam etaw
		alpha0 alphax1 alphax2 alphaz alphax3 alphax4 alphat
        delta1 delta2 deltam
		log_varv log_vare log_varu;

	run;

data _null_;
	set est1;
	file "&out.\partIII_estI_real_cd4100_quantile.txt";
	put Parameter Estimate;
	run;

ODS RESULTS OFF;
ODS LISTING CLOSE;
ODS EXCEL file="&out.\PartIII estimation setting I _cd4100_quantile.xlsx"
    options (start_at="B1" tab_color="red" absolute_row_height="15" embedded_titles="yes" embedded_footnotes="yes");
	proc report data=est1 nowindows style(Header)=[background=white foreground=black];
	ods Excel;
	column _all_;

	title j=L h=2.5 "Estimation setting I_quantile";

	run;

ODS Excel CLOSE;
ODS LISTING;
ODS RESULTS ON;


ODS RESULTS OFF;
ODS LISTING CLOSE;
ODS EXCEL file="&out.\Part III covariance setting I _cd4100_quantile.xlsx"
    options (start_at="B1" tab_color="red" absolute_row_height="15" embedded_titles="yes" embedded_footnotes="yes");
	proc report data=cov1 nowindows style(Header)=[background=white foreground=black];
	ods Excel;
	column _all_;

	title j=L h=2.5 "covariance Setting I_quantile";

	run;

ODS Excel CLOSE;
ODS LISTING;
ODS RESULTS ON;

