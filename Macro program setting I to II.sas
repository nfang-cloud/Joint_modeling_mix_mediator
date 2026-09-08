proc datasets lib=work kill; run; quit;
DM 'log;clear;out;clear;odsresults;clear;';
***********************************************************************************
*************************Model I***************************************************;

%macro jointmodel1(Data=,mod=);

%do ii=1 %to 200;
title "Iteration &ii.";
*Import the dataset;
*&InD is the path for import the datasets;
PROC IMPORT OUT= WORK.sim1 
            DATAFILE= "&InD.\&Data.&ii..csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;
*Import the jump points;
data quant_r;
infile "&InD.\quant.txt";
input qr0 qr2 qr4 qr6 qr10 aa;
run;

data quant_d;
infile "&InD.\quant.txt";
input qd0 qd2 qd4 qd6 qd10 aa;
run;

data quant_w;
infile "&InD.\quant_W.txt";
input qw0 qw2 qw4 qw6 qw8 qw10 aa;
run;

data four;
set sim1;
aa=1;
run;

data four2;
merge four quant_r quant_d quant_w;
by aa;
if event=3 then event=4;
if event=0 then event=3;
run;
* Calculate W for each stoptime;
Data four3; set four2; by id;
retain W_new;
if first.id then W_new=W;
else if W^=. then W_new=W;
drop W;
rename W_new=W;
run;

/***************Need create a new variable W_dur for the last W value in the duration************/

proc sort data=four3; by event ID;run;

data four3_2; set four3;
by event id;
retain last_W;
if first.id then do;
	   W_dur=0;
	   
	   last_W=W;	
end;
else do;
	W_dur=last_W;
	last_W=W;
end;

if event^=4 then W_dur=W;
drop last_W;
run;

/*
Data W; set four3; if event=4;keep id W stoptime; run;
proc sort data=W; by id stoptime;run;
proc transpose data=W out=W_t; by id; VAR W; ID stoptime;run;
data W_t;
set W_t;
array a(*) _character_;
do i=1 to dim(a);
if a(i) =" " then a(i) = 0;
end;
drop i _NAME_;
rename _0=W0 _2=W2 _4=W4 _6=W6 _8=W8;
run;

Data four4; merge four3 W_t; by id;run;
*/
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
array quant_r {5} qr0 qr2 qr4 qr6 qr10;
array quant_d {5} qd0 qd2 qd4 qd6 qd10;

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
array quant {5} qd0 qd2 qd4 qd6 qd10;
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


proc nlmixed data=five2 qpoints=5 MAXITER=5000 corr cov;

parms log_lam_m1=-1.6 log_lam_m2=-0.7 log_lam_m3=-0.5 log_lam_m4=-0.8
		log_lam_y1=-3.4 log_lam_y2=-2.5 log_lam_y3=-2.0 log_lam_y4=-1.6
		betax=0.2 betaz=0.35 betaw=0.1
		etaz=0.35 etax=0.15 etam=0.25 etaw=0.2
		alpha0=0 alphaz=0.6 alphax=0.2 alphat=0.2
        delta1=1 delta2=-0.5 deltam=0.5
		log_varv=-0.7 log_vare=0 log_varu=0;
covab=0;
* dur{k} is the duration of nevent in each of the death quantiles, where nevent is a constant in each dur{k};
array dur {4} dur1-dur4;
array baseh {4} log_lam_y1-log_lam_y4;
base_haz_r=exp(log_lam_m1) * event_r1 + exp(log_lam_m2)* event_r2 
+ exp(log_lam_m3) * event_r3 + exp(log_lam_m4) * event_r4;

****need modify this one since here assume the jumppoint for W and M are same;
cum_base_haz_r=exp(log_lam_m1) * dur_r1 * exp(betaw*W_dur) + exp(log_lam_m2) * dur_r2 * exp(betaw*W_dur)
+ exp(log_lam_m3) *dur_r3*exp(betaw*W_dur) + exp(log_lam_m4) * dur_r4 * exp(betaw*W_dur);

mu1_1= betaz * X1 + betax * X2 + betaw * W + vi + deltam*ui;	/* for recurrent event */
mu1_2= betaz * X1 + betax * X2 + vi + deltam*ui;	/* for recurrent event */

loglik1=-exp(mu1_2) * cum_base_haz_r;

sum2=0;
do k=1 to 4 ;
	/* cumulative baseline hazard for time dependent measure */
	sum2=sum2 + exp(baseh{k}) * dur{k} * exp(etam * nevent + etaw * W_dur);
	
end;

mu2= etaz * X1 + etax * X2 + delta1 * vi + delta2 * ui;	/* for death event */

loglik0=-exp(mu2) * sum2;


if event=2 then do;  /* for death event */
    base_haz_d=exp(log_lam_y1) * event_d1 + exp(log_lam_y2) * event_d2 + exp(log_lam_y3) * event_d3 + exp(log_lam_y4) * event_d4;
	
	mu4= etaz * X1 + etax * X2 + etam * nevent  + etaw * W + delta1 * vi + delta2 * ui;	
   
end;
if event=4 then do;
elpson = W - alpha0 - alphaz * X1 - alphax * X2 - alphat * stoptime - ui;
loglik2=-0.5*log_vare-0.5*(elpson*elpson/exp(log_vare));
end;

if event=1 then loglik= log(base_haz_r) + mu1_1 + loglik0 + loglik1; 			/*log likelihood for recurrent event */
if event=2 then loglik= loglik0 + log(base_haz_d) + mu4 + loglik1;	/*log likelihood for death */
if event=3 then loglik= loglik0 + loglik1;							/*log likelihood for censoring */
if event=4 then loglik= loglik0 + loglik1 + loglik2;

model id ~ general(loglik);

random vi ui ~ normal([0,0],  [exp(log_varv),covab,exp(log_varu)]) subject=id;

ods output ParameterEstimates=est&mod&ii FitStatistics=fit&mod&ii CorrMatParmEst=corr&mod&ii CovMatParmEst=cov&mod&ii; 
run;

%end;

%mend;

%jointmodel1(Data=part3_sim1_alldata,mod=1);




/* &out is the saving file path*/
*Output the covariance;
%macro outcov(mod=,filename=);
%do ii=1 %to 200;
data _null_;
	set cov&mod&ii;
	file "&out.\&filename&ii..txt";
	put  parameter 
        log_lam_m1 log_lam_m2 log_lam_m3 log_lam_m4
		log_lam_y1 log_lam_y2 log_lam_y3 log_lam_y4
		betax betaz betaw
		etaz etax etam etaw
		alpha0 alphaz alphax alphat
        delta1 delta2 deltam
		log_varv log_vare log_varu;
	run;
%end;
%mend;

*Output the estimation;
%macro outest(mod=,filename=);
%do ii=1 %to 200;
data _null_;
	set est&mod&ii;
	file "&out.\&filename&ii..txt";
	put Parameter Estimate;
	run;
%end;
%mend;

%outest(mod=1,filename=estI);
%outcov(mod=1,filename=covI);

***************Simulation 2********************;
*reset input and output path as needed;

%jointmodel1(Data=part3_sim2_alldata,mod=2);


/* &out is the saving file path*/

%macro outcov(mod=,filename=);
%do ii=1 %to 200;
data _null_;
	set cov&mod&ii;
	file "&out.\&filename&ii..txt";
	put parameter
log_lam_m1 log_lam_m2 log_lam_m3 log_lam_m4
		log_lam_y1 log_lam_y2 log_lam_y3 log_lam_y4
		betax betaz betaw
		etaz etax etam etaw
		alpha0 alphaz alphax alphat
        delta1 delta2 deltam
		log_varv log_vare log_varu;
	run;
%end;
%mend;

*Output the estimation;
%macro outest(mod=,filename=);
%do ii=1 %to 200;
data _null_;
	set est&mod&ii;
	file "&out.\&filename&ii..txt";
	put Parameter Estimate;
	run;
%end;
%mend;

%outest(mod=2,filename=estII);
%outcov(mod=2,filename=covII);

