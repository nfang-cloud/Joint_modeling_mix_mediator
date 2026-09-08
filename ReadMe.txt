**************************************************************************************************************************
*****************************************Programs for the Study*********************************************************

*******************Part I Simulation Setting I to Setting II*************************
1. Generate Datasets for setting I to II. R
    - Used to generate 200 datasets for each setting

2. Data Generation and Calculate the True Value for NDE NIE.R
   -Used to generate 1000 datasets to compute the true value for NDE NIE_m NIE_w
   

3. Macro program setting I to II. sas
   - Estimate the model parameters under settings I to II

4. Summary the fitting results from setting I to II.R
   -Summary the estimation of model parameters under settings I to II
   - Resutls shown in Table S1-S2



*******************Part II Self Defined Function*************************
1. Lambdainv.R
    -Function used to find vector T such that Lambda(T)=s for a vector input s

2. myprod.R
    -Fucntion used to product of two step function

3. mysimRec.R
   -Function using Cinlar's inversion Method to generate Non-homogeneous Poisson process

4. datagenT01_M_CD4.R
   - Function used to generate recurrent events

5. datagen_X.R
   -Function used to generate datasets with simulated X, Z and W

6. mymed.R
   -Function used to compute the NDE/NIE_m/NIE_w for simulation settings

7. mymed_T01.R
   -Function used to compute the NDE/NIE_m/NIE_w for plotting

8. S.R
    - Used to compute the survival functions

 

*******************Part III Estimate NDE/NIE Under Simulation Setting I to II*************************

1. Folder "NDE NIE from Simulation I"
    *SimulationI_parallel.R
              - Function used to estimate NDE/NIE_m/NIE_w under simulation setting I
    *SimulationI_parallel_boot.R
             -Function of the bootstrap for NDE/NIE_m/NIE_w under simulation setting I
     *SimulationI_summary.R
             -Function used to summary and estimate bias, SD, Mese and CR for NDE/NIE under simulation setting I


2. Folder "NDE NIE from Simulation II"
    *SimulationII_parallel.R
              - Function used to estimate NDE/NIE_m/NIE_w under simulation setting II
    *SimulationII_parallel_boot.R
             -Function of the bootstrap for NDE/NIE_m/NIE_w under simulation setting II
    *SimulationII_summary.R
             -Function used to summary and estimate bias, SD, Mese and CR for NDE/NIE under simulation setting II





*******************Part IV Real Data Analysis*************************

1. CPCRA study analysis CD4 and recur OI model fitting.sas
    -Estimates the parameters under simulation setting I with five OIs for CPCRA study

2.Folder "CPCRA NDE NIE Plot"

   * NDE_NIE_PrevOI.R 
      -Estimate NDE/NIE_OI/NIE_CD4 for prevOI

   * NDE_NIE_CI_PrevOI.R
     -Boostrap of NDE/NIE_OI/NIE_CD4 for prevOI

  * NDE_NIE_TRT.R 
      -Estimate NDE/NIE_OI/NIE_CD4 for treatment

   * NDE_NIE_CI_TRT.R
     -Boostrap of NDE/NIE_OI/NIE_CD4 for treatment

   * Summary PrevOI and Treatment.R
    -Summary the results and plot the estimates of NDE/NIE_OI/NIE_CD4 with bootstraped 95% CI for prevOI and Treatment
