# Acknowledgements
This work was done in collaboration with Public Health Seattle & King County as part of a CSTE/CDC - supported initiative, "Development of forecast, analytic, and visualization tools to improve outbreak response and support public health decision making."

This model is an adaptation of earlier work, please see:
Pitzer VE, Viboud C, Alonso WJ, et al. "Environmental Drivers of the Spatiotemporal Dynamics of Respiratory Syncytial Virus in the United States". PLOS Pathogens. 2015. https://doi.org/10.1371/journal.ppat.1004591

Zheng Z, Weinberger DM, Pitzer VE. "Predicted effectiveness of vaccines and extended half-life monoclonal antibodies against RSV hospitalizations in children". NJP Vaccines. 2022. https://doi.org/10.1038/s41541-022-00550-5

# Model Structure 
 <img src="https://github.com/chelsea-hansen/RSV-Interventions/assets/81387982/c53f4f2a-3a92-4ce3-8204-bafdbb18b74c" width="40%" height="40%" align="left">
 
The model assumes that all infants are born into an "M" compartment (representing maternally derived immunity) with partial immunity against infection and hospitalization given infection. After this protection wanes, infants become fully susceptible (S0). Following the first infection (I1) individuals have a short period of immunity from infection (R1). After this immunity wanes, individuals are susceptible again, but with a lower relative risk of infection. Following each infection the duration of infectiousness becomes shorter, the duration of immunity increases, and the relative risk of future infections becomes less. As mentioned above, this model is an adaptation of earlier work. The model has been modified to include a recovered "R" compartment following each infectious "I" compartment. A list of parameters is provided below. Parameters marked with * have been adopted from the earlier work by (Pitzer et al, 2015).

|Parameter|Fixed Value|
|---|---|
|*Duration of infectiousness - first infection (1/&gamma;<sub>1</sub>)|10 days|
|*Duration of infectiousness - second infection (1/&gamma;<sub>2</sub>)|7 days|
|*Duration of infectiousness - third or later infection (&gamma;<sub>3</sub>)|5 days|
|*Relative risk of infection following first infection (&sigma;<sub>1</sub>)|0.76|
|*Relative risk of infection following second infection (&sigma;<sub>2</sub>)|0.6|
|*Relative risk of infection following third or later infection (1/&sigma;<sub>3</sub>)|0.4|
|Relative risk of infection with maternal immunity (1/&sigma;<sub>3</sub>)|0.4|
|Relative risk of hospitalization given infection for M compartment|0.7|
|*Duration of maternal immunity (1/&omega;<sub>1</sub>)|112 days|
|Duration of immunity following first and second infections (1/&omega;<sub>2</sub>)|182.625 days|
|Duration of immunity following third or later infections (1/&omega;<sub>3</sub>)|365.25 days|
|*Relative infectiousness - second infections (&rho;<sub>1</sub>)|0.75|
|*Relative infectiousness - third or later infections (&rho;<sub>2</sub>)|0.51|
|Baseline transmission rate (&beta;)|Fitted|
|Amplitude of seasonal forcing (*b*1)|Fitted|
|Phase of seasonal forcing (&phi;)|Fitted|
|Infections that lead to reported hospitalizations (&theta;)|Fitted|


# Step 1 - Data 
The data needed to run the model can be found in the ```1. Data``` folder. Example datasets from King County, Washington are provided. Please note, in the sample datasets values between 1-9 have been supressed and reinterpolated. The data folder is further divided into 2 subfolders: ```RSV Data``` and ```Demographic Data```. Details for each subfolder are provided below

## RSV Data
To run the model you will need to have a weekly time series of RSV hospitalizations (or ED visits) and an age distribution of RSV hospitalizations/ED visits. For the weekly time series it is best if you can have at least 3 years prior to the COVID-19 pandemic, however the code should work with a slightly shorter time series. The sample dataset is from January 2017 - November 2023. An example is provided below. Note, a 3-week moving average has been applied to the time series, and values have been rounded to the nearest whole number. The model fitting procedure uses a Poisson regression and you will get an error if the RSV time series has not been rounded to whole numbers. 

The age distribution is divided into 2 time periods: pre-pandemic (Janaury 2017 - March 2020) and post-pandemic (April 2020 - November 2023). The example uses 5 age groups (<6 months, 6-11 months, 1-4 years, 5-59 years, 60+ years). See below. The code could be modified to use different age groups. The code can also be run without the age distribution, however if this is done the projections based on the intervention scenarios will only be for all ages, not age-specific estimates. 

## Demographic data
The code will also require birth rates and the age-specific population distribution. The ```data_prep.R``` R script will pull and format all of the necessary data. 

The first dataset you will create is ```yinit.rds```. This dataset divides the 13 age groups (<2m, 2-3m, 4-5m, 6-7, 8-9, 10-11m, 1y, 2-4y, 5-9y, 10-19y, 20-39y, 40-59y, 60+y) into the model compartments, starting with the M and S0 compartments. Note: these age groups do not need to be the same as the age groups from the RSV age distributions. One infection is seeded into each age group >6m in the I1 compartment. See below. 

The model will initiate in January 1995 and "burn-in" until your RSV time series begins (in the sample this is January 2017). 

The code will also save another format of this dataset ```yinit.vector.rds``` and versions with the additional compartments for the immunizations (Mn,Mv,N,Si,Vs1,Vs2), ```yinit_interventions.rds``` and ```yinit.vector_interventions.rds```. These versions will be used during the Interventions step later. 

Birth rate data has already been pulled from CDC Wonder and saved as ```birth_rates_by_state.rds``` and ```birth_rates_by_county.rds```. Use the provided code to format correctly for the model. 

The ```data_prep.R``` script will also save a few additional parameters as ```other_parms.rds```.

To run the next step (calibration you will need 7 datasets). 
1. The RSV time series (example saved as ```rsv_ts.rds```)
2. The RSV age distribution (example saved as ```age_distribution.rds```)
3. ```yinit.rds``` - created with the ```data_prep.R``` R script 
4. ```yinit.vector.rds``` - created with the ```data_prep.R``` R script
5. ```other_parms.rds``` - created with the ```data_prep.R``` R script
6. The birth rate dataset created with the ```data_prep.R``` R script (example provided as ```births_kingcounty.rds```
7. ```contact_POLYMOD.rds``` - already provided in the correct format.
Please make sure you have prepared all of these datasets before moving to the next step. 
   
# Step 2 - Calibration
This step is the trickiest step in the process and likely where you will have the most difficulty. But the good news is you only need to run this once and then you can update change your immunozation scenarios as you like! 

This step uses the R scripts in the ```2. Calibration``` folder. ```model_dynamics.R``` is the function which runs the model equations. ```model_calibration.R``` is the script which uploads the data, calls the rsv_dynamics function, and uses maximum likelihood estimation to estimate the model parameters. You will use the ```model_calibration.R``` script. As shown in the table at the top of the page, many of the parameters are fixed based on earlier papers. 

The initial step uses Maximil Likelihood Estimation to fit 4 parameters. The initial step is fitting only to the pre-pandemic time period (before March 2020). 
 1.  The baseline transmission rate (&beta;), bounded between 6 and 9.
 2.  The amplitude of seasonal forcing (*b*1).
 3.  The timing/phase of seasonal forcing (&phi;), bounded between 0 and 2&pi;.
 4.  The proportion of infections leading to reported hospitalizations - reporting rate (&theta;), bounded to be between 0 and 1. 

The second step uses these fitted parameters to fit the pandemic through the 2022-23 rebound season. To fit this time period the code will reduce the number of contacts by different amoungs at different points. In addition to reducing the number of contacts, exteral seeing of infections is reduced to zero from April 2020 - February 2021, and then gradually returns to normal by May 2021. 
Fitting contact reductions: 
 1. It is assumed that there is a reduction in contacts from April - June 2020, corresponding with stay-at-home-order. Model will fit the amount of reduction (proportion between 0-1, with 1 representing pre-pandemic contact patterns)
 2. It is assumed that following the end of stay-at-home orders, contacts remained lower into 2021. The model will fit the amount of reduction and the number of weeks contacts remained supressed.
 3. It is assumed that contacts began to return to normal in 2021. The model will estimate how long this ramp-up took and what level it reached by the end of 2021.
 4. It is assumed that there was another drop in contacts coinciding with the emergence of the Omicron variant in winter 2021-2022. The model will fit the amount and duration of this reduction.
 5. It is assumed that following the initial Omicron wave, contacts gradually increased through 2022, returning to pre-pandemic levels by the fall of 2022. The model will fit the number of weeks this return to normal took.

In addition to fitting the reductions in contacts, the model will fit a new reporting rate to account for increases in testing and RSV case reporting. The new reporting rate is bounded between the pre-pandemic reporting rate, and a 100% increase. 

After completing parts 1 and 2 of the calibration step, you will have a figure that looks like this: 


# Step 3 - Interventions 
This step uses the R scripts in the ```3. Interventions``` folder. ```MSIRS_scenarios.R``` is the script that runs the scenarios, ```MSIRS_intervention_models.R``` is teh script which runs the model equations and ```MSIRS_scenario_functions.R``` is an R script that interfaces between the other 2 scripts. You will open and run the ```MSIRS_scenarios.R``` script. The script is designed to run 8 scenarios with a combination of optimistic and pessimistic coverage and effectiveness for the three interventions (RSV vaccination for adults >60 years, RSV vaccination for pregnant women, RSV monoclonal antibodies for infants <8 months). See details in the table below. Many of the parameters for the vaccinations have been fixed based on clinical trials (see tables). In the script the user is able to easily set an optimistic and pessimistic scenario for cumulative coverage of the intervention. 

Scenario Overview: 
|Senior Immunization|Infant Immunization||||
|-|-|-|-|-|
|---|Monoclonal Antibodies - Optimistic|Monoclonal Antibodies - Optimistic|Monoclonal Antibodies - Pessimistic|Monoclonal Antibodies - Pessimistic|
|---|Maternal Vaccination - Optimistic|Maternal Vaccination - Pessimistic|Maternal Vaccination - Optimistic|Maternal Vaccination - Pessimistic|
|Senior Vaccination - Optimistic|A|B|C|D|
|Senior Vaccination - Pessimistic|E|F|G|H|






This folder also includes a file ```coverage_curves_2023_24.rds``` which includes coverage curves for influenza, scaled between 0 and 1. 


function that runs the model equations (with compartments for vaccination in seniors and monoclonal antibodies in infants). intervention_scenarios.R uses the parameters from Step 2 and estimates the impact of interventions under different scenarios for intervention coverage, effectiveness, and timing. Sample scenarios are provided, but these can be modified to answer your research questions. 

### Notes on interventions 
The model assumes that interventions are providing protection against severe disease (hospitalization) but not against infection. The model also assumes that the interventions do not impact an individual's infectiousness if they become infected. 

### Notes on fitting to the pandemic period 
The interventions model uses the model parameters fit during the pre-pandemic period and then projects forward. To account for distruptions during the COVID-19 pandemic we reduced the baseline transmission rate (beta) by 25% from April - June 2020 (representing the stay-at-home order) and then allowed beta to increase linearly to pre-pandemic levels between July 2020 and March 2021. We included an additional 25% beta reduction from mid-December 2021 through February 2022 to account for the emergence of the Omicron BA.1 variant. We did this because it is likely that people adopted protective measures during this time, reducing RSV transmission. Furthermore, it is also possible that viral interactions may have resulted in a decline in RSV transmission during this period. Together with these reductions in beta, we reduced imported infections (seeding) to zero from April 2020 - January 2021 and allowed for a linear increase to normal levels between February 2021 - May 2021. This approach is based on the paper by Zheng et al. "Estimation of the Timing and Intensity of Reemergence of Respiratory Syncytial Virus Following the COVID-19 Pandemic in the US." JAMA Netw Open. 2021. 
