# Acknowledgements
This work was done in collaboration with Public Health Seattle & King County as part of a CSTE/CDC - supported initiative, "Development of forecast, analytic, and visualization tools to improve outbreak response and support public health decision making."

This model is an adaptation of earlier work, please see:
Pitzer VE, Viboud C, Alonso WJ, et al. "Environmental Drivers of the Spatiotemporal Dynamics of Respiratory Syncytial Virus in the United States". PLOS Pathogens. 2015. https://doi.org/10.1371/journal.ppat.1004591

Zheng Z, Weinberger DM, Pitzer VE. "Predicted effectiveness of vaccines and extended half-life monoclonal antibodies against RSV hospitalizations in children". NJP Vaccines. 2022. https://doi.org/10.1038/s41541-022-00550-5

# Model Structure 
 <img src="https://github.com/chelsea-hansen/RSV-Interventions/assets/81387982/c53f4f2a-3a92-4ce3-8204-bafdbb18b74c" width="35%" height="40%" align="left">
 
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


### birth 
This dataset is the birth rate per 1000 population for each week you will be running the model (from initiation of burn-in to projections). You can use the annual birth rate for each week or set the first week of the year to the annual birth rate and interpolate missing weeks between years. Later in the code this will be divided by 52.17 to convert from an annual to a weekly scale. Note, there is a column for each age group but only the first column (infants <2m) has data. The other columns are set to zero. 
### contact 
This model uses the contact martix from the POLYMOD study (Mossong et al. 2008. "Social Contacts and Mixing Patterns Relevant to the Spread of Infectious Diseases." PLOS Medicine. https://doi.org/10.1371/journal.pmed.0050074.) The contact matrix has been aggregated to match the age groups used in this model. 
### rsv_ts
This is a weekly count of RSV-coded hospitalizations (and optionally, ED visits) for all age groups combined. Values from 1-9 in the sample dataset have been suppressed and are interpolated. It is recommended to have at least 3 years of data prior to the COVID-19 pandemic, but more years of data is better. 
### age_distribution
This is the proporiton of RSV-coded hospitalizations (and optionally, ED visits) in each age group (with middle age groups collapsed in to 1-4yrs and 5-64yrs). The proportions are separated into the pre-pandemic time period (pre-March 2020) and the post-pandemic (post-March 2020) time period. 
### scale_ed_to_hosp (optional)
The model is calibrated to RSV-coded hospitalizations. This dataset is used to scale hospitalizations to ED visits using the ratio of hospitalizations to ED visits in each age group during the pre- and post- pandemic time periods. This dataset is only necessary if you would like to produce estimates for ED visits in addition to hospitalizations. 

# Step 2 - Parameter Estimation 
This step uses the R scripts in the "Calibration" folder. rsv_dynamics.R is the function which runs the model equations. model_calibration.R is the script which uploads the data, calls the rsv_dynamics function, and uses maximum likelihood estimation to estimate the model parameters. You will use the model_calibration.R script. Many parameters are set based on the literature. 
The model fits 8 parameters: 
 - The baseline transmission rate (beta)
 - The amplitude of seasonal forcing
 - The timing/phase of seasonal forcing
 - The duration of transplacental immunity
 - The relative risk of hospitalization given infection among infants protected by transplacental immunity
 - The number of weekly external introductions of infection (seeding)
 - The proportion of infections leading to a reported RSV hospitalization in seniors
   - The proportion of infections leading to hospitalizations in the 1-4yrs age group is set to the same value as the proportion in adults 65+ yrs (van Boven et al. "Estimating Transmission Parameters for Respiratory Syncytial Virus and Predicting the Impact of Maternal and Pediatric Vaccination." JID. 2020.)
   - The proportion of infections leading to hospitalizations in the 5-64 yrs age group is .06 times the proportion in adults 65+ yrs, based on influenza (Biggerstaff et al. "Estimating the Potential Effects of a Vaccine Program Against an Emerging Influenza Pandemic—United States." CID. 2015.)
 - The proportion of infections leading to a reported RSV hospitalization in infants <2 months
   - The proportion of infections in infants >2 months (in 2-month age bands) is set relative to those <2 months based on the paper by Zheng et al referenced at the top of this document.
  
Save the estimated parameters for use in modeling the impact of interventions

# Step 3 - Model the Impact of Interventions 
This step uses the R scripts in the "Interventions" folder. intervention_models.R is the function that runs the model equations (with compartments for vaccination in seniors and monoclonal antibodies in infants). intervention_scenarios.R uses the parameters from Step 2 and estimates the impact of interventions under different scenarios for intervention coverage, effectiveness, and timing. Sample scenarios are provided, but these can be modified to answer your research questions. 

### Notes on interventions 
The model assumes that interventions are providing protection against severe disease (hospitalization) but not against infection. The model also assumes that the interventions do not impact an individual's infectiousness if they become infected. 

### Notes on fitting to the pandemic period 
The interventions model uses the model parameters fit during the pre-pandemic period and then projects forward. To account for distruptions during the COVID-19 pandemic we reduced the baseline transmission rate (beta) by 25% from April - June 2020 (representing the stay-at-home order) and then allowed beta to increase linearly to pre-pandemic levels between July 2020 and March 2021. We included an additional 25% beta reduction from mid-December 2021 through February 2022 to account for the emergence of the Omicron BA.1 variant. We did this because it is likely that people adopted protective measures during this time, reducing RSV transmission. Furthermore, it is also possible that viral interactions may have resulted in a decline in RSV transmission during this period. Together with these reductions in beta, we reduced imported infections (seeding) to zero from April 2020 - January 2021 and allowed for a linear increase to normal levels between February 2021 - May 2021. This approach is based on the paper by Zheng et al. "Estimation of the Timing and Intensity of Reemergence of Respiratory Syncytial Virus Following the COVID-19 Pandemic in the US." JAMA Netw Open. 2021. 
