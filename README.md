
**We will be hosting a webinar through CSTE on September, 2024 12-1pm Eastern Time to explain the model and code posted here. Please register using this link:https://cste-org.zoom.us/webinar/register/WN_ZFPLhr8qQOOTkaiYE8XUOQ#/registration**

In 2023 several new immunizations for RSV were approved. These include 2 vaccines for adults >=60 years, a vaccination for pregnant women (to protect newborns), and an extended half-life monoclonal antibody for infants <8 months. The R scripts and datasets provided here fit a deterministic MSIRS model to RSV hospitalizations and provide projections for the impact of these new interventions under optimistic and pessimistic scenarios for coverage and effectiveness. 

For more information about recommendations for these new immunizations please see: https://www.cdc.gov/vaccines/vpd/rsv/index.html


# Acknowledgements
This work was done in collaboration with Public Health Seattle & King County as part of a CSTE/CDC - supported initiative, "Development of forecast, analytic, and visualization tools to improve outbreak response and support public health decision making."

This model is an adaptation of earlier work, please see: Pitzer et al. and Zheng et al.

# Model Structure 
 <img src="https://github.com/chelsea-hansen/RSV-Interventions/assets/81387982/c53f4f2a-3a92-4ce3-8204-bafdbb18b74c" width="35%" height="35%" align="left">
 
The model assumes that all infants are born into an "M" compartment (representing maternally derived immunity) with partial immunity against infection. After this protection wanes, infants become susceptible (S0). Following the first infection (I1) individuals have a short period of immunity from infection (R1). After this immunity wanes, individuals are susceptible again, but with a lower relative risk of infection (S1). Following each infection the duration of infectiousness becomes shorter, the duration of immunity increases, and the relative risk of future infections is lower. As mentioned above, this model is an adaptation of earlier work. The model has been modified to include a recovered "R" compartment following each infectious "I" compartment. A list of parameters is provided below. 

|Parameter|Fixed Value|
|---|---|
|<sup>1</sup>Duration of infectiousness - first infection (1/&gamma;<sub>1</sub>)|10 days|
|<sup>1</sup>Duration of infectiousness - second infection (1/&gamma;<sub>2</sub>)|7 days|
|<sup>1</sup>Duration of infectiousness - third or later infection (1/&gamma;<sub>3</sub>)|5 days|
|<sup>3</sup>Relative risk of infection following first infection (&sigma;<sub>1</sub>)|0.89|
|<sup>3</sup>Relative risk of infection following second infection (&sigma;<sub>2</sub>)|0.72|
|<sup>3</sup>Relative risk of infection following third or later infection (&sigma;<sub>3</sub>)|0.24|
|Relative risk of infection with maternal immunity (same as RR following third infection) (&sigma;<sub>3</sub>)|0.24|
|<sup>1</sup>Duration of maternal immunity (1/&omega;<sub>1</sub>)|112 days|
|<sup>4</sup>Duration of immunity following first and second infections (1/&omega;<sub>2</sub>)|182.625 days|
|<sup>3</sup>Duration of immunity following third or later infections (1/&omega;<sub>3</sub>)|358.9 days|
|<sup>1</sup>Relative infectiousness - second infections (&rho;<sub>1</sub>)|0.75|
|<sup>1</sup>Relative infectiousness - third or later infections (&rho;<sub>2</sub>)|0.51|
|Baseline transmission rate (&beta;)|Fitted|
|Amplitude of seasonal forcing (*b*1)|Fitted|
|Phase of seasonal forcing (&phi;)|Fitted|
|Infections that lead to reported hospitalizations (<2m, 2-11 months fixed relative to this)|Fitted|
|Infections that lead to reported hospitalizations (1-4yrs)|Fitted|
|Infections that lead to reported hospitalizations (5-59 yrs)|Fitted|
|Infections that lead to reported hospitalizations (60+ yrs)|Fitted|
|Proportion of contacts relative to the pre-pandemic period, April 2020 - June 2020 |Fitted|
|Proportion of contacts relative to the pre-pandemic period, July 2020 - March 2021|Fitted|
|Proportion of contacts relative to the pre-pandemic period, October 2021 - December 2021|Fitted|
|Proportion of contacts relative to the pre-pandemic period, January 2022 - June 2022|Fitted|

# Step 1 - Data 
The data needed to run the model can be found in the ```1. Data``` folder. The data folder is further divided into 2 subfolders: ```RSV Data``` and ```Demographic Data```. Details for each subfolder are provided below

## RSV Data
To run the model you will need to have a weekly (or monthly) time series of RSV hospitalizations (or ED visits) and an age distribution of RSV hospitalizations (or ED visits if that is the time series you are using). Example datasets from King County, Washington are provided. **Please note, in the sample datasets values between 1-9 have been supressed and reinterpolated. Additionally, these sample data have been provided here for educational purposes, but are not intended to be redistributed or used for other purposes**. For the weekly time series it is best if you can have at least 3 years of data prior to the COVID-19 pandemic, however the code should work with a slightly shorter time series. The sample dataset is from January 2017 - November 2023. Note, a 3-week moving average has been applied to the time series, and values have been rounded to the nearest whole number. The model fitting procedure uses a Poisson distribution and you will get an error if the RSV time series has not been rounded to whole numbers. 

The age distribution is divided into 2 time periods: pre-pandemic (Janaury 2017 - March 2020) and post-pandemic (April 2020 - November 2023). The example uses 5 age groups (<6 months, 6-11 months, 1-4 years, 5-59 years, 60+ years). The code could be modified to use different age groups.


## Demographic data
The code will also require birth rates and the age-specific population distribution. The ```data_prep.R``` R script will pull and format all of the necessary data using the ```tidycensus``` R package. 

The first dataset you will create is ```yinit.rds```. This dataset divides the population from 13 age groups (<2m, 2-3m, 4-5m, 6-7, 8-9, 10-11m, 1y, 2-4y, 5-9y, 10-19y, 20-39y, 40-59y, 60+y) into the model compartments, starting with the M and S0 compartments. Note: these age groups do not need to be the same as the age groups from the RSV age distributions. One 
infection is seeded into each age group >6m in the I1 compartment. See below. The model will initiate in January 1995 and "burn-in" until your RSV time series begins (in the sample this is January 2017). 

<img src="https://github.com/chelsea-hansen/RSV-Interventions/assets/81387982/4f5419ff-a5d0-464b-b241-b7c7e32824e8" width="80%" height="80%" align="center">

The code will also save another format of this dataset ```yinit.vector.rds``` and versions with the additional compartments for the immunizations (Mn,Mv,N,Si,Vs1,Vs2), ```yinit_interventions.rds``` and ```yinit.vector_interventions.rds```. These versions will be used during the Interventions step later. 

Birth rate data has already been pulled from CDC Wonder and saved as ```birth_rates_by_state.rds``` and ```birth_rates_by_county.rds```. Use the provided code to format correctly for the model. 

The ```data_prep.R``` script will also save a few additional parameters as ```other_parms.rds```.

To run the next step (Calibration) you will need 7 datasets. 
1. The RSV time series (example saved as ```rsv_ts.rds```)
2. The RSV age distribution (example saved as ```age_distribution.rds```)
3. ```yinit.rds``` - created with the ```data_prep.R``` R script 
4. ```yinit.vector.rds``` - created with the ```data_prep.R``` R script
5. ```other_parms.rds``` - created with the ```data_prep.R``` R script
6. The birth rate dataset created with the ```data_prep.R``` R script (example provided as ```births_kingcounty.rds```
7. ```contact_POLYMOD.rds``` - already provided in the correct format.

Please make sure you have prepared all of these datasets before moving to the next step. 
   
# Step 2 - Calibration
Part 1: This step is the trickiest step in the process and likely where you will have the most difficulty. But the good news is you only need to run this once and then you can update and change your immunization scenarios as much as you like! 

This step uses the R scripts in the ```2. Calibration``` folder. ```model_dynamics.R``` is the function which runs the model equations. ```model_calibration.R``` is the script which uploads the data, calls the model_dynamics function, and uses Maximum Likelihood Estimation to estimate the model parameters. You will use the ```model_calibration.R``` script. You do not need to open the ```model_dynamics.R``` script to run the code. As shown in the table at the top of the page, many of the parameters are fixed based on earlier papers. 

The initial step uses Maximum Likelihood Estimation to fit 4 parameters. The initial step is fitting only to the pre-pandemic time period (before March 2020). 
 1.  The baseline transmission rate (&beta;), bounded between 6 and 9.
 2.  The amplitude of seasonal forcing (*b*1).
 3.  The timing/phase of seasonal forcing (&phi;), bounded between 0 and 2&pi;.
 4.  The proportion of infections leading to reported hospitalizations - reporting rate (&theta;), bounded to be between 0 and 1. 

Part 2: The second part of this step uses these fitted parameters to fit data during the pandemic through the 2022-23 rebound season. To fit this time period the code will reduce the number of contacts by different amounts at different time points. In addition to reducing the number of contacts, external seeing of infections is reduced to zero from April 2020 - February 2021, and then gradually returns to normal by May 2021. 

Assumptions when fitting contact reductions: 
 1. It is assumed that there is a reduction in contacts from April - June 2020, corresponding with stay-at-home-orders. The model will fit the magnitude of the reduction (proportion between 0-1, with 1 representing pre-pandemic contact patterns)
 2. It is assumed that following the end of stay-at-home orders, contacts remained lower into 2021. The model will fit the magnitude of the reduction and the number of weeks contacts remained supressed.
 3. It is assumed that contacts began to return to normal in 2021. The model will estimate the level of contacts by the end of 2021.
 4. It is assumed that there was another drop in contacts coinciding with the emergence of the Omicron variant in winter 2021-2022. The model will fit the magnitude and duration of this reduction.
 5. It is assumed that following the initial Omicron wave, contacts gradually increased through 2022, returning to pre-pandemic levels by the fall of 2022. The model will fit the number of weeks over which this return to normal occurs.

In addition to fitting the reductions in contacts, the model will fit a new reporting rate to account for increases in testing and RSV case reporting. The new reporting rate is bounded between the pre-pandemic reporting rate, and a 100% increase. 

After completing parts 1 and 2 of the calibration step, you will have a figure that looks like this: 

<img src="https://github.com/chelsea-hansen/RSV-Interventions/assets/81387982/08c3d85f-78ec-4203-a02a-68013e36ae75" width="80%" height="80%" align="center">

The last part of the ```Calibration.R``` script uses Latin Hypercube Sampling to sample parameter values from a plausible range around the fitted parameters (+/- 3%). Two versions are saved, one with 100 samples and the other with 1000 samples. These sets of parameter values will be used for estimating confidence intervals around the model projections. 

Note: The saved parameters provided here were estimated using a burn-in period starting in 1980. The code has been updated to start in 1995 to match the publicly available demographic data for other states and counties provided in the ```1. Data``` folder. Therefore, when working through this example the estimated parameters you get might be slightly different from the ones saved here. 

# Step 3 - Interventions 
This step uses the R scripts in the ```3. Interventions``` folder. ```MSIRS_scenarios.R``` is the script that runs the scenarios, ```MSIRS_intervention_models.R``` is the script which runs the model equations and ```MSIRS_scenario_functions.R``` is an R script that interfaces between the other 2 scripts. You will open and run the ```MSIRS_scenarios.R``` script. The script is designed to run 9 scenarios: 1 counterfactual (no interventions) plus 8 combinations of optimistic and pessimistic coverage and effectiveness for the three interventions (RSV vaccination for adults >60 years, RSV vaccination for pregnant women, RSV monoclonal antibodies for infants <8 months). See details in the table below. Many of the parameters for the immunizations have been fixed based on clinical trials (see tables). This folder also includes a file ```coverage_curves_2023_24.rds``` which includes coverage curves scaled between 0 and 1. In the script the user is able to easily set an optimistic and pessimistic scenario for cumulative coverage of the intervention using these curves. 

For each scenario, the counterfactual and A-G, there is an option to run just the point estimates (very fast) and an option to run the confidence intervals using the parameters sampled with Latin Hypercube Sampling. Running the confidence intervals with 1000 replicates is recommended, but will take around 2 hours for each scenario. Running the confidence intervals with 100 replicates will take 15-30 minutes per scenario. 

Scenario Overview: 
|Senior Immunization|Infant Immunization||||
|-|-|-|-|-|
||Monoclonal Antibodies - Optimistic|Monoclonal Antibodies - Optimistic|Monoclonal Antibodies - Pessimistic|Monoclonal Antibodies - Pessimistic|
||Maternal Vaccination - Optimistic|Maternal Vaccination - Pessimistic|Maternal Vaccination - Optimistic|Maternal Vaccination - Pessimistic|
|Senior Vaccination - Optimistic|A|B|C|D|
|Senior Vaccination - Pessimistic|E|F|G|H|

## Infant Immunizations
The model now assumes that a proportion of infants are born to vaccinated mothers (Mv Compartment) and a proportion of infants receive monoclonal antibodies at or shortly after birth (Mn Compartment). These infants retain the same protection against infection as the M compartment, but have a higher protection against hospitalization given infection. Additionally, some infants receive a monoclonal antibody a few months after birth (N Compartment). These infants do not have any protection against infection, but they do have protection against hospitalization given infection. When the protection wanes from these compartments (Mv, Mn, N), infants move to the Si compartment. This compartment is functionally the same as the S0 compartment, but ensures that infants do not receive two interventions. 

<img src="https://github.com/chelsea-hansen/RSV-Interventions/assets/81387982/8c830c59-3d7d-4100-8a7d-f4f728656f62" width="50%" height="50%" align="left">


|Parameter|Optimistic Value|Pessimistic Value|
|---------|----------------|-----------------|
|Duration of protection from monoclonal antibodies (days)|180|180|
|Effectiveness of monoclonal antibodies against hospitalization|80%|60%|
|Cumulative coverage of monoclonal antibodies|user defined|user defined|
|Duration or protection from maternal vaccination (days)|180|180|
|Effectiveness of maternal vaccination against hospitalization|70%|50%|
|Cumulative coverage of maternal vaccination|user defined|user defined| 

Clinical trials

Monoclonal antibodies: (Hammitt et al, 2022) https://www.nejm.org/doi/full/10.1056/NEJMoa2110275

Maternal Vaccination: (Kampmann et al, 2023) https://www.nejm.org/doi/full/10.1056/NEJMoa2216480

## Senior Vaccination 
The vaccination compartment (Vs1) draws seniors from the S3 and R4 compartments. Current data suggests that the vaccine is effective for at least 2 seasons. Seniors spend approximately 1 year in the Vs1 compartment before waning to the Vs2 compartment for another year and then returning to the S3 compartment. 

<img src="https://github.com/chelsea-hansen/RSV-Interventions/assets/81387982/a569dc65-8398-4059-8782-917b27c42041" width="50%" height="50%" align="left">


|Parameter|Optimistic Value|Pessimistic Value|
|---------|----------------|-----------------|
|Duration of protection from vaccination (days)|730.5|730.5|
|Effectiveness of vaccine against hospitalization|90%|70%|
|Cumulative coverage of vaccine|user defined|user defined|

Clinical trials

(Walsh et al, 2023) https://www.nejm.org/doi/full/10.1056/NEJMoa2213836

(Papi et al, 2023) https://www.nejm.org/doi/full/10.1056/NEJMoa2209604


### Notes about interventions
The model assumes that interventions are providing protection against severe disease (hospitalization) but not against infection. The model also assumes that the interventions do not impact an individual's infectiousness if they become infected. 


# Step 4 - Visualization 
Visualize your results in the Shiny App! Code is provided in the ```RSV-Scenarios-Shiny-App``` folder and a working example is available [here](https://chelsea-doing-epi.shinyapps.io/rsv-app/).
After you have run the code in the ```3. Interventions``` folder and saved the results, launch the Shiny App by opening the ```app.R``` script in the ```RSV-Scenarios-Shiny-App``` folder and selecting the green triangle "Run App" in the top right corner of the R Script window. 

# References
Pitzer VE, Viboud C, Alonso WJ, et al. "Environmental Drivers of the Spatiotemporal Dynamics of Respiratory Syncytial Virus in the United States". PLOS Pathogens. 2015. https://doi.org/10.1371/journal.ppat.1004591
Zheng Z, Weinberger DM, Pitzer VE. "Predicted effectiveness of vaccines and extended half-life monoclonal antibodies against RSV hospitalizations in children". NJP Vaccines. 2022. https://doi.org/10.1038/s41541-022-00550-5

# Notes 

For questions or assistance, please contact chelsea.hansen@nih.gov
