
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

# Step 1 - Data Requirements 
## RSV Data

To run the model you will need to have a weekly (or monthly) time series of RSV hospitalizations (or ED visits) and an age distribution of RSV hospitalizations (or ED visits if that is the time series you are using). These should becount data not rates. We used data from King County, Washington. We hope to post these as example datasets soon. For the weekly time series it is best if you can have at least 3 years of data prior to the COVID-19 pandemic, however the code should work with a slightly shorter time series.  Note, if using a smoothed time series or concverting from hospitalization rates, you must rounded to the nearest whole number. The model fitting procedure uses a Poisson distribution and you will get an error if the RSV time series has not been rounded to whole numbers. 

The age distribution should be disaggregated into 5 age groups (<6 months, 6-11 months, 1-4 years, 5-59 years, 60+ years). The code could be modified to use different age groups.

## Demographic data
The code will also require birth rates, net migration rates, and the age-specific population distribution. The ```1.data_prep.R``` R script will pull and format all of the necessary data using the ```tidycensus``` R package. This code will create 3 datasets and save them together as a list. 1. The first dataset is all of the fixed parameter values. 
2. The second dataset you will create is ```yinit.rds```. This dataset divides the population from 13 age groups (<2m, 2-3m, 4-5m, 6-7, 8-9, 10-11m, 1y, 2-4y, 5-9y, 10-19y, 20-39y, 40-59y, 60+y) into starting values for the model compartments. Note: these age groups are not the same as the age groups from the RSV age distributions. 
<img src="https://github.com/chelsea-hansen/RSV-Interventions/assets/81387982/4f5419ff-a5d0-464b-b241-b7c7e32824e8" width="80%" height="80%" align="center">
3.The code will also save another format of this dataset ```yinit.vector.rds``` 

The annual birth rate is converted to a weekly number of births and is used to introduce new individuals into the <2m age class. The model assumes that individuals age exponentially into the next age class with the rate of aging equal to the inverse of the time spent in each age class. The duration of the oldest age class was set to 20 years. The net migration rate was applied uniformly across age classes. We used an expanded version of the contact matrix described by Mossong et al to define contacts between age classes. 

To run the next step (Calibration) you will need 3 things. 
1. The RSV time series 
2. The RSV age distribution 
3. The list you saved (containing 3 datasets) using the ```1.data_prep.R``` R script

   
# Step 2 - Calibration
Part 1: This step is the trickiest step in the process and likely where you will have the most difficulty. But the good news is you only need to run this once and then you can update and change your immunization scenarios as much as you like! 
In this step you will use the ```2. Calibration``` R script. This script uses Maximum Likelihood Estimation to fit 11 parameters (see table above). These include: 

 1.  The baseline transmission rate (&beta;).
 2.  The amplitude of seasonal forcing (*b*1).
 3.  The timing/phase of seasonal forcing (&phi;), bounded between 0 and 2&pi;.
 4-7.  Age-specific reporting rates for the proportion of infections leading to reported hospitalizations, bounded between 0 and 1.
 8-11. Reductions in contacts for 4 periods during the 2020-2022, bounded between 0 and 1.

Assumptions when fitting contact reductions: 
 1. It is assumed that there is a reduction in contacts from April - June 2020, corresponding with stay-at-home-orders. The model will fit the magnitude of the reduction.
 2. It is assumed that following the end of stay-at-home orders, contacts remained lower into 2021. The model will fit the magnitude of the reduction. 
 3. It is assumed that contacts began to return to normal in 2021. The model will estimate the level of contacts by the end of 2021.
 4. It is assumed that there was another drop in contacts coinciding with the emergence of the Omicron variant in winter 2021-2022. The model will fit the magnitude of this reduction.
 5. It is assumed that following the initial Omicron wave, contacts gradually increased through 2022, returning to pre-pandemic levels by the fall of 2022. 

After completing parts 1 and 2 of the calibration step, you will have a figure that looks like this: 

<img src="https://github.com/chelsea-hansen/RSV-Interventions/assets/81387982/08c3d85f-78ec-4203-a02a-68013e36ae75" width="80%" height="80%" align="center">

After fitting these 11 parameters the script uses Latin Hypercube Sampling to sample parameter values from a plausible range around the fitted parameters (+/- 10%) for all parameters except (&phi;) (+/- 2%), which is more sensitive to small changes. Two versions are saved, one with 100 samples and the other with 1000 samples. These sets of parameter values will be used for estimating confidence intervals around the model projections. 

Note: Although example datasets are not posted at this time, the fitted parameter values are available in the ```Data``` folder. You can proceed with Step 3 using these fitted parameters, you do not need the original time series data. 

# Step 3 - Scenario Projections 
This step uses the ```3.scenario_projections``` R script.. The script runs 5 scenarios for the 2024-25 RSV season: 1 counterfactual (no interventions) plus 4 combinations of optimistic and pessimistic coverage and for the three immunization strategies (RSV vaccination for adults >60 years, RSV vaccination for pregnant women, RSV monoclonal antibodies for infants <8 months). See details in the table below. The parameters for the immunizations have been fixed based on clinical trials (see tables). in the ```Data``` folder there is a dataset with the weekly number of administered immunization doses for each scenario. It is suggested to start with the version of the fitted parameters that runs 100 trajectories as this will be quicker. 

Scenario Overview: 
|Infant Immunization|Senior Immunization||
|-|-|-|
||Senior Vax = Optimistic|Senior Vax = Pessimistic|
|Monoclonal Abs & Maternal Vax =  Optimistic|A|B|
|Monoclonal Abs & Maternal Vax = Pessimistic|C|D|

## Infant Immunizations
The model now assumes that a proportion of infants are born to vaccinated mothers (Mv Compartment) and a proportion of infants receive monoclonal antibodies at or shortly after birth (Mn Compartment). These infants retain the same protection against infection as the M compartment, but have additional protection against hospitalization given infection. Additionally, some infants receive a monoclonal antibody a few months after birth (N Compartment). These infants do not have any protection against infection, but they do have protection against hospitalization given infection. When the protection wanes from these compartments (Mv, Mn, N), infants move to the Si compartment. This compartment is functionally the same as the S0 compartment, but ensures that infants do not receive two interventions. 

<img src="https://github.com/chelsea-hansen/RSV-Interventions/assets/81387982/8c830c59-3d7d-4100-8a7d-f4f728656f62" width="50%" height="50%" align="left">


|Parameter|Optimistic Value|Pessimistic Value|
|---------|----------------|----|
|Duration of protection from monoclonal antibodies (days)|150|150|
|Effectiveness of monoclonal antibodies against hospitalization|80%|80%|
|Cumulative coverage of monoclonal antibodies|75%|25%|
|Duration or protection from maternal vaccination (days)|180|180|
|Effectiveness of maternal vaccination against hospitalization|55%|55%|
|Cumulative coverage of maternal vaccination|50%|30%|

## Senior Vaccination 
The vaccination compartment (Vs1) draws seniors from the S3 compartment. Current data suggests that the vaccine is effective for at least 2 seasons. Seniors spend approximately 1 year in the Vs1 compartment before waning to the Vs2 compartment for another year and then returning to the S3 compartment. Because the vaccine lasts for 2 years the optimistic and pessimistic coverage for the 2024-25 season includes the observed 25% coverage during the 2023-24 season. 

<img src="https://github.com/chelsea-hansen/RSV-Interventions/assets/81387982/a569dc65-8398-4059-8782-917b27c42041" width="50%" height="50%" align="left">


|Parameter|Optimistic Value|Pessimistic Value|
|---------|----------------|-----------------|
|Duration of protection from vaccination (days)|730.5|730.5|
|Effectiveness of vaccine against hospitalization|80%|80%|
|Cumulative coverage of vaccine|40%|30%|


### Notes about interventions
The model assumes that interventions are providing protection against severe disease (hospitalization) but not against infection. The model also assumes that the interventions do not impact an individual's infectiousness if they become infected. 


# Step 4 - Visualization 
Visualize your results in the Shiny App! Code is provided in the ```RSV-Scenarios-Shiny-App``` folder and a working example is available [here](https://chelsea-doing-epi.shinyapps.io/RSV-Scenarios-Shiny-App/).
After you have run the ```3.scenario_projections``` R script, save the results in the ```RSV-Scenarios-Shiny-App``` folder (examples already included). Launch the Shiny App by opening the ```app.R``` script in the ```RSV-Scenarios-Shiny-App``` folder and selecting the green triangle "Run App" in the top right corner of the R Script window. 

# References
Pitzer VE, Viboud C, Alonso WJ, et al. "Environmental Drivers of the Spatiotemporal Dynamics of Respiratory Syncytial Virus in the United States". PLOS Pathogens. 2015. https://doi.org/10.1371/journal.ppat.1004591

Zheng Z, Weinberger DM, Pitzer VE. "Predicted effectiveness of vaccines and extended half-life monoclonal antibodies against RSV hospitalizations in children". NJP Vaccines. 2022. https://doi.org/10.1038/s41541-022-00550-5

# Notes 

For questions or assistance, please contact chelsea.hansen@nih.gov
