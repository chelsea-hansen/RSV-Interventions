# Note
This model buildings on earlier work, please see:

Pitzer VE, Viboud C, Alonso WJ, et al. "Environmental Drivers of the Spatiotemporal Dynamics of Respiratory Syncytial Virus in the United States". PLOS Pathogens. 2015. https://doi.org/10.1371/journal.ppat.1004591

Zheng Z, Weinberger DM, Pitzer VE. "Predicted effectiveness of vaccines and extended half-life monoclonal antibodies against RSV hospitalizations in children". NJP Vaccines. 2022. https://doi.org/10.1038/s41541-022-00550-5
# Model Structure 
![RSV model diagram](https://github.com/chelsea-hansen/RSV-Interventions/assets/81387982/ebacd6f0-da6d-456f-8b46-002d0ebfd627)

# Step 1 - Data Requirements 
The data needed to run the model can be found in the Data folder. Example datasets have been provided where possible. Where data is not publicly available "dummy" datasets have been provided to demonstrate the structure of the data. Simulated data may be provided in the future
### yinit
This dataset is what you will use to initiate the burn-in for the model. The 13 age groups (<2m, 2-3m, 4-5m, 6-7, 8-9, 10-11m, 1y, 2-4y, 5-9y, 10-19y, 20-39y, 40-64y, 65+y) are divided into the model compartments based on the age distribution and size of the population when you are starting your burn-in period. This assumes the age distribution in the population is relatively stable over time. You seed 1 infection in each age group >6 months. Note - some of the middle age groups are collapsed later in the code (1-4y, 5-64y). 
### birth 
This dataset is the birth rate per 1000 population for each week you will be running the model (from initiation of burn-in to projections). You can use the annual birth rate for each week or set the first week of the year to the annual birth rate and interpolate missing weeks between years. Later in the code this will be divded by 52.17 to convert from annual to weekly scale. Note, there is a column for each age group but only the first column has data. The other columns are set to zero. 
### contact 
This model uses the contact martix from the POLYMOD study (Mossong et al. 2008. "Social Contacts and Mixing Patterns Relevant to the Spread of Infectious Diseases." PLOS Medicine. https://doi.org/10.1371/journal.pmed.0050074.) The contact matrix has been aggregated to match the age groups used in this model. 
### rsv_ts_dummy
This is a weekly count of RSV-coded hospitalizations (and optionally, ED visits) for all age groups combined. These data are not publicly available so a "dummy" dataset has been provided as a formatting example. (Simulated data may be provided in the future). It is recommended to have at least 3 years of data prior to the COVID-19 pandemic, more data is better. 
### age_distribution_dummy
This is the proporiton of RSV-coded hospitalizations (and optionally, ED visits) in each age group (with middle age groups collapsed). The proportions are separated into the pre-pandemic time period and the post-pandemic (rebound) time period. 
### scale_ed_to_hosp_dummy (optional)
The model is calibrated to RSV-coded hospitalizations. This dataset is used to scale hospitalizations to ED visits using the ratio of hospitalizations to ED visits in each age group during the pre- and post- pandemic time periods. This dataset is only necessary if you would like to produce estimates for ED visits in addition to hospitalizations. 

# Step 2 - Parameter Estimation 
This step uses the R scripts in the "Calibration" folder. rsv_dynamics.R is the function which runs the ODEs. model_calibration.R is the script which uploads the data, calls the rsv_dynamics function, and uses maximum likelihood estimation to estimate the model parameters. You will use the model_calibration.R script. Many parameters are set based on the literature. 
The model fits 8 parameters: 
 - The baseline transmission rate (beta)
 - The amplitude of seasonal forcing
 - The timing/phase of seasonal forcing
 - The duration of transplacental immunity
 - The relative risk of hospitalization among infants protected by transplacental immunity
 - The number of weekly external introductions of infection (seeding)
 - The proportion of infections leading to a reported RSV hospitalization in seniors
   - The proportion of infections leading to hospitalizations in the 1-4yrs age group is set to the same value (van Boven et al. "Estimating Transmission Parameters for Respiratory Syncytial Virus and Predicting the Impact of Maternal and Pediatric Vaccination." JID. 2020.)
   - The proportion of infections leading to hospitalizations in the 5-64 yrs age group is .06 times the proportion in adults 65+ yrs, based on influenza (Biggerstaff et al. "Estimating the Potential Effects of a Vaccine Program Against an Emerging Influenza Pandemic—United States." CID. 2015.)
 - The proportion of infections leading to a reported RSV hospitalization in infants <2 months
   - The proportion of infections in older infants is set relative to the <2 months based on the paper by Zheng et al referenced at the top of this document.
  
Save the estimated parameters for use in modeling the impact of interventions

# Step 3 - Model the Impact of Interventions 
This step uses the R scripts in the "Interventions" folder. intervention_models.R is the function that runs the ODEs (with compartments for vaccination in seniors and monoclonal antibodies in infants). intervention_scenarios.R uses the parameters from Step 2 and estimates the impact of interventions under different scenarios for intervention coverage, effectiveness, and timing. Sample scenarios are provided but these can be modified to answer your research questions. 

### Notes on interventions 
The model assumes that interventions are providing protection against severe disease (hospitalization) but not against infection. The model also assumes that the interventions do not impact an individual's infectiousness if they become infected. 

### Notes on fitting to the pandemic period 
The interventions model uses the model parameters fit during the pre-pandemic period and then projects forward. To account for distruptions during the COVID-19 pandemic we reduced the baseline transmission rate (beta) by 25% from April - June 2020 (representing the stay-at-home order) and then allowed beta to increase linearly to pre-pandemic levels between July 2020 and March 2021. We included an additional 25% beta reduction from mid-December 2021 through February 2022 to account for the emergence of the Omicron BA.1 variant. We did this because our approach fitting NPI periods suggested it was necessary and because it is likely that people adopted protective measures during this time, reducing RSV transmission. Furthermore, it is also possible that viral interactions may have resulted in a decline in RSV transmission during this period. Together with these reductions in beta, we reduced imported infections to zero from April 2020 - January 2021 and allowed for a linear increase to normal levels between February 2021 - May 2021. This approach is based on the approach taken by Zheng et al. "Estimation of the Timing and Intensity of Reemergence of Respiratory Syncytial Virus Following the COVID-19 Pandemic in the US." JAMA Netw Open. 2021. 
