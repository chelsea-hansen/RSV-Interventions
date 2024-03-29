
rm(list=ls())
library(deSolve)
library(tidyverse)
library(ggplot2)
library(imputeTS)
library(zoo)


"%notin%" = Negate('%in%')

#This part is similar to the first part of the model calibration code
## reading in and preparing data 


# Read in data and model set-up -------------------------------------------
#Compartment initiations - this time use the interventions version 
yinit = readRDS("1. Data/Demographic Data/yinit_interventions.rds") 
yinit=as.matrix(yinit)
yinit.vector = readRDS("1. Data/Demographic Data/yinit.vector_interventions.rds")

#birth rates 
birth <-readRDS('1. Data/Demographic Data/births_kingcounty.rds') %>% select(-date)
birth = as.matrix(birth)

#vector of data from burn-in through projection period 
dates = data.frame(date=seq(from=as.Date("1995-01-07"), to=as.Date("2024-06-01"), by="weeks"))
tmax0 = 1150 #starting point for when we have data (2017-01-14)
tmax1 = 1317 #end of pre-pandemic period (2020-03-28) 
tmax2 = 1496 #start of projection - first week of September 
tmax3 = 1535#projection for the 2023-2024 season (2024-06-01)
#time steps for ODE 
times <- seq(1, tmax3, by =1) 


#Upload the POLYMOD contact matrix, already aggregated to relevant age groups 
contact <- readRDS('1. Data/Demographic Data/contact_POLYMOD.rds')

#upload the other parameters - weekly seeding and death/migration rates 
seed = readRDS("1. Data/Demographic Data/other_parms.rds")[1]
um = readRDS("1. Data/Demographic Data/other_parms.rds")[2]*-1
pop = readRDS("1. Data/Demographic Data/other_parms.rds")[3]


# Define coverage scenarios ------------------------------------------
senior_coverO = .3
senior_coverP = 0.15

monoclonal_coverO = 0.4
monoclonal_coverP = 0.15

maternal_coverO = 0.25
maternal_coverP = 0.1


# Upload coverage curves and apply scenarios  -----------------------------

curves = readRDS("3. Interventions/coverage_curves_2023_24.rds")

senior_curveO = curves$senior_weekly*senior_coverO
senior_curveP = curves$senior_weekly*senior_coverP
senior_counter = rep(0,tmax3)

monoclonal_curveO = curves$monoclonal_weekly*monoclonal_coverO
monoclonal_curveP = curves$monoclonal_weekly*monoclonal_coverP
monoclonal_counter = rep(0,tmax3)

maternal_curveO = curves$maternal_weekly*maternal_coverO
maternal_curveP = curves$maternal_weekly*maternal_coverP
maternal_counter = rep(0,tmax3)

# model parameters and Latin hypercube sampling ----------------------------------

#upload parameters fit with MLE 
fitLL = readRDS("2. Calibration/parameters_6Mar24.rds")
baseline.txn.rate=6+(3*(exp(fitLL$par[1]))) / (1+exp(fitLL$par[1]))
b1=exp(fitLL$par[2])
phi=(2*pi*(exp(fitLL$par[3]))) / (1+exp(fitLL$par[3]))
reporting_rate = 1/(1+exp(-fitLL$par[4]))


pandLL = readRDS("2. Calibration/NPI_6Mar24.rds")

time1 = round(exp(pandLL$par[1]))
time2 = round(exp(pandLL$par[2]))
time3 = round(exp(pandLL$par[3]))
time4 = round(exp(pandLL$par[4]))
time5 = round(exp(pandLL$par[5]))

npi1 = 1/(1+exp(-pandLL$par[6]))
npi2 = 1/(1+exp(-pandLL$par[7]))
npi3 = 1/(1+exp(-pandLL$par[8]))
npi4 = 1/(1+exp(-pandLL$par[9]))

reporting_rate2 = reporting_rate+(reporting_rate*(exp(pandLL$par[10]))) / (1+exp(pandLL$par[10]))

introductions = data.frame(intros=c(rep(seed,tmax1),rep(0,44),rep(NA,24),rep(seed,179))) %>% 
  mutate(intros = na_interpolation(intros, method="linear"))
introductions = introductions$intros

npi = data.frame(npis=c(rep(1,tmax1),rep(npi1,13),rep(npi2,time1),rep(NA,time2),rep(npi3,time3),rep(npi4,time4),rep(NA,time5),rep(1,160)))%>% 
  mutate(npis= na_interpolation(npis, method="linear"))
npi = npi$npis

reporting = readRDS("2. Calibration/age_specific_reporting_rates.rds")
age_reporting=c(reporting[3,])#mean value to use for point estimates 

#lhs datasets
new_parms = readRDS("lhs_resampling100.rds")
rep_num=100 #which version you are using (100 replicates or 1000 replicates).
#version with 1000 replicates is recommended but will take ~2 hours to run for each scenario. Can be done late after everything is running smoothly 


# function to fit intervention scenarios (no projection intervals) ------------------------------------------
## Read in transmission dynamic model

source("3. Interventions/MSIRS_intervention_model.R")
source("3. Interventions/MSIRS_scenario_functions.R")

# Counterfactual  ---------------------------------------------------------

counterfact = interventions(birth_dose=monoclonal_counter, #vector of infants receiving birth doses of nirsevimab
                             cover_n=monoclonal_counter, #vector of catch up nirsevimab doses 
                             waningN=180, #duration of nirsevimab protection (days)
                             RRIn=1, #relative risk of infection while protected by nirsevimab
                             RRHn=0.2, #relative risk of hospitalization when protected by nirsevimab (scenario based, value doesn't matter for counterfactual because there is no coverage)
                             maternal_dose=maternal_counter, #vector of infants protected by maternal vacciantion 
                             waningV=180, #duration of protection from maternal vaccine (days_)
                             RRIv=1, #relative risk of infection when protected by maternal vaccine
                             RRHv=0.3,#realtive risk of hospitalization given infection when protected by maternal vaccine 
                             cover_s=senior_counter, #vector of seniors receiving the vaccine 
                             waningS=730.5, #duration of vaccine protection in seniors (days)
                             RRIs=1, #relative risk of infection when protected by vaccine 
                             RRHs=.1) %>%  #relative risk of hospitalization given infection when protected by vaccine 
 mutate(scenario='Counterfactual')
write.csv(counterfact,"Interventions/Counterfactual.csv")

#version with confidence intervals - will take more time to run 
counterfact_PI = interventions_PI(birth_dose=monoclonal_counter, 
                                  cover_n=monoclonal_counter, 
                                  waningN=180, 
                                  RRIn=1, 
                                  RRHn=0.2,
                                  maternal_dose=maternal_counter, 
                                  waningV=180, 
                                  RRIv=1, 
                                  RRHv=0.3,
                                  cover_s=senior_counter, 
                                  waningS=730.5, 
                                  RRIs=1, 
                                  RRHs=.1) %>% 
  mutate(scenario="Counterfactual")
write.csv(counterfact,"Interventions/Counterfactual_withProjectionsIntervals.csv")


# Scenario A  -------------------------------------------------------------
#monoclonals = optimistic, maternal vaccination = optimistic, senior vaccination = optimistic 
ScenarioA = interventions(birth_dose=monoclonal_curveO, 
                            cover_n=monoclonal_curveO, 
                            waningN=180, 
                            RRIn=1, 
                            RRHn=0.2,
                            maternal_dose=maternal_curveO, 
                            waningV=180, 
                            RRIv=1, 
                            RRHv=0.3,
                            cover_s=senior_curveO, 
                            waningS=730.5, 
                            RRIs=1, 
                            RRHs=.1) %>% 
  mutate(scenario="Scenario A")
write.csv(ScenarioA,"Interventions/ScenarioA.csv")


#version with confidence intervals - will take more time to run 
ScenarioA_PI = interventions_PI(birth_dose=monoclonal_curveO, 
                                cover_n=monoclonal_curveO, 
                                waningN=180, 
                                RRIn=1, 
                                RRHn=0.2,
                                maternal_dose=maternal_curveO, 
                                waningV=180, 
                                RRIv=1, 
                                RRHv=0.3,
                                cover_s=senior_curveO, 
                                waningS=730.5, 
                                RRIs=1, 
                                RRHs=.1)%>% 
  mutate(scenario="Scenario A")
write.csv(ScenarioA_PI,"Interventions/ScenarioA_withProjectionsIntervals.csv")

# Scenario B  -------------------------------------------------------------
#monoclonals = optimistic, maternal vaccination = pessimistic, senior vaccination = optimistic 
ScenarioB = interventions(birth_dose=monoclonal_curveO, 
                          cover_n=monoclonal_curveO, 
                          waningN=180, 
                          RRIn=1, 
                          RRHn=0.2,
                          maternal_dose=maternal_curveP, 
                          waningV=180, 
                          RRIv=1, 
                          RRHv=0.5,
                          cover_s=senior_curveO, 
                          waningS=730.5, 
                          RRIs=1, 
                          RRHs=.1)%>% 
  mutate(scenario="Scenario B")
write.csv(ScenarioB,"Interventions/ScenarioB.csv")

#version with confidence intervals - will take more time to run 
ScenarioB_PI = interventions_PI(birth_dose=monoclonal_curveO, 
                                cover_n=monoclonal_curveO, 
                                waningN=180, 
                                RRIn=1, 
                                RRHn=0.2,
                                maternal_dose=maternal_curveP, 
                                waningV=180, 
                                RRIv=1, 
                                RRHv=0.5,
                                cover_s=senior_curveO, 
                                waningS=730.5, 
                                RRIs=1, 
                                RRHs=.1)%>% 
  mutate(scenario="Scenario B")
write.csv(ScenarioB_PI,"Interventions/ScenarioB_withProjectionsIntervals.csv")


# Scenario C --------------------------------------------------------------
#monoclonals = pessimistic, maternal vaccination = optimistic, senior vaccination = optimistic 
ScenarioC = interventions(birth_dose=monoclonal_curveP, 
                          cover_n=monoclonal_curveP, 
                          waningN=180, 
                          RRIn=1, 
                          RRHn=0.4,
                          maternal_dose=maternal_curveO, 
                          waningV=180, 
                          RRIv=1, 
                          RRHv=0.3,
                          cover_s=senior_curveO, 
                          waningS=730.5, 
                          RRIs=1, 
                          RRHs=.1)%>% 
  mutate(scenario="Scenario C")
write.csv(ScenarioC,"Interventions/ScenarioC.csv")

#version with confidence intervals - will take more time to run 
ScenarioC_PI = interventions_PI(birth_dose=monoclonal_curveP, 
                                cover_n=monoclonal_curveP, 
                                waningN=180, 
                                RRIn=1, 
                                RRHn=0.4,
                                maternal_dose=maternal_curveO, 
                                waningV=180, 
                                RRIv=1, 
                                RRHv=0.3,
                                cover_s=senior_curveO, 
                                waningS=730.5, 
                                RRIs=1, 
                                RRHs=.1)%>% 
  mutate(scenario="Scenario C")
write.csv(ScenarioC_PI,"Interventions/ScenarioC_withProjectionsIntervals.csv")


# Scenario D --------------------------------------------------------------
#monoclonals = pessimistic, maternal vaccination = pessimistic, senior vaccination = optimistic 
ScenarioD = interventions(birth_dose=monoclonal_curveP, 
                          cover_n=monoclonal_curveP, 
                          waningN=180, 
                          RRIn=1, 
                          RRHn=0.4,
                          maternal_dose=maternal_curveP, 
                          waningV=180, 
                          RRIv=1, 
                          RRHv=0.5,
                          cover_s=senior_curveO, 
                          waningS=730.5, 
                          RRIs=1, 
                          RRHs=.1)%>% 
  mutate(scenario="Scenario D")
write.csv(ScenarioD,"Interventions/ScenarioD.csv")

#version with confidence intervals - will take more time to run 
ScenarioD_PI = interventions_PI(birth_dose=monoclonal_curveP, 
                                cover_n=monoclonal_curveP, 
                                waningN=180, 
                                RRIn=1, 
                                RRHn=0.4,
                                maternal_dose=maternal_curveP, 
                                waningV=180, 
                                RRIv=1, 
                                RRHv=0.5,
                                cover_s=senior_curveO, 
                                waningS=730.5, 
                                RRIs=1, 
                                RRHs=.1)%>% 
  mutate(scenario="Scenario D")
write.csv(ScenarioD_PI,"Interventions/ScenarioD_withProjectionsIntervals.csv")



# Scenario E --------------------------------------------------------------
#monoclonals = optimistic, maternal vaccination = optimistic, senior vaccination = pessimistic
ScenarioE = interventions(birth_dose=monoclonal_curveO, 
                          cover_n=monoclonal_curveO, 
                          waningN=180, 
                          RRIn=1, 
                          RRHn=0.2,
                          maternal_dose=maternal_curveO, 
                          waningV=180, 
                          RRIv=1, 
                          RRHv=0.3,
                          cover_s=senior_curveP, 
                          waningS=730.5, 
                          RRIs=1, 
                          RRHs=.3)%>% 
  mutate(scenario="Scenario E")
write.csv(ScenarioE,"Interventions/ScenarioE.csv")

#version with confidence intervals - will take more time to run 
ScenarioE_PI = interventions_PI(birth_dose=monoclonal_curveO, 
                                cover_n=monoclonal_curveO, 
                                waningN=180, 
                                RRIn=1, 
                                RRHn=0.2,
                                maternal_dose=maternal_curveO, 
                                waningV=180, 
                                RRIv=1, 
                                RRHv=0.3,
                                cover_s=senior_curveP, 
                                waningS=730.5, 
                                RRIs=1, 
                                RRHs=.3)%>% 
  mutate(scenario="Scenario E")
write.csv(ScenarioE_PI,"Interventions/ScenarioE_withProjectionsIntervals.csv")

# Scenario F --------------------------------------------------------------
#monoclonals = optimistic, maternal vaccination = pessimistic, senior vaccination = pessimistic

ScenarioF = interventions(birth_dose=monoclonal_curveO, 
                          cover_n=monoclonal_curveO, 
                          waningN=180, 
                          RRIn=1, 
                          RRHn=0.2,
                          maternal_dose=maternal_curveP, 
                          waningV=180, 
                          RRIv=1, 
                          RRHv=0.5,
                          cover_s=senior_curveP, 
                          waningS=730.5, 
                          RRIs=1, 
                          RRHs=.3)%>% 
  mutate(scenario="Scenario F")
write.csv(ScenarioF,"Interventions/ScenarioF.csv")

#version with confidence intervals - will take more time to run 
ScenarioF_PI = interventions_PI(birth_dose=monoclonal_curveO, 
                                cover_n=monoclonal_curveO, 
                                waningN=180, 
                                RRIn=1, 
                                RRHn=0.2,
                                maternal_dose=maternal_curveP, 
                                waningV=180, 
                                RRIv=1, 
                                RRHv=0.5,
                                cover_s=senior_curveP, 
                                waningS=730.5, 
                                RRIs=1, 
                                RRHs=.3)%>% 
  mutate(scenario="Scenario F")
write.csv(ScenarioF_PI,"Interventions/ScenarioF_withProjectionsIntervals.csv")


# Scenario G --------------------------------------------------------------
#monoclonals = pessimistic, maternal vaccination = optimistic, senior vaccination = pessimistic
ScenarioG = interventions(birth_dose=monoclonal_curveO, 
                          cover_n=monoclonal_curveO, 
                          waningN=180, 
                          RRIn=1, 
                          RRHn=0.4,
                          maternal_dose=maternal_curveP, 
                          waningV=180, 
                          RRIv=1, 
                          RRHv=0.3,
                          cover_s=senior_curveP, 
                          waningS=730.5, 
                          RRIs=1, 
                          RRHs=.3)%>% 
  mutate(scenario="Scenario G")
write.csv(ScenarioG,"Interventions/ScenarioG.csv")

#version with confidence intervals - will take more time to run 
ScenarioG_PI = interventions_PI(birth_dose=monoclonal_curveP, 
                                cover_n=monoclonal_curveP, 
                                waningN=180, 
                                RRIn=1, 
                                RRHn=0.4,
                                maternal_dose=maternal_curveO, 
                                waningV=180, 
                                RRIv=1, 
                                RRHv=0.3,
                                cover_s=senior_curveP, 
                                waningS=730.5, 
                                RRIs=1, 
                                RRHs=.3)%>% 
  mutate(scenario="Scenario G")
write.csv(ScenarioG_PI,"Interventions/ScenarioG_withProjectionsIntervals.csv")




# Scenario H --------------------------------------------------------------
#monoclonals = pessimistic, maternal vaccination = pessimistic, senior vaccination = pessimistic
ScenarioH = interventions(birth_dose=monoclonal_curveP, 
                          cover_n=monoclonal_curveP, 
                          waningN=180, 
                          RRIn=1, 
                          RRHn=0.4,
                          maternal_dose=maternal_curveP, 
                          waningV=180, 
                          RRIv=1, 
                          RRHv=0.5,
                          cover_s=senior_curveP, 
                          waningS=730.5, 
                          RRIs=1, 
                          RRHs=.3)%>% 
  mutate(scenario="Scenario H")
write.csv(ScenarioH,"Interventions/ScenarioH.csv")

#version with confidence intervals - will take more time to run 
ScenarioH_PI = interventions_PI(birth_dose=monoclonal_curveP, 
                                cover_n=monoclonal_curveP, 
                                waningN=180, 
                                RRIn=1, 
                                RRHn=0.4,
                                maternal_dose=maternal_curveP, 
                                waningV=180, 
                                RRIv=1, 
                                RRHv=0.5,
                                cover_s=senior_curveP, 
                                waningS=730.5, 
                                RRIs=1, 
                                RRHs=.3)%>% 
  mutate(scenario="Scenario H")
write.csv(ScenarioH_PI,"Interventions/ScenarioH_withProjectionsIntervals.csv")


# Combine all scenarios ---------------------------------------------------
#combine all point estimates 
all_scenarios=rbind(counterfact,ScenarioA,ScenarioB,ScenarioC,ScenarioD,
                    ScenarioE,ScenarioF,ScenarioG,ScenarioH)
saveRDS(all_scenarios,"scenario_point_estimates.rds")
#combine all confidence intervals 
all_scenarios_PI=rbind(counterfact_PI,ScenarioA_PI,ScenarioB_PI,ScenarioC_PI,ScenarioD_PI,
                    ScenarioE_PI,ScenarioF_PI,ScenarioG_PI,ScenarioH_PI) %>% 
  filter(!is.na(sample))
saveRDS(all_scenarios_PI,"scenarios_with_PI.rds")


