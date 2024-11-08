rm(list=ls())
library(deSolve)
library(tidyverse)
library(zoo)
library(imputeTS)
library(cowplot)

source("R/MSIRS_immunization_dynamics.R")
source("R/projection_function.R")



dates1 = seq(from=as.Date('2017-07-02'),to=as.Date('2024-04-28'),by='week')


data = readRDS("DATA/fixed_parameters.rds")
parmset = data[[1]]
yinit=data[[2]]
yinit.vector=data[[3]]


#add shared parameters 
fitted_params = readRDS("DATA/fitted_parameters_100.rds")

fit_times = seq(1, length(dates1)+104,by=1)

# observed data 
time_series = readRDS("DATA/time_series.rds") %>% 
  mutate(rsv_rate = hosp_rate/100000*2267000)


pop = rowSums(parmset$yinit.matrix)
popsize = c(sum(pop[1:3]),
            sum(pop[4:6]),
            sum(pop[7:8]),
            sum(pop[9:12]),
            sum(pop[13]),
            sum(pop[14]),
            sum(pop))
Ages = c("<6m","6-11m","1-4yrs","5-59yrs","60-74yrs","75+yrs","All")
pops = data.frame(pop=popsize, Age = Ages)

# Observed scenario 2023-24 -----------------------------------------------
#convert coverage to estimated doses
immu = readRDS("DATA/weekly_immunizations_public.rds") %>% 
  mutate(week_maternal = week_maternal/100*11900,
         week_monoclonal_birth = week_monoclonal_birth/100*11900,
         week_monoclonal_catchup = week_monoclonal_catchup/100*13800,
         week_seniors_60_74 = week_seniors_60_74/100*313754,
         week_seniors_75 = week_seniors_75/100*132708)
         

change_recs = c(rep(.33,23),rep(.25,13)) #divide between 3 from Oct - Dec because changes in recs 
change_recs2 = c(rep(.34,23),rep(.25,13))
change_recs3 = c(rep(0,23),rep(.25,13)) # in January go back to dividing by 4 

obs_senior75_vax <- c(rep(0, length(dates1) + 104-40),immu$week_seniors_75*.87,rep(0,4))#13% of doses are "wasted" because given to people who already have immunity (in R compartment)
obs_senior60_vax <- c(rep(0, length(dates1) + 104-40),immu$week_seniors_60_74*.8,rep(0,4))#20% of doses are "wasted" because given to people who already have immunity (in R compartment)
obs_maternal_vax <- c(rep(0, length(dates1) + 104-40),immu$week_maternal,rep(0,4))

obs_monoclonal_01 =  c(rep(0, length(dates1) + 104-40),immu$week_monoclonal_catchup*change_recs,rep(0,4))
obs_monoclonal_23 =  c(rep(0, length(dates1) + 104-40),immu$week_monoclonal_catchup*change_recs,rep(0,4))
obs_monoclonal_45 = c(rep(0, length(dates1) + 104-40),immu$week_monoclonal_catchup*change_recs2,rep(0,4))
obs_monoclonal_67 =  c(rep(0, length(dates1) + 104-40),immu$week_monoclonal_catchup*change_recs3,rep(0,4))
obs_monoclonal_birth = c(rep(0, length(dates1) + 104-40),immu$week_monoclonal_birth,rep(0,4))


observed_main = projection_function(RRHn1 = 0.01,
                                    RRHn2 = 0.01,
                                    RRHv1 = 0.3,
                                    RRHv2 = 0.3,
                                    RRHs = 0.25,
                                    RRIn = 1,
                                    RRIv = 1,
                                    RRIs =1,
                                    waningN1 = 90,
                                    waningN2 = 90,
                                    waningV1 = 90,
                                    waningV2 = 90,
                                    waningS = 730.5,
                                    monoclonal_birth = obs_monoclonal_birth,
                                    monoclonal_01 = obs_monoclonal_01,
                                    monoclonal_23 = obs_monoclonal_23,
                                    monoclonal_45 = obs_monoclonal_45,
                                    monoclonal_67 = obs_monoclonal_67,
                                    maternal_vax = obs_maternal_vax, 
                                    senior_vax_75 = obs_senior75_vax,
                                    senior_vax_60_74 = obs_senior60_vax,
                                    parmset = parmset,
                                    lhs_parms = fitted_params,
                                    min_date = '2017-07-01',
                                    max_date = '2024-04-28',
                                    fit_times = fit_times) %>% 
  mutate(scenario = "observed_main")

ggplot(data=observed_main %>% filter(Age=="All") %>% mutate(date=as.Date(date)))+
  geom_line(aes(x=date, y=value, group=sample))


# Counterfactual  ---------------------------------------------------------

counter_monoclonal_01 =  rep(0, length(dates1) + 104)
counter_monoclonal_23 =  rep(0, length(dates1) + 104)
counter_monoclonal_45 =  rep(0, length(dates1) + 104)
counter_monoclonal_67 =  rep(0, length(dates1) + 104)
counter_monoclonal_birth =  rep(0, length(dates1) + 104)
counter_maternal_vax <- rep(0, length(dates1) + 104)
counter_senior_vax75 <- rep(0, length(dates1) + 104)
counter_senior_vax60 <- rep(0, length(dates1) + 104)


counterfactual = projection_function(RRHn1 = 0.01,
                                     RRHn2 = 0.01,
                                     RRHv1 = 0.3,
                                     RRHv2 = 0.3,
                                     RRHs = 0.25,
                                     RRIn = 1,
                                     RRIv = 1,
                                     RRIs =1,
                                     waningN1 = 90,
                                     waningN2 = 90,
                                     waningV1 = 90,
                                     waningV2 = 90,
                                     waningS = 730.5,
                                     monoclonal_birth = counter_monoclonal_birth,
                                     monoclonal_01 = counter_monoclonal_01,
                                     monoclonal_23 = counter_monoclonal_23,
                                     monoclonal_45 = counter_monoclonal_45,
                                     monoclonal_67 = counter_monoclonal_67,
                                     maternal_vax = counter_maternal_vax, 
                                     senior_vax_75 = counter_senior_vax75,
                                     senior_vax_60 = counter_senior_vax60,
                                     parmset = parmset,
                                     lhs_parms = fitted_params,
                                     min_date = '2017-07-01',
                                     max_date = '2024-04-28',
                                     fit_times = fit_times) %>% 
  mutate(scenario = "counterfactual")

ggplot(data=counterfactual %>% filter(Age=="All") %>% mutate(date=as.Date(date)))+
  geom_line(aes(x=date, y=value, group=sample))


results = rbind(observed_main, counterfactual)
saveRDS(results,"RESULTS/results 100 23-24 check.rds")
# Monoclonal_only ---------------------------------------------------------
monoclonal_only = projection_function(RRHn1 = 0.01,
                                    RRHn2 = 0.01,
                                    RRHv1 = 0.3,
                                    RRHv2 = 0.3,
                                    RRHs = 0.25,
                                    RRIn = 1,
                                    RRIv = 1,
                                    RRIs =1,
                                    waningN1 = 90,
                                    waningN2 = 90,
                                    waningV1 = 90,
                                    waningV2 = 90,
                                    waningS = 730.5,
                                    monoclonal_birth = obs_monoclonal_birth,
                                    monoclonal_01 = obs_monoclonal_01,
                                    monoclonal_23 = obs_monoclonal_23,
                                    monoclonal_45 = obs_monoclonal_45,
                                    monoclonal_67 = obs_monoclonal_67,
                                    maternal_vax = counter_maternal_vax, 
                                    senior_vax_75 = obs_senior75_vax,
                                    senior_vax_60_74 = obs_senior60_vax,
                                    parmset = parmset,
                                    lhs_parms = fitted_params,
                                    min_date = '2017-07-01',
                                    max_date = '2024-04-28',
                                    fit_times = fit_times) %>% 
  mutate(scenario = "monoclonal_only")

# Maternal Vax only  ------------------------------------------------------
maternal_only = projection_function(RRHn1 = 0.01,
                                    RRHn2 = 0.01,
                                    RRHv1 = 0.3,
                                    RRHv2 = 0.3,
                                    RRHs = 0.25,
                                    RRIn = 1,
                                    RRIv = 1,
                                    RRIs =1,
                                    waningN1 = 90,
                                    waningN2 = 90,
                                    waningV1 = 90,
                                    waningV2 = 90,
                                    waningS = 730.5,
                                    monoclonal_birth = counter_monoclonal_birth,
                                    monoclonal_01 = counter_monoclonal_01,
                                    monoclonal_23 = counter_monoclonal_23,
                                    monoclonal_45 = counter_monoclonal_45,
                                    monoclonal_67 = counter_monoclonal_67,
                                    maternal_vax = obs_maternal_vax, 
                                    senior_vax_75 = obs_senior75_vax,
                                    senior_vax_60_74 = obs_senior60_vax,
                                    parmset = parmset,
                                    lhs_parms = fitted_params,
                                    min_date = '2017-07-01',
                                    max_date = '2024-04-28',
                                    fit_times = fit_times) %>% 
  mutate(scenario = "maternal_only")


results = rbind(observed_main,counterfactual,
                monoclonal_only, maternal_only)
table(results$scenario)
saveRDS(results,"RESULTS/results_23-24_100replicates.rds")

