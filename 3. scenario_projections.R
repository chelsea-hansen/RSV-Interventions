rm(list=ls())
library(deSolve)
library(tidyverse)
library(imputeTS)
library(cowplot)


source("MSIRS_immunization_dynamics.R")
source("projection_function.R")

#time_series = readRDS("DATA/rsv time series.rds") %>% mutate(rsv_rate = rsv_rate/100000*2267000)
dates1 = seq(from=as.Date('2017-07-01'),to=as.Date('2024-05-05'),by='week')
dates2 = seq(from=as.Date('2017-07-01'),to=as.Date('2025-05-05'),by='week')


data = readRDS("DATA/fixed_parameters.rds")
parmset = data[[1]]
yinit=data[[2]]
yinit.vector=data[[3]]

#add shared parameters 
fitted_params = readRDS("DATA/fitted_parameters_100.rds")

fit_times = seq(1, length(dates2)+104,by=1)

#vectors of immunization coverage (in doses administered)
#must be the same length as fit_times 
immu = readRDS("DATA/immunization_scenarios_2024-25.rds") 

# Scenario E (counterfactual) ----------------------------------------------

scenarioE = projection_function(RRHn = 0.2,
                                RRHmn = 0.2,
                                RRHv = 0.45,
                                RRHs = 0.2,
                                RRIn = 1,
                                RRIv = 1,
                                RRIs =1,
                                waningN = 150,
                                waningV = 180,
                                waningS = 730.5,
                                monoclonal_01 = immu$counter_monoclonal_01,
                                monoclonal_23 = immu$counter_monoclonal_23,
                                monoclonal_45 = immu$counter_monoclonal_45,
                                monoclonal_67 = immu$counter_monoclonal_67,
                                maternal_vax = immu$counter_maternal_vax, 
                                senior_vax = immu$counter_senior_vax,
                                parmset = parmset,
                                lhs_parms = fitted_params,
                                min_date = '2017-07-01',
                                max_date = '2025-05-05',
                                fit_times = fit_times) %>% 
  mutate(scenario = "scenarioE")



# ScenarioD ---------------------------------------------------------------

scenarioD = projection_function(RRHn = 0.2,
                                       RRHmn = 0.2,
                                       RRHv = 0.45,
                                       RRHs = 0.2,
                                       RRIn = 1,
                                       RRIv = 1,
                                       RRIs =1,
                                       waningN = 150,
                                       waningV = 180,
                                       waningS = 730.5,
                                       monoclonal_01 = immu$pes_monoclonal_01,
                                       monoclonal_23 = immu$pes_monoclonal_23,
                                       monoclonal_45 = immu$pes_monoclonal_45,
                                       monoclonal_67 = immu$pes_monoclonal_67,
                                       maternal_vax = immu$pes_maternal_vax, 
                                       senior_vax = immu$pes_senior_vax,
                                       parmset = parmset,
                                       lhs_parms = fitted_params,
                                       min_date = '2017-07-01',
                                       max_date = '2025-05-05',
                                       fit_times = fit_times) %>% 
  mutate(scenario = "scenarioD")



# Scenario C --------------------------------------------------------------

scenarioC = projection_function(RRHn = 0.2,
                                RRHmn = 0.2,
                                RRHv = 0.45,
                                RRHs = 0.2,
                                RRIn = 1,
                                RRIv = 1,
                                RRIs =1,
                                waningN = 150,
                                waningV = 180,
                                waningS = 730.5,
                                monoclonal_01 = immu$pes_monoclonal_01,
                                monoclonal_23 = immu$pes_monoclonal_23,
                                monoclonal_45 = immu$pes_monoclonal_45,
                                monoclonal_67 = immu$pes_monoclonal_67,
                                maternal_vax = immu$pes_maternal_vax, 
                                senior_vax = immu$opt_senior_vax,
                                parmset = parmset,
                                lhs_parms = fitted_params,
                                min_date = '2017-07-01',
                                max_date = '2025-05-05',
                                fit_times = fit_times) %>% 
  mutate(scenario = "scenarioC")





# Scenario B --------------------------------------------------------------

scenarioB = projection_function(RRHn = 0.2,
                                RRHmn = 0.2,
                                RRHv = 0.45,
                                RRHs = 0.2,
                                RRIn = 1,
                                RRIv = 1,
                                RRIs =1,
                                waningN = 150,
                                waningV = 180,
                                waningS = 730.5,
                                monoclonal_01 = immu$opt_monoclonal_01,
                                monoclonal_23 = immu$opt_monoclonal_23,
                                monoclonal_45 = immu$opt_monoclonal_45,
                                monoclonal_67 = immu$opt_monoclonal_67,
                                maternal_vax = immu$opt_maternal_vax, 
                                senior_vax = immu$pes_senior_vax,
                                parmset = parmset,
                                lhs_parms = fitted_params,
                                min_date = '2017-07-01',
                                max_date = '2025-05-05',
                                fit_times = fit_times) %>% 
  mutate(scenario = "scenarioB")


# ScenarioA  --------------------------------------------------------------
scenarioA = projection_function(RRHn = 0.2,
                                RRHmn = 0.2,
                                RRHv = 0.45,
                                RRHs = 0.2,
                                RRIn = 1,
                                RRIv = 1,
                                RRIs =1,
                                waningN = 150,
                                waningV = 180,
                                waningS = 730.5,
                                monoclonal_01 = immu$opt_monoclonal_01,
                                monoclonal_23 = immu$opt_monoclonal_23,
                                monoclonal_45 = immu$opt_monoclonal_45,
                                monoclonal_67 = immu$opt_monoclonal_67,
                                maternal_vax = immu$opt_maternal_vax, 
                                senior_vax = immu$opt_senior_vax,
                                parmset = parmset,
                                lhs_parms = fitted_params,
                                min_date = '2017-07-01',
                                max_date = '2025-05-05',
                                fit_times = fit_times) %>% 
  mutate(scenario = "scenarioA")




# Save All Data  ----------------------------------------------------------
results = rbind(scenarioE,scenarioD,scenarioC,scenarioB,scenarioA) %>% 
  mutate(date = as.Date(date))
saveRDS(results,"Results/results_2024-25.rds")

# Convert everything to rates  --------------------------------------------
pop = rowSums(parmset$yinit.matrix)
popsize = c(sum(pop[1:3]),
            sum(pop[4:6]),
            sum(pop[7:8]),
            sum(pop[9:12]),
            sum(pop[13]),
            sum(pop))
Ages = c("<6m","6-11m","1-4yrs","5-59yrs","60+yrs","All")
pops = data.frame(pop=popsize, Age = Ages)

rates = results %>% left_join(pops, by="Age") %>% 
  mutate(rate = value/pop*100000) %>% 
  group_by(Age, scenario, date) %>% 
  summarize(median = median(rate),
            lower_95 = quantile(rate,probs=0.025),
            upper_95 = quantile(rate,probs=0.975),
            lower_50 = quantile(rate,probs=0.25),
            upper_50 = quantile(rate,probs=0.75)) %>% 
  filter(date>='2024-10-01')
saveRDS(rates,"RSV-Scenarios-Shiny-App/results_rates_2024-25.rds")


counter = results %>% 
  filter(scenario=="scenarioE",date>='2024-10-01') %>%
  group_by(Age,sample) %>% 
  summarize(counterfactual = sum(value))
  
diffs = results %>% 
  filter(date>='2024-10-01') %>% 
  group_by(Age,scenario,sample) %>%
  summarize(total = sum(value)) %>% 
  left_join(counter, by=c("Age","sample")) %>%
  mutate(diff = (counterfactual - total)/counterfactual*100) %>% 
  group_by(Age,scenario) %>% 
  summarize(median = median(diff),
         lower_95 = quantile(diff,probs=0.025),
            upper_95 = quantile(diff,probs=0.975),
            lower_50 = quantile(diff,probs=0.25),
            upper_50 = quantile(diff,probs=0.75))
saveRDS(diffs,"RSV-Scenarios-Shiny-App/percentage_reduction_2024-25.rds")

# Save data for plotting cumulative coverage ------------------------------
cov_dates = seq(from=as.Date('2024-07-01'),to=as.Date('2025-05-01'),by="weeks")
coverage = immu %>% 
  mutate(pes_maternal_vax = cumsum(pes_maternal_vax)-3333,
         opt_maternal_vax = cumsum(opt_maternal_vax)-3333,
         pes_senior_vax = cumsum(pes_senior_vax),
         opt_senior_vax = cumsum(opt_senior_vax),
         pes_monoclonal = cumsum(pes_monoclonal_01+pes_monoclonal_23+pes_monoclonal_45+pes_monoclonal_67)-5336,
         opt_monoclonal = cumsum(opt_monoclonal_01+opt_monoclonal_23+opt_monoclonal_45+opt_monoclonal_67)-5336) %>% 
  select(pes_maternal_vax,opt_maternal_vax, pes_senior_vax,opt_senior_vax, pes_monoclonal,opt_monoclonal)  
coverage = coverage[471:514,]
coverage$date =cov_dates

coverage_long = coverage %>% 
  pivot_longer(cols=pes_maternal_vax:opt_monoclonal,names_to="immunization",values_to="doses") %>% 
  mutate(population = case_when(grepl("maternal",immunization)~"Maternal Vaccine",
                                grepl("senior",immunization)~"Senior Vaccine",
                                grepl("monoclonal",immunization)~"Monoclonal Abs."),
         scenario = case_when(grepl("pes",immunization)~"Pessimistic",
                              grepl("opt",immunization)~"Optimistic"))
saveRDS(coverage_long,"RSV-Scenarios-Shiny-App/coverage for figures.rds")
