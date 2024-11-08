rm(list=ls())
library(deSolve)
library(tidyverse)
library(imputeTS)
library(cowplot)

'%notin%' = Negate('%in%')
source("R/MSIRS_immunization_dynamics.R")
source("R/projection_function.R")

time_series = readRDS("DATA/time_series_public.rds") %>% mutate(hosp_rate = hosp_rate/100000*2267000)
dates1 = seq(from=as.Date('2017-07-01'),to=as.Date('2024-04-28'),by='week')
dates2 = seq(from=as.Date('2017-07-01'),to=as.Date('2025-04-28'),by='week')


data = readRDS("DATA/fixed_parameters.rds")
parmset = data[[1]]
yinit=data[[2]]
yinit.vector=data[[3]]

#add shared parameters 
fitted_params = readRDS("DATA/fitted_parameters_100.rds")

fit_times = seq(1, length(dates2)+104,by=1)

immu = readRDS("DATA/weekly_immunizations_public.rds") %>% 
  mutate(week_maternal = week_maternal/100*11900,
         week_monoclonal_birth = week_monoclonal_birth/100*11900,
         week_monoclonal_catchup = week_monoclonal_catchup/100*13800,
         week_seniors_60_74 = week_seniors_60_74/100*313754,
         week_seniors_75 = week_seniors_75/100*132708,
         cum_maternal = cum_maternal/100*11900,
         cum_monoclonal_birth = cum_monoclonal_birth/100*11900,
         cum_monoclonal_catchup = cum_monoclonal_catchup/100*13800,
         cum_seniors_60_74 = cum_seniors_60_74/100*313754,
         cum_seniors_75 = cum_seniors_75/100*132708,)

#Move maternal vaccination up by 1 month and stop all immunizations at end of March

immu2 = immu  %>% mutate(cum_maternal=lead(cum_maternal,4),
                         week_maternal = lead(week_maternal,4)) %>% 
  replace(is.na(.),0) %>% 
  mutate(scale_mat = week_maternal/ max(cum_maternal),
         scale_birth = week_monoclonal_birth/max(cum_monoclonal_birth),
         scale_catchup = week_monoclonal_catchup/max(cum_monoclonal_catchup))

cum_sen_75 = max(immu2$cum_seniors_75)
cum_sen_60 = max(immu2$cum_seniors_60_74)



# Optimistic and Pessimistic Coverage  ------------------------------------

#Counterfactuals 
counter_monoclonal_01 =  rep(0, length(dates2) + 104)
counter_monoclonal_23 =  rep(0, length(dates2) + 104)
counter_monoclonal_45 =  rep(0, length(dates2) + 104)
counter_monoclonal_67 =  rep(0, length(dates2) + 104)
counter_monoclonal_birth =  rep(0, length(dates2) + 104)
counter_maternal_vax <- rep(0, length(dates2) + 104)
counter_senior_vax75 <- rep(0, length(dates2) + 104)
counter_senior_vax60 <- rep(0, length(dates2) + 104)


#Optimistic 
opt_sen_75 = c(rep(0, length(dates2) + 104-41),cum_sen_75*.87,immu2$week_seniors_75*.87,rep(0,4))
opt_sen_60_74 = c(rep(0, length(dates2) + 104-41),cum_sen_60*.8,immu2$week_seniors_60_74*.8,rep(0,4))
opt_maternal_vax = c(rep(0, length(dates2) + 104-40),immu2$scale_mat*4800,rep(0,4))
opt_monoclonal_birth = c(rep(0, length(dates2) + 104-31),rep(5550/27,27),rep(0,4))
opt_monoclonal_01 = c(rep(0, length(dates2) + 104-40),immu2$scale_catchup*10200*.25,rep(0,4))
opt_monoclonal_23 = c(rep(0, length(dates2) + 104-40),immu2$scale_catchup*10200*.25,rep(0,4))
opt_monoclonal_45 = c(rep(0, length(dates2) + 104-40),immu2$scale_catchup*10200*.25,rep(0,4))
opt_monoclonal_67 = c(rep(0, length(dates2) + 104-40),immu2$scale_catchup*10200*.25,rep(0,4))


#Pessimistic 
pes_sen_75 = c(rep(0, length(dates2) + 104-41),cum_sen_75*.87,immu2$week_seniors_75*.87*.5,rep(0,4))
pes_sen_60_74 = c(rep(0, length(dates2) + 104-41),cum_sen_60*.8,immu2$week_seniors_60_74*.8*.5,rep(0,4))
pes_maternal_vax = c(rep(0, length(dates2) + 104-40),immu2$scale_mat*3200,rep(0,4))
pes_monoclonal_birth = c(rep(0, length(dates2) + 104-31),rep(3700/27,27),rep(0,4))
pes_monoclonal_01 = c(rep(0, length(dates2) + 104-40),immu2$scale_catchup*3400*.25,rep(0,4))
pes_monoclonal_23 = c(rep(0, length(dates2) + 104-40),immu2$scale_catchup*3400*.25,rep(0,4))
pes_monoclonal_45 = c(rep(0, length(dates2) + 104-40),immu2$scale_catchup*3400*.25,rep(0,4))
pes_monoclonal_67 = c(rep(0, length(dates2) + 104-40),immu2$scale_catchup*3400*.25,rep(0,4))


# Counterfactual for 2024-25 ----------------------------------------------


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
                                     max_date = '2025-04-28',
                                     fit_times = fit_times) %>% 
  mutate(scenario = "counterfactual")

ggplot(data=counterfactual %>% filter(Age=="All") %>% mutate(date=as.Date(date)))+
  geom_line(aes(x=date, y=value, group=sample))



# ScenarioD ---------------------------------------------------------------


scenarioD = projection_function(RRHn1 = 0.01,
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
                                     monoclonal_birth = pes_monoclonal_birth,
                                     monoclonal_01 = pes_monoclonal_01,
                                     monoclonal_23 = pes_monoclonal_23,
                                     monoclonal_45 = pes_monoclonal_45,
                                     monoclonal_67 = pes_monoclonal_67,
                                     maternal_vax = pes_maternal_vax, 
                                     senior_vax_75 = pes_sen_75,
                                     senior_vax_60 = pes_sen_60_74,
                                     parmset = parmset,
                                     lhs_parms = fitted_params,
                                     min_date = '2017-07-01',
                                     max_date = '2025-04-28',
                                     fit_times = fit_times) %>% 
  mutate(scenario = "Scenario D")

ggplot(data=scenarioD %>% filter(Age=="All") %>% mutate(date=as.Date(date)))+
  geom_line(aes(x=date, y=value, group=sample))




# Scenario C --------------------------------------------------------------


scenarioC = projection_function(RRHn1 = 0.01,
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
                                monoclonal_birth = pes_monoclonal_birth,
                                monoclonal_01 = pes_monoclonal_01,
                                monoclonal_23 = pes_monoclonal_23,
                                monoclonal_45 = pes_monoclonal_45,
                                monoclonal_67 = pes_monoclonal_67,
                                maternal_vax = pes_maternal_vax, 
                                senior_vax_75 = opt_sen_75,
                                senior_vax_60 = opt_sen_60_74,
                                parmset = parmset,
                                lhs_parms = fitted_params,
                                min_date = '2017-07-01',
                                max_date = '2025-04-28',
                                fit_times = fit_times) %>% 
  mutate(scenario = "Scenario C")

ggplot(data=scenarioC %>% filter(Age=="All") %>% mutate(date=as.Date(date)))+
  geom_line(aes(x=date, y=value, group=sample))


# Scenario B --------------------------------------------------------------
scenarioB = projection_function(RRHn1 = 0.01,
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
                                monoclonal_birth = opt_monoclonal_birth,
                                monoclonal_01 = opt_monoclonal_01,
                                monoclonal_23 = opt_monoclonal_23,
                                monoclonal_45 = opt_monoclonal_45,
                                monoclonal_67 = opt_monoclonal_67,
                                maternal_vax = opt_maternal_vax, 
                                senior_vax_75 = pes_sen_75,
                                senior_vax_60 = pes_sen_60_74,
                                parmset = parmset,
                                lhs_parms = fitted_params,
                                min_date = '2017-07-01',
                                max_date = '2025-04-28',
                                fit_times = fit_times) %>% 
  mutate(scenario = "Scenario B")

ggplot(data=scenarioB %>% filter(Age=="All") %>% mutate(date=as.Date(date)))+
  geom_line(aes(x=date, y=value, group=sample))



# ScenarioA  --------------------------------------------------------------
scenarioA = projection_function(RRHn1 = 0.01,
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
                                monoclonal_birth = opt_monoclonal_birth,
                                monoclonal_01 = opt_monoclonal_01,
                                monoclonal_23 = opt_monoclonal_23,
                                monoclonal_45 = opt_monoclonal_45,
                                monoclonal_67 = opt_monoclonal_67,
                                maternal_vax = opt_maternal_vax, 
                                senior_vax_75 = opt_sen_75,
                                senior_vax_60 = opt_sen_60_74,
                                parmset = parmset,
                                lhs_parms = fitted_params,
                                min_date = '2017-07-01',
                                max_date = '2025-04-28',
                                fit_times = fit_times) %>% 
  mutate(scenario = "Scenario A")

ggplot(data=scenarioA %>% filter(Age=="All") %>% mutate(date=as.Date(date)))+
  geom_line(aes(x=date, y=value, group=sample))



#Averted just from monoclonals 
scenarioA_monoclonal = projection_function(RRHn1 = 0.01,
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
                                monoclonal_birth = opt_monoclonal_birth,
                                monoclonal_01 = opt_monoclonal_01,
                                monoclonal_23 = opt_monoclonal_23,
                                monoclonal_45 = opt_monoclonal_45,
                                monoclonal_67 = opt_monoclonal_67,
                                maternal_vax = counter_maternal_vax, 
                                senior_vax_75 = opt_sen_75,
                                senior_vax_60 = opt_sen_60_74,
                                parmset = parmset,
                                lhs_parms = fitted_params,
                                min_date = '2017-07-01',
                                max_date = '2025-04-28',
                                fit_times = fit_times) %>% 
  mutate(scenario = "Scenario A - monoclonal")


scenarioD_monoclonal = projection_function(RRHn1 = 0.01,
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
                                           monoclonal_birth = pes_monoclonal_birth,
                                           monoclonal_01 = pes_monoclonal_01,
                                           monoclonal_23 = pes_monoclonal_23,
                                           monoclonal_45 = pes_monoclonal_45,
                                           monoclonal_67 = pes_monoclonal_67,
                                           maternal_vax = counter_maternal_vax, 
                                           senior_vax_75 = pes_sen_75,
                                           senior_vax_60 = pes_sen_60_74,
                                           parmset = parmset,
                                           lhs_parms = fitted_params,
                                           min_date = '2017-07-01',
                                           max_date = '2025-04-28',
                                           fit_times = fit_times) %>% 
  mutate(scenario = "Scenario D - monoclonal")





#Averted just from maternal vax 
scenarioA_maternal = projection_function(RRHn1 = 0.01,
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
                                maternal_vax = opt_maternal_vax, 
                                senior_vax_75 = opt_sen_75,
                                senior_vax_60 = opt_sen_60_74,
                                parmset = parmset,
                                lhs_parms = fitted_params,
                                min_date = '2017-07-01',
                                max_date = '2025-04-28',
                                fit_times = fit_times) %>% 
  mutate(scenario = "Scenario A - maternal")


scenarioD_maternal = projection_function(RRHn1 = 0.01,
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
                                         maternal_vax = pes_maternal_vax, 
                                         senior_vax_75 = pes_sen_75,
                                         senior_vax_60 = pes_sen_60_74,
                                         parmset = parmset,
                                         lhs_parms = fitted_params,
                                         min_date = '2017-07-01',
                                         max_date = '2025-04-28',
                                         fit_times = fit_times) %>% 
  mutate(scenario = "Scenario D - maternal")



# Infant Sensitivity  -----------------------------------------------------

# Nirsevimab lasts 270 days 
scenarioA2 = projection_function(RRHn1 = 0.01,
                                RRHn2 = 0.01,
                                RRHv1 = 0.3,
                                RRHv2 = 0.3,
                                RRHs = 0.25,
                                RRIn = 1,
                                RRIv = 1,
                                RRIs =1,
                                waningN1 = 90,
                                waningN2 = 180,
                                waningV1 = 90,
                                waningV2 = 90,
                                waningS = 730.5,
                                monoclonal_birth = opt_monoclonal_birth,
                                monoclonal_01 = opt_monoclonal_01,
                                monoclonal_23 = opt_monoclonal_23,
                                monoclonal_45 = opt_monoclonal_45,
                                monoclonal_67 = opt_monoclonal_67,
                                maternal_vax = opt_maternal_vax, 
                                senior_vax_75 = opt_sen_75,
                                senior_vax_60 = opt_sen_60_74,
                                parmset = parmset,
                                lhs_parms = fitted_params,
                                min_date = '2017-07-01',
                                max_date = '2025-04-28',
                                fit_times = fit_times) %>% 
  mutate(scenario = "Scenario A2")

ggplot(data=scenarioA2 %>% filter(Age=="All") %>% mutate(date=as.Date(date)))+
  geom_line(aes(x=date, y=value, group=sample))



#Averted just from monoclonals 
scenarioA2_monoclonal = projection_function(RRHn1 = 0.01,
                                           RRHn2 = 0.01,
                                           RRHv1 = 0.3,
                                           RRHv2 = 0.3,
                                           RRHs = 0.25,
                                           RRIn = 1,
                                           RRIv = 1,
                                           RRIs =1,
                                           waningN1 = 90,
                                           waningN2 = 180,
                                           waningV1 = 90,
                                           waningV2 = 90,
                                           waningS = 730.5,
                                           monoclonal_birth = opt_monoclonal_birth,
                                           monoclonal_01 = opt_monoclonal_01,
                                           monoclonal_23 = opt_monoclonal_23,
                                           monoclonal_45 = opt_monoclonal_45,
                                           monoclonal_67 = opt_monoclonal_67,
                                           maternal_vax = counter_maternal_vax, 
                                           senior_vax_75 = opt_sen_75,
                                           senior_vax_60 = opt_sen_60_74,
                                           parmset = parmset,
                                           lhs_parms = fitted_params,
                                           min_date = '2017-07-01',
                                           max_date = '2025-04-28',
                                           fit_times = fit_times) %>% 
  mutate(scenario = "Scenario A2 - monoclonal")


# All catch-up given in October - November
early_monoclonal_01 = c(rep(0, length(dates2) + 104-31),rep(10200*.25/9,9),rep(0,22))
early_monoclonal_23 = c(rep(0, length(dates2) + 104-31),rep(10200*.25/9,9),rep(0,22))
early_monoclonal_45 = c(rep(0, length(dates2) + 104-31),rep(10200*.25/9,9),rep(0,22))
early_monoclonal_67 = c(rep(0, length(dates2) + 104-31),rep(10200*.25/9,9),rep(0,22))
scenarioA3 = projection_function(RRHn1 = 0.01,
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
                                monoclonal_birth = opt_monoclonal_birth,
                                monoclonal_01 = early_monoclonal_01,
                                monoclonal_23 = early_monoclonal_23,
                                monoclonal_45 = early_monoclonal_45,
                                monoclonal_67 = early_monoclonal_67,
                                maternal_vax = opt_maternal_vax, 
                                senior_vax_75 = opt_sen_75,
                                senior_vax_60 = opt_sen_60_74,
                                parmset = parmset,
                                lhs_parms = fitted_params,
                                min_date = '2017-07-01',
                                max_date = '2025-04-28',
                                fit_times = fit_times) %>% 
  mutate(scenario = "Scenario A3")

ggplot(data=scenarioA3 %>% filter(Age=="All") %>% mutate(date=as.Date(date)))+
  geom_line(aes(x=date, y=value, group=sample))

#Averted just from monoclonals 
scenarioA3_monoclonal = projection_function(RRHn1 = 0.01,
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
                                           monoclonal_birth = opt_monoclonal_birth,
                                           monoclonal_01 = early_monoclonal_01,
                                           monoclonal_23 = early_monoclonal_23,
                                           monoclonal_45 = early_monoclonal_45,
                                           monoclonal_67 = early_monoclonal_67,
                                           maternal_vax = counter_maternal_vax, 
                                           senior_vax_75 = opt_sen_75,
                                           senior_vax_60 = opt_sen_60_74,
                                           parmset = parmset,
                                           lhs_parms = fitted_params,
                                           min_date = '2017-07-01',
                                           max_date = '2025-04-28',
                                           fit_times = fit_times) %>% 
  mutate(scenario = "Scenario A3 - monoclonal")


# All birth doses are nirsevimab 
all_monoclonal_birth = c(rep(0, length(dates2) + 104-31),rep(10350/27,27),rep(0,4))
scenarioA4 = projection_function(RRHn1 = 0.01,
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
                                monoclonal_birth = all_monoclonal_birth,
                                monoclonal_01 = opt_monoclonal_01,
                                monoclonal_23 = opt_monoclonal_23,
                                monoclonal_45 = opt_monoclonal_45,
                                monoclonal_67 = opt_monoclonal_67,
                                maternal_vax = counter_maternal_vax, 
                                senior_vax_75 = opt_sen_75,
                                senior_vax_60 = opt_sen_60_74,
                                parmset = parmset,
                                lhs_parms = fitted_params,
                                min_date = '2017-07-01',
                                max_date = '2025-04-28',
                                fit_times = fit_times) %>% 
  mutate(scenario = "Scenario A4")

ggplot(data=scenarioA4 %>% filter(Age=="All") %>% mutate(date=as.Date(date)))+
  geom_line(aes(x=date, y=value, group=sample))

# All birth doses are maternal vaccination 
all_maternal_vax = c(rep(0, length(dates2) + 104-40),immu2$scale_mat*10350,rep(0,4))
scenarioA5 = projection_function(RRHn1 = 0.01,
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
                                monoclonal_01 = opt_monoclonal_01,
                                monoclonal_23 = opt_monoclonal_23,
                                monoclonal_45 = opt_monoclonal_45,
                                monoclonal_67 = opt_monoclonal_67,
                                maternal_vax = all_maternal_vax, 
                                senior_vax_75 = opt_sen_75,
                                senior_vax_60 = opt_sen_60_74,
                                parmset = parmset,
                                lhs_parms = fitted_params,
                                min_date = '2017-07-01',
                                max_date = '2025-04-28',
                                fit_times = fit_times) %>% 
  mutate(scenario = "Scenario A5")

ggplot(data=scenarioA5 %>% filter(Age=="All") %>% mutate(date=as.Date(date)))+
  geom_line(aes(x=date, y=value, group=sample))

#Averted just from maternal vax 
scenarioA5_maternal = projection_function(RRHn1 = 0.01,
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
                                         maternal_vax = all_maternal_vax, 
                                         senior_vax_75 = opt_sen_75,
                                         senior_vax_60 = opt_sen_60_74,
                                         parmset = parmset,
                                         lhs_parms = fitted_params,
                                         min_date = '2017-07-01',
                                         max_date = '2025-04-28',
                                         fit_times = fit_times) %>% 
  mutate(scenario = "Scenario A5 - maternal")


# Seniors Sensitivity  ----------------------------------------------------
pes2_sen_75 = c(rep(0, length(dates2) + 104-41),cum_sen_75*.87*.5,immu2$week_seniors_75*.87*.5,rep(0,4))
pes2_sen_60_74 = c(rep(0, length(dates2) + 104-41),cum_sen_60*.8*.5,immu2$week_seniors_60_74*.8*.5,rep(0,4))
# Assume that protection reduced by 50% in the second year 

waning_pessimistic = projection_function(RRHn1 = 0.01,
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
                                monoclonal_birth = pes_monoclonal_birth,
                                monoclonal_01 = pes_monoclonal_01,
                                monoclonal_23 = pes_monoclonal_23,
                                monoclonal_45 = pes_monoclonal_45,
                                monoclonal_67 = pes_monoclonal_67,
                                maternal_vax = pes_maternal_vax, 
                                senior_vax_75 = pes2_sen_75,
                                senior_vax_60 = pes2_sen_60_74,
                                parmset = parmset,
                                lhs_parms = fitted_params,
                                min_date = '2017-07-01',
                                max_date = '2025-04-28',
                                fit_times = fit_times) %>% 
  mutate(scenario = "Waning Pessimistic")

opt2_sen_75 = c(rep(0, length(dates2) + 104-41),cum_sen_75*.87*.5,immu2$week_seniors_75*.87,rep(0,4))
opt2_sen_60_74 = c(rep(0, length(dates2) + 104-41),cum_sen_60*.8*.5,immu2$week_seniors_60_74*.8,rep(0,4))
waning_optimistic = projection_function(RRHn1 = 0.01,
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
                                monoclonal_birth = opt_monoclonal_birth,
                                monoclonal_01 = opt_monoclonal_01,
                                monoclonal_23 = opt_monoclonal_23,
                                monoclonal_45 = opt_monoclonal_45,
                                monoclonal_67 = opt_monoclonal_67,
                                maternal_vax = opt_maternal_vax, 
                                senior_vax_75 = opt2_sen_75,
                                senior_vax_60 = opt2_sen_60_74,
                                parmset = parmset,
                                lhs_parms = fitted_params,
                                min_date = '2017-07-01',
                                max_date = '2025-04-28',
                                fit_times = fit_times) %>% 
  mutate(scenario = "Waning Optimistic")


# Assume that protection is combination of hospitalization and infection 
 
txn_pessimistic = projection_function(RRHn1 = 0.01,
                                         RRHn2 = 0.01,
                                         RRHv1 = 0.3,
                                         RRHv2 = 0.3,
                                         RRHs = 0.5,
                                         RRIn = 1,
                                         RRIv = 1,
                                         RRIs =0.5,
                                         waningN1 = 90,
                                         waningN2 = 90,
                                         waningV1 = 90,
                                         waningV2 = 90,
                                         waningS = 730.5,
                                         monoclonal_birth = pes_monoclonal_birth,
                                         monoclonal_01 = pes_monoclonal_01,
                                         monoclonal_23 = pes_monoclonal_23,
                                         monoclonal_45 = pes_monoclonal_45,
                                         monoclonal_67 = pes_monoclonal_67,
                                         maternal_vax = pes_maternal_vax, 
                                         senior_vax_75 = pes_sen_75,
                                         senior_vax_60 = pes_sen_60_74,
                                         parmset = parmset,
                                         lhs_parms = fitted_params,
                                         min_date = '2017-07-01',
                                         max_date = '2025-04-28',
                                         fit_times = fit_times) %>% 
  mutate(scenario = "TXN Pessimistic")

txn_optimistic = projection_function(RRHn1 = 0.01,
                                        RRHn2 = 0.01,
                                        RRHv1 = 0.3,
                                        RRHv2 = 0.3,
                                        RRHs = 0.5,
                                        RRIn = 1,
                                        RRIv = 1,
                                        RRIs =0.5,
                                        waningN1 = 90,
                                        waningN2 = 90,
                                        waningV1 = 90,
                                        waningV2 = 90,
                                        waningS = 730.5,
                                        monoclonal_birth = opt_monoclonal_birth,
                                        monoclonal_01 = opt_monoclonal_01,
                                        monoclonal_23 = opt_monoclonal_23,
                                        monoclonal_45 = opt_monoclonal_45,
                                        monoclonal_67 = opt_monoclonal_67,
                                        maternal_vax = opt_maternal_vax, 
                                        senior_vax_75 = opt_sen_75,
                                        senior_vax_60 = opt_sen_60_74,
                                        parmset = parmset,
                                        lhs_parms = fitted_params,
                                        min_date = '2017-07-01',
                                        max_date = '2025-04-28',
                                        fit_times = fit_times) %>% 
  mutate(scenario = "TXN Optimistic")
#protection lower in Year 2 

waning_txn_pessimistic = projection_function(RRHn1 = 0.01,
                                         RRHn2 = 0.01,
                                         RRHv1 = 0.3,
                                         RRHv2 = 0.3,
                                         RRHs = 0.5,
                                         RRIn = 1,
                                         RRIv = 1,
                                         RRIs =0.5,
                                         waningN1 = 90,
                                         waningN2 = 90,
                                         waningV1 = 90,
                                         waningV2 = 90,
                                         waningS = 730.5,
                                         monoclonal_birth = pes_monoclonal_birth,
                                         monoclonal_01 = pes_monoclonal_01,
                                         monoclonal_23 = pes_monoclonal_23,
                                         monoclonal_45 = pes_monoclonal_45,
                                         monoclonal_67 = pes_monoclonal_67,
                                         maternal_vax = pes_maternal_vax, 
                                         senior_vax_75 = pes2_sen_75,
                                         senior_vax_60 = pes2_sen_60_74,
                                         parmset = parmset,
                                         lhs_parms = fitted_params,
                                         min_date = '2017-07-01',
                                         max_date = '2025-04-28',
                                         fit_times = fit_times) %>% 
  mutate(scenario = "Waning TXN Pessimistic")


waning_txn_optimistic = projection_function(RRHn1 = 0.01,
                                        RRHn2 = 0.01,
                                        RRHv1 = 0.3,
                                        RRHv2 = 0.3,
                                        RRHs = 0.5,
                                        RRIn = 1,
                                        RRIv = 1,
                                        RRIs =0.5,
                                        waningN1 = 90,
                                        waningN2 = 90,
                                        waningV1 = 90,
                                        waningV2 = 90,
                                        waningS = 730.5,
                                        monoclonal_birth = opt_monoclonal_birth,
                                        monoclonal_01 = opt_monoclonal_01,
                                        monoclonal_23 = opt_monoclonal_23,
                                        monoclonal_45 = opt_monoclonal_45,
                                        monoclonal_67 = opt_monoclonal_67,
                                        maternal_vax = opt_maternal_vax, 
                                        senior_vax_75 = opt2_sen_75,
                                        senior_vax_60 = opt2_sen_60_74,
                                        parmset = parmset,
                                        lhs_parms = fitted_params,
                                        min_date = '2017-07-01',
                                        max_date = '2025-04-28',
                                        fit_times = fit_times) %>% 
  mutate(scenario = "Waning TXN Optimistic")


results  = rbind(counterfactual, scenarioA, scenarioB, scenarioC,scenarioD,
                 scenarioA_maternal,scenarioA_monoclonal,scenarioD_maternal,scenarioD_monoclonal,
                 scenarioA2, scenarioA2_monoclonal, scenarioA3, scenarioA3_monoclonal,
                 scenarioA4, scenarioA5, scenarioA5_maternal,
                 waning_optimistic, waning_pessimistic,
                 txn_optimistic, txn_pessimistic, 
                 waning_txn_optimistic, waning_txn_pessimistic)
saveRDS(results, "RESULTS/results_24-25_100replicates.rds")

