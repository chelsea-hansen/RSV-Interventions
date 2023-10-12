
rm(list=ls())
library(deSolve)
library(RColorBrewer)
library(reshape2)
library(cdcfluview)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(sjstats)
library(lubridate)
library(imputeTS)
library(zoo)
library(lhs)
library(plotrix)

"%notin%" = Negate('%in%')


#This part is similar to the first part of the model calibration code
## reading in and preparing data 

yinit <- readRDS('Data/yinit.rds') %>% 
  mutate(Mn=0,N=0,Vs1=0,Vs2=0) #add intervention compartments 
yinit=as.matrix(yinit)# set for 2003 population size 
yinit[,1:2]=yinit[,1:2]*0.7 # set for 1980 population size instead 
p <- sum(yinit) 


birth <-readRDS('Data/birth.rds') %>% select(-date, -period)
birth_new = birth[2284:2335,]

#extend dates into 2024-25 season, copy birthrates from last year of data 
birth = as.matrix(rbind(birth,birth_new))

#vector of dates from burn-in to end of 2024-25 season 
dates = data.frame(date=seq(from=as.Date("1980-01-05"), to=as.Date("2025-10-01"), by="weeks"))

#pre-pandemic fitting seasons are line 1934-2126

contact <- readRDS('Data/contact_POLYMOD.rds')

N_ages <- nrow(yinit) 
agenames <- c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y","2-4Y","5-9Y","10-19Y","20-39Y","40-64Y","65Y+")
al <- N_ages


rownames(yinit) <- agenames
yinit.vector <- as.vector(yinit)

name.array <- array(NA, dim=dim(yinit))
for(i in 1:dim(name.array)[1]){
  for(j in 1:dim(name.array)[2]){
    name.array[i,j] <- paste(dimnames(yinit)[[1]][i],dimnames(yinit)[[2]][j]  )
  }
}

name.vector <- as.vector(name.array)
names(yinit.vector) <- name.vector

#ageing 
WidthAgeClassMonth = c(rep(2,times=6), 12,12*3,  60, 120, 240, 300, 180)  #Aging rate=1/width age class (months)

#annual birth rates (converted to weeks in ODE code)
PerCapitaBirthsYear=birth 

#time steps for ODE 
tmax=2387
wa_times <- seq(1, tmax, by =1) 

#upload parameters fit by MLE 
fitLL = readRDS("parameters_21Sep2023.rds")

baseline.txn.rate=exp(fitLL$par[1])
b1=exp(fitLL$par[2])
phi=(2*pi*(exp(fitLL$par[3]))) / (1+exp(fitLL$par[3]))
DurationMatImmunityDays = exp(fitLL$par[4])
report_seniors <-  exp(fitLL$par[5])
report_infants <-  exp(fitLL$par[6])
I_ex <-  exp(fitLL$par[7])
RRHm <-  exp(fitLL$par[8])

#make vector for seasonal catch-up doses 

## Read in transmission dynamic model
source("Interventions/intervention_models.R")


# function to fit intervention scenarios (no projection intervals) ------------------------------------------

interventions = function(birth_dose, #timing and coverage of nirsevimab birth doses
                         cover_n, #timing and coverage of nirsevimab catch-up doses
                         waningN, #duration of nirsevimab protection (in days)
                         RRIn, #relative risk of infection while protected by nirsevimab (default = 1)
                         RRHn,#relative risk of hospitalization with protected bynirsevimab
                         cover_s, #timing and coverage of senior vaccination 
                         waningS, #duration of protection from vaccines (in days) for seniors
                         RRIs, #relative risk of infection after vaccination in seniors (default = 1)
                         RRHs #relative risk of hospitalization after vaccination in seniors
                         ){
  
  #parameters for fitting the pandemic and rebound periods 
  #suppression of seeding
  introductions = data.frame(intros=c(rep(I_ex,2100),rep(0,44),rep(NA,12),rep(I_ex,231))) %>% 
    mutate(intros = na_interpolation(intros, method="linear"))
  introductions = introductions$intros
  #suppression of beta (due to NPIs)
  npi = data.frame(npis=c(rep(1,2100),rep(.75,13),rep(NA,35),rep(1,40),rep(.75,12),rep(1,187)))%>% 
    mutate(npis= na_interpolation(npis, method="linear"))
  npi = npi$npis
  
  parmset<-list(PerCapitaBirthsYear=PerCapitaBirthsYear,
                WidthAgeClassMonth=WidthAgeClassMonth,
                baseline.txn.rate = baseline.txn.rate,
                phi=phi,
                b1=b1,
                DurationMatImmunityDays=DurationMatImmunityDays,
                report_seniors=report_seniors,
                report_infants=report_infants,
                um=-0.0002227,
                rho1=0.75,
                rho2=0.51,
                dur.days1=10,
                dur.days2=7,
                dur.days3=5,
                yinit.matrix=yinit,
                q=1,
                c2=contact,
                sigma1=.76,
                sigma2=.6,
                sigma3=.4,
                birth_dose = birth_dose,
                cover_n = cover_n,
                waningN=waningN,
                RRIn=RRIn,
                RRHn=RRHn,
                RRHm=RRHm,
                cover_s = cover_s,
                waningS=waningS,
                RRIs = RRIs,
                RRHs = RRHs,
                npi=npi,
                introductions = introductions,
                time.step='week'
  )
  
  output <- ode(y=yinit.vector, t=wa_times,
                method = "ode45",
                func=intervention_models, 
                parms=c(parmset))
  
  
  t0 <- 455 #just save results after burn-in period
  al <- nrow(yinit)
  output <- tail(output,t0)
  St <- output[,-1]
  
  I1 <- St[,grep('I1', colnames(St))]
  I2 <- St[,grep('I2', colnames(St))]
  I3 <- St[,grep('I3', colnames(St))]
  I4 <- St[,grep('I4', colnames(St))]
  S0 <- St[,grep('S0', colnames(St))]
  S1 <- St[,grep('S1', colnames(St))]
  S2 <- St[,grep('S2', colnames(St))]
  S3 <- St[,grep('S3', colnames(St))]
  
  M <- St[,grep('M', colnames(St))]
  Mn <- St[,grep('Mn', colnames(St))]
  N <- St[,grep('N', colnames(St))]
  
  Vs1 <- St[,grep('Vs1', colnames(St))]
  Vs2 <- St[,grep('Vs2', colnames(St))]
  
  contact2 = npi[1933:2387]
  intro2 = introductions[1933:2387]
  
  lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
  for (t in 1:t0){
    beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(parmset$yinit.matrix)^(1-parmset$q))*parmset$c2*contact2[t]
    lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+intro2[t])%*%beta)/sum(St[t,]))}
  
  hosp = c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2,report_infants*0.15, report_infants*0.15, report_seniors,report_seniors,report_seniors*0.06,report_seniors*0.06,report_seniors*0.06,report_seniors*0.06,report_seniors)
  
  H1=matrix(0,nrow=t0,ncol=al)#Number of hospitalizations by age
  for (i in 1:al){
    H1[,i]=RRHm*hosp[i]*M[,i]*lambda1[,i]+
      RRHn*((RRHm*2)/3)*hosp[i]*RRIn*Mn[,i]*lambda1[,i]+ #note, maternal protection is weighted to reflect that 1/3 of infants in this compartment are still protected by maternal immunity (nirsevimab is added to this protection)
      RRHn*hosp[i]*RRIn*N[,i]*lambda1[,i]+  
      hosp[i]*S0[,i]*lambda1[,i]+
      hosp[i]*parmset$sigma1*S1[,i]*lambda1[,i]+
      hosp[i]*parmset$sigma2*S2[,i]*lambda1[,i]+
      hosp[i]*parmset$sigma3*S3[,i]*lambda1[,i]+
      RRHs*hosp[i]*parmset$RRIs*parmset$sigma3*Vs1[,i]*lambda1[,i]+
      RRHs*hosp[i]*parmset$RRIs*parmset$sigma3*Vs2[,i]*lambda1[,i]
  }
  
  H2=cbind(H1[,1],H1[,2],H1[,3],H1[,4],H1[,5],H1[,6],rowSums(H1[,7:8]),rowSums(H1[,9:12]),H1[,13])

  
  all_ages = data.frame(hosp = rowSums(H2[352:455,]), 
                        date = dates$date[2284:2387], 
                        age="All")# total hospitalizations
  infants = data.frame(hosp = rowSums(H2[352:455,1:6]),
                       date = dates$date[2284:2387], 
                       age="Infants <1y")# intervention given to infants <8m but benefit is intended for all infants 
  seniors = data.frame(hosp = H2[352:455,9],
                       date = dates$date[2284:2387], 
                       age="Seniors")# seniors
 
  results = rbind(all_ages, infants, seniors) %>%
    mutate(season = ifelse(date>='2024-10-01',"season24_25","season23_24"))
  
  return(results)
}


# function to fit intervention scenarios (with projection intervals)------------------------------------------
# this function uses lating hypercube sampling to get the median and 95% projection intervals 
# this version takes much longer to run than the version above 

interventions_PI = function(birth_dose, cover_n, waningN,RRIn, RRHn,
                         cover_s, waningS, RRIs, RRHs){
  
  #parameters for fitting the pandemic and rebound periods 
  introductions = data.frame(intros=c(rep(I_ex,2100),rep(0,44),rep(NA,12),rep(I_ex,231))) %>% 
    mutate(intros = na_interpolation(intros, method="linear"))
  introductions = introductions$intros
  
  npi = data.frame(npis=c(rep(1,2100),rep(.75,13),rep(NA,35),rep(1,40),rep(.75,12),rep(1,187)))%>% 
    mutate(npis= na_interpolation(npis, method="linear"))
  npi = npi$npis
  
  parmset<-list(PerCapitaBirthsYear=PerCapitaBirthsYear,
                WidthAgeClassMonth=WidthAgeClassMonth,
                baseline.txn.rate = baseline.txn.rate,
                phi=phi,
                b1=b1,
                DurationMatImmunityDays=DurationMatImmunityDays,
                report_seniors=report_seniors,
                report_infants=report_infants,
                um=-0.0002227,
                rho1=0.75,
                rho2=0.51,
                dur.days1=10,
                dur.days2=7,
                dur.days3=5,
                yinit.matrix=yinit,
                q=1,
                c2=contact,
                sigma1=.76,
                sigma2=.6,
                sigma3=.4,
                birth_dose = birth_dose,
                cover_n = cover_n,
                waningN=waningN,
                RRIn=RRIn,
                RRHn=RRHn,
                RRHm=RRHm,
                cover_s = cover_s,
                waningS=waningS,
                RRIs = RRIs,
                RRHs = RRHs,
                npi=npi,
                introductions = introductions,
                time.step='week'
  )
  
    output <- ode(y=yinit.vector, t=wa_times,
                  method = "ode45",
                  func=intervention_models, 
                  parms=c(parmset))
    
    
    t0 <- 455
    al <- nrow(yinit)
    output <- tail(output,t0)
    St <- output[,-1]
    
    I1 <- St[,grep('I1', colnames(St))]
    I2 <- St[,grep('I2', colnames(St))]
    I3 <- St[,grep('I3', colnames(St))]
    I4 <- St[,grep('I4', colnames(St))]
    S0 <- St[,grep('S0', colnames(St))]
    S1 <- St[,grep('S1', colnames(St))]
    S2 <- St[,grep('S2', colnames(St))]
    S3 <- St[,grep('S3', colnames(St))]
    
    M <- St[,grep('M', colnames(St))]
    Mn <- St[,grep('Mn', colnames(St))]
    N <- St[,grep('N', colnames(St))]
    
    Vs1 <- St[,grep('Vs1', colnames(St))]
    Vs2 <- St[,grep('Vs2', colnames(St))]
    
    
    contact2 = npi[1933:2387]
    intro2 = introductions[1933:2387]
    
    lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
    for (t in 1:t0){
      beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(parmset$yinit.matrix)^(1-parmset$q))*parmset$c2*contact2[t]
      lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+intro2[t])%*%beta)/sum(St[t,]))}
    
    hosp = c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2,report_infants*0.15, report_infants*0.15, report_seniors,report_seniors,report_seniors*0.06,report_seniors*0.06,report_seniors*0.06,report_seniors*0.06,report_seniors)
    
    H1=matrix(0,nrow=t0,ncol=al)#Number of hospitalizations by age
    for (i in 1:al){
      H1[,i]=RRHm*hosp[i]*M[,i]*lambda1[,i]+
        RRHn*((RRHm+2)/3)*hosp[i]*RRIn*Mn[,i]*lambda1[,i]+#note, maternal protection is weighted to reflect that 1/3 of infants in this compartment are still protected by maternal immunity (nirsevimab is added to this protection)
        RRHn*hosp[i]*RRIn*N[,i]*lambda1[,i]+
        hosp[i]*S0[,i]*lambda1[,i]+
        hosp[i]*parmset$sigma1*S1[,i]*lambda1[,i]+
        hosp[i]*parmset$sigma2*S2[,i]*lambda1[,i]+
        hosp[i]*parmset$sigma3*S3[,i]*lambda1[,i]+
        RRHs*hosp[i]*parmset$RRIs*parmset$sigma3*Vs1[,i]*lambda1[,i]+
        RRHs*hosp[i]*parmset$RRIs*parmset$sigma3*Vs2[,i]*lambda1[,i]
    }
    
    H2=cbind(H1[,1],H1[,2],H1[,3],H1[,4],H1[,5],H1[,6],rowSums(H1[,7:8]),rowSums(H1[,9:12]),H1[,13])
    
    all_ages = data.frame(hosp = rowSums(H2[352:455,]), 
                          date = dates$date[2284:2387], 
                          age="All")# total hospitalizations
    infants = data.frame(hosp = rowSums(H2[352:455,1:6]),
                         date = dates$date[2284:2387], 
                         age="Infants <1y")# intervention given to infants <8m but benefit is intended for all infants 
    seniors = data.frame(hosp = H2[352:455,9],
                         date = dates$date[2284:2387], 
                         age="Seniors")# seniors
    
    projection_intervals = rbind(all_ages, infants, seniors) %>%
      mutate(season = ifelse(date>='2024-10-01',"season24_25","season23_24"),
             lhs = "point estimate")
   
    
   
  #create easy confidence intervals (+/- 5%)
  # does not include variables that scale infections to hospitalizations 
  beta_lower = baseline.txn.rate*.95
  beta_upper = baseline.txn.rate*1.05
  b1_lower = b1*.95
  b1_upper = b1*1.05
  phi_lower = phi*.95
  phi_upper = phi*1.05
  duration_lower = 28
  duration_upper = 112
  I_ex_lower = I_ex*.95
  I_ex_upper = I_ex*1.05
  
  #set-up LHS
  h <- 1000         
  set.seed(123)
  lhs<-maximinLHS(h,5)
  
  new_parms <- cbind(
    baseline.txn.rate = lhs[,1]*(beta_upper-beta_lower)+beta_lower,
    b1 = lhs[,2]*(b1_upper-b1_lower)+b1_lower,
    phi = lhs[,3]*(phi_upper-phi_lower)+phi_lower,
    DurationMatImmunityDays = lhs[,4]*(duration_upper-duration_lower)+duration_lower,
    I_ex = lhs[,5]*(I_ex_upper-I_ex_lower)+I_ex_lower)
  
  #parameters for fitting the pandemic and rebound periods 
  for(l in 1:1000){
   introductions = data.frame(intros=c(rep(new_parms[l,"I_ex"],2100),rep(0,44),rep(NA,12),rep(new_parms[l,"I_ex"],231))) %>% 
    mutate(intros = na_interpolation(intros, method="linear"))
  introductions = introductions$intros
  
  parmset<-list(PerCapitaBirthsYear=PerCapitaBirthsYear,
              WidthAgeClassMonth=WidthAgeClassMonth,
              baseline.txn.rate = new_parms[l,"baseline.txn.rate"],
              phi=new_parms[l,"phi"],
              b1=new_parms[l,"b1"],
              DurationMatImmunityDays=new_parms[l,"DurationMatImmunityDays"],
              report_seniors=report_seniors,
              report_infants=report_infants,
              um=-0.0002227,
              rho1=0.75,
              rho2=0.51,
              dur.days1=10,
              dur.days2=7,
              dur.days3=5,
              yinit.matrix=yinit,
              q=1,
              c2=contact,
              sigma1=.76,
              sigma2=.6,
              sigma3=.4,
              birth_dose = birth_dose,
              cover_n = cover_n,
              waningN=waningN,
              RRIn=RRIn,
              RRHn=RRHn,
              RRHm=RRHm,
              cover_s = cover_s,
              waningS=waningS,
              RRIs = RRIs,
              RRHs = RRHs,
              npi=npi,
              introductions = introductions,
              time.step='week'
)

  
output <- ode(y=yinit.vector, t=wa_times,
              method = "ode45",
              func=intervention_models, 
              parms=c(parmset))


t0 <- 455
al <- nrow(yinit)
output <- tail(output,t0)
St <- output[,-1]

I1 <- St[,grep('I1', colnames(St))]
I2 <- St[,grep('I2', colnames(St))]
I3 <- St[,grep('I3', colnames(St))]
I4 <- St[,grep('I4', colnames(St))]
S0 <- St[,grep('S0', colnames(St))]
S1 <- St[,grep('S1', colnames(St))]
S2 <- St[,grep('S2', colnames(St))]
S3 <- St[,grep('S3', colnames(St))]

M <- St[,grep('M', colnames(St))]
Mn <- St[,grep('Mn', colnames(St))]
N <- St[,grep('N', colnames(St))]

Vs1 <- St[,grep('Vs1', colnames(St))]
Vs2 <- St[,grep('Vs2', colnames(St))]


contact2 = npi[1933:2387]
intro2 = introductions[1933:2387]

lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
for (t in 1:t0){
  beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(parmset$yinit.matrix)^(1-parmset$q))*parmset$c2*contact2[t]
  lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+intro2[t])%*%beta)/sum(St[t,]))}

hosp = c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2,report_infants*0.15, report_infants*0.15, report_seniors,report_seniors,report_seniors*0.06,report_seniors*0.06,report_seniors*0.06,report_seniors*0.06,report_seniors)

H1=matrix(0,nrow=t0,ncol=al)#Number of hospitalizations by age
for (i in 1:al){
  H1[,i]=RRHm*hosp[i]*M[,i]*lambda1[,i]+
    RRHn*((RRHm+2)/3)*hosp[i]*RRIn*Mn[,i]*lambda1[,i]+#note, maternal protection is weighted to reflect that 1/3 of infants in this compartment are still protected by maternal immunity (nirsevimab is added to this protection)
    RRHn*hosp[i]*RRIn*N[,i]*lambda1[,i]+
    hosp[i]*S0[,i]*lambda1[,i]+
    hosp[i]*parmset$sigma1*S1[,i]*lambda1[,i]+
    hosp[i]*parmset$sigma2*S2[,i]*lambda1[,i]+
    hosp[i]*parmset$sigma3*S3[,i]*lambda1[,i]+
    RRHs*hosp[i]*parmset$RRIs*parmset$sigma3*Vs1[,i]*lambda1[,i]+
    RRHs*hosp[i]*parmset$RRIs*parmset$sigma3*Vs2[,i]*lambda1[,i]
  }

H2=cbind(H1[,1],H1[,2],H1[,3],H1[,4],H1[,5],H1[,6],rowSums(H1[,7:8]),rowSums(H1[,9:12]),H1[,13])

all_ages = data.frame(hosp = rowSums(H2[352:455,]), 
                      date = dates$date[2284:2387], 
                      age="All")# total hospitalizations
infants = data.frame(hosp = rowSums(H2[352:455,1:6]),
                     date = dates$date[2284:2387], 
                     age="Infants <1y")# intervention given to infants <8m but benefit is intended for all infants 
seniors = data.frame(hosp = H2[352:455,9],
                     date = dates$date[2284:2387], 
                     age="Seniors")# seniors

estimates = rbind(all_ages, infants, seniors) %>%
  mutate(season = ifelse(date>='2024-10-01',"season24_25","season23_24"),lhs = l)

projection_intervals = rbind(projection_intervals, estimates)

  } 
  return(projection_intervals)
  }




# flu curves for  to base timing of vaccination (from Flu Scenario Modeling Hub) --------------------
# flu curves give timing, scale them down to match the level of coverage for scenarios 

#timing of senior vaccines 
#76% cumulative coverage
senior_vaccination = read_csv("https://raw.githubusercontent.com/midas-network/flu-scenario-modeling-hub_resources/main/Rd2_datasets/Age_Specific_Coverage_Flu_RD2_2022_23_Sc_A_B_C_D.csv") %>% 
  filter(Geography=="Washington", Age=="65+ Years") %>% 
  select("date"=Week_Ending_Sat, "flu_coverage"=flu.coverage.rd2.sc_A_B_C_D) %>% 
  mutate(flu_cov_cum = flu_coverage/100,
         flu_cov_week = flu_cov_cum - lag(flu_cov_cum,1)) %>% 
  replace(is.na(.),0)

#timing of vaccination for catch-up vaccines (based on flu in 6mo-4yrs)
#75% cumulative coverage 
child_vaccination = read_csv("https://raw.githubusercontent.com/midas-network/flu-scenario-modeling-hub_resources/main/Rd2_datasets/Age_Specific_Coverage_Flu_RD2_2022_23_Sc_A_B_C_D.csv") %>% 
  filter(Geography=="Washington", Age=="6 Months - 4 Years") %>% 
  select("date"=Week_Ending_Sat, "flu_coverage"=flu.coverage.rd2.sc_A_B_C_D) %>% 
  mutate(flu_cov_cum = flu_coverage/100,
         flu_cov_week = flu_cov_cum - lag(flu_cov_cum,1),
         flu_lag_cum = lag(flu_cov_cum,7), #start coverage first week of October 
         flu_lag_week = lag(flu_cov_week,7)) %>% 
  replace(is.na(.),0)


# No Interventions - baseline scenario  --------------------------------------
#coverage is zero for full time period
no_interventions = interventions(birth_dose=c(rep(0,2387)), 
                                 cover_n=c(rep(0,2387)), 
                                 waningN=36500,
                                 RRIn=1, 
                                 RRHn=1, 
                                 cover_s=c(rep(0,2387)), 
                                 waningS=36500, 
                                 RRIs=1, 
                                 RRHs=1)

no_intervention_totals = no_interventions %>% 
  group_by(age, season) %>% 
  summarize(hosp = sum(hosp))

no_interventions_PI = interventions_PI(birth_dose=c(rep(0,2387)), 
                                 cover_n=c(rep(0,2387)), 
                                 waningN=36500,
                                 RRIn=1, 
                                 RRHn=1, 
                                 cover_s=c(rep(0,2387)), 
                                 waningS=36500, 
                                 RRIs=1, 
                                 RRHs=1)


no_interventions_timeseries = no_interventions_PI %>% filter(lhs!="point estimate") %>% 
  group_by(age, date) %>% 
  summarize(median = median(hosp),
         lower = quantile(hosp,probs=0.025),
         upper = quantile(hosp,probs=0.975))

no_interventions_sum = no_interventions_PI %>% filter(lhs!="point estimate") %>% 
  group_by(age, season, lhs) %>% 
  summarize(hosp = sum(hosp)) %>% 
  ungroup() %>% 
  group_by(age, season) %>% 
  summarize(median = median(hosp),
            lower = quantile(hosp,probs=0.025),
            upper = quantile(hosp,probs=0.975))

# realistic/pessimistic interventions  ------------------------------------------------
#senior vaccination = 15% 
#senior vaccine effectiveness = 70% (RRHs = .3)


#nirsevimab coverage = 50% (2023-24); 75% (2024-25)
#nirsevimab effectiveness = 60% (RRHn=.4)
#timing of birth doses = October to March (6 months)


#scaling 75% (flu coverage) to 50% for 2023-24 (*.67)
#scaling 75% (flu coverage) to 75% for 2024-25 (*1)
catch_up.p = c(rep(0,2274),child_vaccination$flu_lag_week[1:44]*.67,rep(0,8),child_vaccination$flu_cov_week[1:44]*1,rep(0,8),child_vaccination$flu_cov_week[1:9]*1)
catch_up.p_cum1 = cumsum(catch_up.p[2275:2327])
plot(catch_up.p_cum1) #season 1
catch_up.p_cum2 = cumsum(catch_up.p[2328:2379])
plot(catch_up.p_cum2) #season 2

#assumes coverage is the same every week (no ramp-up), administered October to March
#alternative is to use the same coverage as the catch-up doses 
birth_dose.p =c(rep(0,2283),rep(.5,26),rep(0,26),rep(.75,26),rep(0,26))

#scaling 76% to 15% in both seasons (*.2)
senior_vax.p = c(rep(0,2274),senior_vaccination$flu_cov_week*.2,rep(0,8),senior_vaccination$flu_cov_week*.2,rep(0,8),senior_vaccination$flu_cov_week[1:9]*.2)
senior_cum.p = cumsum(senior_vax.p[2275:2387])
plot(senior_cum.p)


realistic_interventions = interventions(birth_dose= birth_dose.p,
                                        cover_n=catch_up.p, 
                                        waningN=180,
                                        RRIn=1, 
                                        RRHn=.4, #60% VE
                                        cover_s=senior_vax.p, 
                                        waningS=730.5, 
                                        RRIs=1, 
                                        RRHs=.3) #70% VE
  
  
 
 realistic_intervention_totals = realistic_interventions %>% 
   group_by(age, season) %>% 
   summarize(hosp = sum(hosp))
 
 
 realistic_interventions_PI = interventions_PI(birth_dose=birth_dose.p,
                                         cover_n=catch_up.p, 
                                         waningN=180,
                                         RRIn=1, 
                                         RRHn=.4, 
                                         cover_s=senior_vax.p, 
                                         waningS=730.5, 
                                         RRIs=1, 
                                         RRHs=.3)


  realistic_interventions_timeseries = realistic_interventions_PI %>% filter(lhs!="point estimate") %>% 
    group_by(age, date) %>% 
    summarize(median = median(hosp),
              lower = quantile(hosp,probs=0.025),
              upper = quantile(hosp,probs=0.975))
  
  realistic_interventions_sum = realistic_interventions_PI %>% filter(lhs!="point estimate") %>% 
    group_by(age, season, lhs) %>% 
    summarize(hosp = sum(hosp)) %>% 
    ungroup() %>% 
    group_by(age, season) %>% 
    summarize(median = median(hosp),
              lower = quantile(hosp,probs=0.025),
              upper = quantile(hosp,probs=0.975))
  
  

# Optimistic scenarios ----------------------------------------------------
 #senior vaccination = 30% 
 #senior vaccine effectiveness = 90% 
 
 #nirsevimab coverage = 75% (2023-24); 90% (2024-25)
 #nirsevimab effectiveness = 80% 
 #timing of birth doses = October to March (6 months )
 
  #scale 75% (flu coverage) to 75% for 2023-2024 (*1)
  #scale 75% (flu coverage) to 90% for 2024-2025 (*1.2)
  catch_up.o = c(rep(0,2274),child_vaccination$flu_lag_week[1:44]*1,rep(0,8),child_vaccination$flu_cov_week[1:44]*1.2,rep(0,8),child_vaccination$flu_cov_week[1:9]*1.2)
  catch_up.o_cum1 = cumsum(catch_up.o[2275:2327])
  plot(catch_up.o_cum1) #season 1
  catch_up.o_cum2 = cumsum(catch_up.o[2328:2379])
  plot(catch_up.o_cum2) #season 2
  
  #assumes coverage is the same every week (no ramp-up), administered October to March 
  #alternative is to use the same coverage as the catch-up doses 
  birth_dose.o =c(rep(0,2283),rep(.75,26),rep(0,26),rep(.9,26),rep(0,26))
  
  
  #scaling 76% (flu coverage) to 30% in both seasons (*.4)
  senior_vax.o = c(rep(0,2274),senior_vaccination$flu_cov_week*.4,rep(0,8),senior_vaccination$flu_cov_week*.4,rep(0,8),senior_vaccination$flu_cov_week[1:9]*.4)
  senior_cum.o = cumsum(senior_vax.o[2275:2387])
  plot(senior_cum.o)
  
  
  
  optimistic_interventions = interventions(birth_dose=birth_dose.o,
                                          cover_n=catch_up.o, 
                                          waningN=180,
                                          RRIn=1, 
                                          RRHn=.2, 
                                          cover_s=senior_vax.o, 
                                          waningS=730.5, 
                                          RRIs=1, 
                                          RRHs=.1)
  
  
  optimistic_intervention_totals = optimistic_interventions %>% 
    group_by(age, season) %>% 
    summarize(hosp = sum(hosp))
  
  optimistic_interventions_PI = interventions_PI(birth_dose=birth_dose.o,
                                                cover_n=catch_up.o, 
                                                waningN=180,
                                                RRIn=1, 
                                                RRHn=.2, 
                                                cover_s=senior_vax.o, 
                                                waningS=730.5, 
                                                RRIs=1, 
                                                RRHs=.1)
  
  
 
  optimistic_interventions_timeseries = optimistic_interventions_PI %>% filter(lhs!="point estimate") %>% 
    group_by(age, date) %>% 
    summarize(median = median(hosp),
              lower = quantile(hosp,probs=0.025),
              upper = quantile(hosp,probs=0.975))
  
  optimistic_interventions_sum = optimistic_interventions_PI %>% filter(lhs!="point estimate") %>% 
    group_by(age, season, lhs) %>% 
    summarize(hosp = sum(hosp)) %>% 
    ungroup() %>% 
    group_by(age, season) %>% 
    summarize(median = median(hosp),
              lower = quantile(hosp,probs=0.025),
              upper = quantile(hosp,probs=0.975))

# Optimistic Scenario with different timing  ------------------------------
 
# Combine into single files  ----------------------------------------------

point_estimates = rbind(no_interventions %>% mutate(interventions = "None"),
                        realistic_interventions %>% mutate(interventions = "Realistic"),
                        optimistic_interventions %>% mutate(interventions = "Optimistic"))
 write.csv(point_estimates, "point_estimates.csv")

 
point_estimates_totals = rbind(no_intervention_totals %>% mutate(interventions = "None"),
                          realistic_intervention_totals %>% mutate(interventions = "Realistic"),
                          optimistic_intervention_totals %>% mutate(interventions = "Optimistic"))
  write.csv(point_estimates_totals, "point_estimates_totals.csv")
  
  projection_intervals = rbind(no_interventions_PI %>% mutate(interventions = "None"),
                          realistic_interventions_PI %>% mutate(interventions = "Realistic"),
                          optimistic_interventions_PI %>% mutate(interventions = "Optimistic"))
  write.csv(projection_intervals, "projection_intervals.csv")
  
  projection_intervals_timeseries = rbind(no_interventions_timeseries %>% mutate(interventions = "None"),
                               realistic_interventions_timeseries %>% mutate(interventions = "Realistic"),
                               optimistic_interventions_timeseries %>% mutate(interventions = "Optimistic"))
  
  write.csv(projection_intrevals_timeseries, "projection_intervals_timeseries.csv")
  
  projection_intervals_sum = rbind(no_interventions_sum %>% mutate(interventions = "None"),
                                          realistic_interventions_sum %>% mutate(interventions = "Realistic"),
                                          optimistic_interventions_sum %>% mutate(interventions = "Optimistic"))
  write.csv(projection_intrevals_sum, "projection_intervals_sum.csv")
  