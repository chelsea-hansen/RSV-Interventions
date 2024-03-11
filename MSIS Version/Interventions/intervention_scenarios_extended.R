
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
      RRHn*RRHm*hosp[i]*RRIn*Mn[,i]*lambda1[,i]+ #note, this is adding nirsevimab impact on top of maternal immunity (might be too high?)
      #RRHn*hosp[i]*RRIn*Mn[,i]*lambda1[,i]+ #alternative approach is to add nirsevimab impact on it's own
      #RRHn*RRHm*hosp[i]*RRIn*N[,i]*lambda1[,i]+
      RRHn*hosp[i]*RRIn*N[,i]*lambda1[,i]+# see note above 
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
        RRHn*hosp[i]*RRIn*Mn[,i]*lambda1[,i]+
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
    RRHn*hosp[i]*RRIn*Mn[,i]*lambda1[,i]+
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



# No Interventions - baseline scenario  --------------------------------------
#coverage is zero for full time period
no_interventions = interventions_PI(birth_dose=c(rep(0,2387)), 
                                 cover_n=c(rep(0,2387)), 
                                 waningN=36500,
                                 RRIn=1, 
                                 RRHn=1, 
                                 cover_s=c(rep(0,2387)), 
                                 waningS=36500, 
                                 RRIs=1, 
                                 RRHs=1) %>% 
  mutate(coverage="None",
         VE="None")

write.csv(no_interventions, "no_interventions_full.rds")

# Vaccinations for seniors ------------------------------------------------

senior_vaccination = read_csv("https://raw.githubusercontent.com/midas-network/flu-scenario-modeling-hub_resources/main/Rd2_datasets/Age_Specific_Coverage_Flu_RD2_2022_23_Sc_A_B_C_D.csv") %>% 
  filter(Geography=="Washington", Age=="65+ Years") %>% 
  select("date"=Week_Ending_Sat, "flu_coverage"=flu.coverage.rd2.sc_A_B_C_D) %>% 
  mutate(flu_cov_cum = flu_coverage/100,
         flu_cov_week = flu_cov_cum - lag(flu_cov_cum,1)) %>% 
  replace(is.na(.),0)

#pessimistic/ realistic vaccination coverage (15%) 
senior_vax_pess = c(rep(0,2274),senior_vaccination$flu_cov_week*.2,rep(0,7),senior_vaccination$flu_cov_week*.2,rep(0,7),senior_vaccination$flu_cov_week[1:11]*.2)

#optimistic vaccination coverage (30%)
senior_vax_opt = c(rep(0,2274),senior_vaccination$flu_cov_week*.4,rep(0,7),senior_vaccination$flu_cov_week*.4,rep(0,7),senior_vaccination$flu_cov_week[1:11]*.4)


#pessimistic/ realistic vaccine effectiveness (70%)
VE_pessimistic = 0.3
#optimistic vaccine effectiveness (90%) 
VE_optimistic = 0.1

#vectorize scenairos 
senior_coverage = list("pessimistic"=senior_vax_pess, "optimistic"=senior_vax_opt)
senior_ve = c("pessimistic"=VE_pessimistic, "optimistic"=VE_optimistic)

#because we are assuming no impact on transmission, the infant parameters can remain equal to no interventions
senior_interventions = data.frame(hosp=NA, date=NA,age=NA,season=NA, lhs=NA,coverage=NA, VE=NA)

for(c in 1:2){
  for(v in 1:2){
    intervention_results = interventions_PI(birth_dose=c(rep(0,2387)), 
                                           cover_n=c(rep(0,2387)), 
                                           waningN=36500,
                                           RRIn=1, 
                                           RRHn=1, 
                                           cover_s=senior_coverage[[c]], 
                                           waningS=730.5, 
                                           RRIs=1, 
                                           RRHs=senior_ve[v]) %>% 
      filter(age=="Seniors") %>% 
      mutate(coverage = names(senior_coverage[c]),
             VE = names(senior_ve[v]))
    senior_interventions = rbind(senior_interventions, intervention_results)
  }
}

saveRDS(senior_interventions, "senior_interventions_full.rds")

# Infant interventions ----------------------------------------------------
#timing of vaccination for catch-up vaccines (based on flu in 6mo-4yrs)
#assuming catch-up doses end after December (#week 22)
#64% cumulative coverage by end of Dec - scale for RSV coverage scenarios
child_vaccination = read_csv("https://raw.githubusercontent.com/midas-network/flu-scenario-modeling-hub_resources/main/Rd2_datasets/Age_Specific_Coverage_Flu_RD2_2022_23_Sc_A_B_C_D.csv") %>% 
  filter(Geography=="Washington", Age=="6 Months - 4 Years") %>% 
  select("date"=Week_Ending_Sat, "flu_coverage"=flu.coverage.rd2.sc_A_B_C_D) %>% 
  mutate(flu_coverage = ifelse(date=="2022-08-13",0,flu_coverage),
         flu_cov_cum = flu_coverage/100,
         flu_cov_week = flu_cov_cum - lag(flu_cov_cum,1),
         week_lag = lag(flu_cov_week,6),
         cum_lag = lag(flu_coverage,6)) %>% 
  replace(is.na(.),0)


#(25%/50%)
child_pess = c(rep(0,2274),child_vaccination$week_lag[1:22]*.5,rep(0,30),child_vaccination$flu_cov_week[1:22]*.78,rep(0,30),child_vaccination$flu_cov_week[1:9]*.78)
birth_pess = c(rep(0,2283),rep(.25,26),rep(0,26),rep(.5,26),rep(0,26))
child_pess_cum = cumsum(child_pess[2327:2387])
plot(child_pess_cum)

#(50%/75%)
child_real = c(rep(0,2274),child_vaccination$week_lag[1:22]*.9,rep(0,30),child_vaccination$flu_cov_week[1:22]*1.17,rep(0,30),child_vaccination$flu_cov_week[1:9]*1.17)
birth_real = c(rep(0,2283),rep(.5,26),rep(0,26),rep(.75,26),rep(0,26))
child_real_cum = cumsum(child_real[2327:2387])
plot(child_real_cum)


#(75%/90%)
child_opt = c(rep(0,2274),child_vaccination$week_lag[1:22]*1.4,rep(0,30),child_vaccination$flu_cov_week[1:22]*1.4,rep(0,30),child_vaccination$flu_cov_week[1:9]*1.4)
birth_opt = c(rep(0,2283),rep(.75,26),rep(0,26),rep(.9,26),rep(0,26))
child_opt_cum = cumsum(child_opt[2327:2379])
plot(child_opt_cum)

#pessimistic/ realistic vaccine effectiveness (60%)
VE_pessimistic_infant = 0.4
#optimistic vaccine effectiveness (80%) 
VE_optimistic_infant = 0.2

#vectorize scenarios
infant_coverage = list("pessimistic"=child_pess, "realistic"=child_real, "optimistic"=child_opt)
birth_coverage = list("pessimistic"=birth_pess, "realistic"=birth_real, "optimistic"=birth_opt)
infant_ve = c("pessimistic"=VE_pessimistic_infant, "optimistic"=VE_optimistic_infant)

infant_interventions = data.frame(hosp=NA, date=NA,age=NA,season=NA, lhs=NA,coverage=NA, VE=NA)
for(c in 1:3){
    for(v in 1:2){
intervention_results = interventions_PI(birth_dose=birth_coverage[[c]], 
                                       cover_n=infant_coverage[[c]], 
                                       waningN=180,
                                       RRIn=1, 
                                       RRHn=infant_ve[v], 
                                       cover_s=c(rep(0,2387)), 
                                       waningS=36500, 
                                       RRIs=1, 
                                       RRHs=1) %>% 
  filter(age=="Infants <1y") %>% 
  mutate(coverage = names(infant_coverage[c]),
         VE = names(infant_ve[v]))
infant_interventions = rbind(infant_interventions, intervention_results)

}}


infant_interventions = data.frame(hosp=NA, date=NA,age=NA,season=NA, coverage=NA, VE=NA)
for(c in 1:3){
  for(v in 1:2){
    intervention_results = interventions(birth_dose=birth_coverage[[c]], 
                                            cover_n=infant_coverage[[c]], 
                                            waningN=180,
                                            RRIn=1, 
                                            RRHn=infant_ve[v], 
                                            cover_s=c(rep(0,2387)), 
                                            waningS=36500, 
                                            RRIs=1, 
                                            RRHs=1) %>% 
      filter(age=="Infants <1y") %>% 
      mutate(coverage = names(infant_coverage[c]),
             VE = names(infant_ve[v]))
    infant_interventions = rbind(infant_interventions, intervention_results)
    
  }}


saveRDS(infant_interventions,"infant_interventions_full.rds")

# Timing of birth dose administration -------------------------------------
#Oct - Mar
birth_dose_t1 = c(rep(0,2283),rep(.5,26),rep(0,26),rep(.75,26),rep(0,26))
#Nov - Mar
birth_dose_t2 = c(rep(0,2287),rep(.5,22),rep(0,30),rep(.75,22),rep(0,26))
#Oct - Feb
birth_dose_t3 = c(rep(0,2283),rep(.5,22),rep(0,30),rep(.75,22),rep(0,30))
#Nov - Feb
birth_dose_t4 = c(rep(0,2287),rep(.5,18),rep(0,34),rep(.75,18),rep(0,30))

birth_dose_timing = list("Oct_Mar"=birth_dose_t1,"Nov_Mar" = birth_dose_t2,
                      "Oct_Feb"=birth_dose_t3,"Nov_Feb"=birth_dose_t4)

timing_interventions = data.frame(hosp=NA, date=NA,age=NA,season=NA, lhs=NA,timing=NA)
for(t in 1:4){
    intervention_results = interventions_PI(birth_dose=birth_dose_timing[[t]], 
                                            cover_n=child_real, 
                                            waningN=180,
                                            RRIn=1, 
                                            RRHn=.2, 
                                            cover_s=c(rep(0,2387)), 
                                            waningS=36500, 
                                            RRIs=1, 
                                            RRHs=1) %>% 
      filter(age=="Infants <1y") %>% 
      mutate(timing = names(birth_dose_timing[t]))
    timing_interventions = rbind(timing_interventions, intervention_results)
    
  }


timing_interventions = data.frame(hosp=NA, date=NA,age=NA,season=NA,timing=NA)
for(t in 1:4){
  intervention_results = interventions(birth_dose=birth_dose_timing[[t]], 
                                          cover_n=child_real, 
                                          waningN=180,
                                          RRIn=1, 
                                          RRHn=.2, 
                                          cover_s=c(rep(0,2387)), 
                                          waningS=36500, 
                                          RRIs=1, 
                                          RRHs=1) %>% 
    filter(age=="Infants <1y") %>% 
    mutate(timing = names(birth_dose_timing[t]))
  timing_interventions = rbind(timing_interventions, intervention_results)
  
}

saveRDS(timing_interventions,"timing_interventions_full.rds")


# Get confidence intervals  -----------------------------------------------

baseline_seniors = no_interventions_full %>% 
  filter(age=="Seniors",lhs != "point estimate") %>% 
  group_by(season, coverage, VE, lhs) %>% 
  summarize(hosp = sum(hosp)) %>% 
  ungroup() %>% 
  group_by(season,coverage, VE) %>% 
  summarize(median = median(hosp),
            lower = quantile(hosp,probs=0.025),
            upper = quantile(hosp,probs=0.975))
  

seniors = readRDS("senior_interventions_full.rds") %>% 
  mutate(scenario = paste(coverage, VE, sep="_"),
         date = as.Date(date)) %>% 
  filter(lhs != "point estimate") %>% 
  group_by(season, coverage, VE, lhs) %>% 
  summarize(hosp = sum(hosp)) %>% 
  ungroup() %>% 
  group_by(season,coverage, VE) %>% 
  summarize(median = median(hosp),
            lower = quantile(hosp,probs=0.025),
            upper = quantile(hosp,probs=0.975)) #%>% 
  #ungroup()%>% 
 # group_by(season, scenario) %>% 
  #summarize(total = sum(hosp),
         #lower = sum(hosp),
         #upper = sum(hosp))

combined = rbind(baseline_seniors, seniors) %>% 
  mutate(coverage = factor(coverage, levels=c("None","pessimistic","optimistic"),
                           labels=c("None","Pessimistic","Optimistic")),
         VE = factor(VE, levels=c("None","pessimistic","optimistic"),
                           labels=c("None","Pessimistic","Optimistic")),
         season = factor(season,levels=c("season23_24","season24_25"),
                         labels=c("2023-24","2024-25")))

percent_diff = seniors %>% 
  left_join(baseline_seniors %>% dplyr::select(season, "baseline"=median), by="season") %>% 
  mutate(diff = round((1-median/baseline)*100),
         season = factor(season,levels=c("season23_24","season24_25"),
                         labels=c("2023-24","2024-25")),
         coverage=coverage.x,
         coverage = factor(coverage, levels=c("None","pessimistic","optimistic"),
                           labels=c("None","Pessimistic","Optimistic")),
         VE = factor(VE, levels=c("None","pessimistic","optimistic"),
                     labels=c("None","Pessimistic","Optimistic")),) %>% 
  select(-coverage.y, -coverage.x)

senior_results = ggplot()+
  theme_bw()+
  geom_bar(data=combined,aes(x=coverage, y=median, group=VE, fill=VE),stat="identity",position=position_dodge(preserve="single"))+
  facet_grid(rows=vars(season))+
  geom_errorbar(data=combined,aes(x=coverage, ymin=lower, ymax=upper, group=VE),width=0.1,position = position_dodge(0.9, preserve="single"))+
  geom_text(data=percent_diff, aes(x=coverage, y=median+50, group=VE, label=paste0("-",diff)),position = position_dodge(.9),hjust=0,size=5)+
  labs(x="Vaccine Coverage", y="RSV Hospitalizations")+
  scale_fill_manual(name="VE",values=c("maroon4","goldenrod","lightblue"))

ggsave("senior_results.png",height=5.2, width=6.2,units="in")


baseline_infants = no_interventions_full %>% 
  filter(age=="Infants <1y",lhs != "point estimate") %>% 
  group_by(season, coverage, VE, lhs) %>% 
  summarize(hosp = sum(hosp)) %>% 
  ungroup() %>% 
  group_by(season,coverage, VE) %>% 
  summarize(median = median(hosp),
            lower = quantile(hosp,probs=0.025),
            upper = quantile(hosp,probs=0.975))


infants = readRDS("infant_interventions_full.rds") %>% 
  mutate(scenario = paste(coverage, VE, sep="_"),
         date = as.Date(date)) %>% 
  filter(lhs != "point estimate") %>% 
  group_by(season, coverage, VE, lhs) %>% 
  summarize(hosp = sum(hosp)) %>% 
  ungroup() %>% 
  group_by(season,coverage, VE) %>% 
  summarize(median = median(hosp),
            lower = quantile(hosp,probs=0.025),
            upper = quantile(hosp,probs=0.975)) #%>% 
#ungroup()%>% 
# group_by(season, scenario) %>% 
#summarize(total = sum(hosp),
#lower = sum(hosp),
#upper = sum(hosp))

infants = infant_interventions %>% 
  mutate(scenario = paste(coverage, VE, sep="_"),
         date = as.Date(date)) %>% 
  filter(#lhs != "point estimate"
    !is.na(season)) %>% 
  group_by(season, coverage, VE) %>% 
  summarize(median = sum(hosp)) %>% 
  ungroup() %>% 
 # group_by(season,coverage, VE) %>% 
  mutate(lower = median*.85,
            upper = median*1.2) #%>% 
#ungroup()%>% 
# group_by(season, scenario) %>% 
#summarize(total = sum(hosp),
#lower = sum(hosp),
#upper = sum(hosp))

timing = timing_interventions %>% 
  mutate(#scenario = paste(coverage, VE, sep="_"),
         date = as.Date(date)) %>% 
  filter(#lhs != "point estimate"
    !is.na(season)) %>% 
  group_by(season, timing) %>% 
  summarize(median = sum(hosp)) %>% 
  ungroup() %>% 
  # group_by(season,coverage, VE) %>% 
  mutate(lower = median*.85,
         upper = median*1.2) #

combined2 = rbind(baseline_infants, infants) %>% 
  mutate(coverage = factor(coverage, levels=c("None","pessimistic","realistic","optimistic"),
                           labels=c("None","Pessimistic","Realistic","Optimistic")),
         VE = factor(VE, levels=c("None","pessimistic","optimistic"),
                     labels=c("None","Pessimistic","Optimistic")),
         season = factor(season,levels=c("season23_24","season24_25"),
                         labels=c("2023-24","2024-25")))

percent_diff = infants %>% 
  left_join(baseline_infants %>% dplyr::select(season, "baseline"=median), by="season") %>% 
  mutate(diff = round((1-median/baseline)*100),
         season = factor(season,levels=c("season23_24","season24_25"),
                         labels=c("2023-24","2024-25")),
         coverage=coverage.x,
         coverage = factor(coverage, levels=c("None","pessimistic","realistic","optimistic"),
                           labels=c("None","Pessimistic","Realistic","Optimistic")),
         VE = factor(VE, levels=c("None","pessimistic","optimistic"),
                     labels=c("None","Pessimistic","Optimistic")),) %>% 
  select(-coverage.y, -coverage.x)

infant_results = ggplot()+
  theme_bw()+
  geom_bar(data=combined2,aes(x=coverage, y=median, group=VE, fill=VE),stat="identity",position=position_dodge(preserve="single"))+
  facet_grid(rows=vars(season))+
  geom_errorbar(data=combined2,aes(x=coverage, ymin=lower, ymax=upper, group=VE),width=0.1,position = position_dodge(0.9, preserve="single"))+
  geom_text(data=percent_diff, aes(x=coverage, y=median+50, group=VE, label=paste0("-",diff)),position = position_dodge(.9),hjust=0,size=5)+
  labs(x="Vaccine Coverage", y="RSV Hospitalizations")+
  scale_fill_manual(name="VE",values=c("maroon4","goldenrod","lightblue"))
infant_results

ggsave("infant_results.png",height=5.2, width=7,units="in")

#timing scenarios 

timing = readRDS("timing_interventions_full.rds") %>% 
  filter(lhs!="point estimate") %>% 
  group_by(season, timing, lhs) %>% 
  summarize(hosp = sum(hosp)) %>% 
  ungroup() %>% 
  group_by(season,timing) %>% 
  summarize(median = median(hosp),
            lower = quantile(hosp,probs=0.025),
            upper = quantile(hosp,probs=0.975)) #%>% 



# Testing infant scenarios ------------------------------------------------
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
      RRHn*((RRHm+2)/3)*hosp[i]*RRIn*Mn[,i]*lambda1[,i]+ #note, this is adding nirsevimab impact on top of maternal immunity (might be too high?)
      #RRHn*hosp[i]*RRIn*Mn[,i]*lambda1[,i]+ #alternative approach is to add nirsevimab impact on it's own
      RRHn*((RRHm+3)/4)*hosp[i]*RRIn*N[,i]*lambda1[,i]+
      #RRHn*hosp[i]*RRIn*N[,i]*lambda1[,i]+# see note above 
      hosp[i]*S0[,i]*lambda1[,i]+
      hosp[i]*parmset$sigma1*S1[,i]*lambda1[,i]+
      hosp[i]*parmset$sigma2*S2[,i]*lambda1[,i]+
      hosp[i]*parmset$sigma3*S3[,i]*lambda1[,i]+
      RRHs*hosp[i]*parmset$RRIs*parmset$sigma3*Vs1[,i]*lambda1[,i]+
      RRHs*hosp[i]*parmset$RRIs*parmset$sigma3*Vs2[,i]*lambda1[,i]
  }
  
  H2=cbind(H1[,1],H1[,2],H1[,3],H1[,4],H1[,5],H1[,6],rowSums(H1[,7:8]),rowSums(H1[,9:12]),H1[,13])
  
  
  infants = data.frame(hosp = rowSums(H2[352:455,1:6]),
                       date = dates$date[2284:2387], 
                       age="Infants <1y")# intervention given to infants <8m but benefit is intended for all infants 
  
  results = infants %>%
    mutate(season = ifelse(date>='2024-10-01',"season24_25","season23_24")) %>% 
    group_by(season) %>% 
    summarize(hosp=sum(hosp))
  
  return(results)
}

no_interventions = interventions(birth_dose=c(rep(0,2387)), 
                                    cover_n=c(rep(0,2387)), 
                                    waningN=36500,
                                    RRIn=1, 
                                    RRHn=1, 
                                    cover_s=c(rep(0,2387)), 
                                    waningS=36500, 
                                    RRIs=1, 
                                    RRHs=1) %>% 
  mutate(coverage="None",
         VE="None")

child_vaccination = read_csv("https://raw.githubusercontent.com/midas-network/flu-scenario-modeling-hub_resources/main/Rd2_datasets/Age_Specific_Coverage_Flu_RD2_2022_23_Sc_A_B_C_D.csv") %>% 
  filter(Geography=="Washington", Age=="6 Months - 4 Years") %>% 
  select("date"=Week_Ending_Sat, "flu_coverage"=flu.coverage.rd2.sc_A_B_C_D) %>% 
  mutate(flu_coverage = ifelse(date=="2022-08-13",0,flu_coverage),
         flu_cov_cum = flu_coverage/100,
         flu_cov_week = flu_cov_cum - lag(flu_cov_cum,1),
         week_lag = lag(flu_cov_week,6),
         cum_lag = lag(flu_coverage,6)) %>% 
  replace(is.na(.),0)


#(25%/50%)
child_pess = c(rep(0,2274),child_vaccination$week_lag[1:22]*.5,rep(0,30),child_vaccination$flu_cov_week[1:22]*.78,rep(0,29),child_vaccination$flu_cov_week[1:9]*.78)
birth_pess = c(rep(0,2283),rep(.25,26),rep(0,26),rep(.5,26),rep(0,26))
#(50%/75%)
child_real = c(rep(0,2274),child_vaccination$week_lag[1:22]*.9,rep(0,30),child_vaccination$flu_cov_week[1:22]*1.17,rep(0,29),child_vaccination$flu_cov_week[1:9]*1.17)
birth_real = c(rep(0,2283),rep(.5,26),rep(0,26),rep(.75,26),rep(0,26))
#(75%/90%)
child_opt = c(rep(0,2274),child_vaccination$week_lag[1:22]*1.4,rep(0,30),child_vaccination$flu_cov_week[1:22]*1.3,rep(0,30),child_vaccination$flu_cov_week[1:9]*1.3)
birth_opt = c(rep(0,2283),rep(.75,26),rep(0,26),rep(.9,26),rep(0,26))

option1 = interventions(birth_dose=birth_real, 
                                        cover_n=child_real, 
                                        waningN=180,
                                        RRIn=1, 
                                        RRHn=.2, 
                                        cover_s=c(rep(0,2387)), 
                                        waningS=36500, 
                                        RRIs=1, 
                                        RRHs=1) 
