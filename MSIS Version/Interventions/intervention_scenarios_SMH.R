
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
dates = data.frame(date=seq(from=as.Date("1980-01-05"), to=as.Date("2024-06-01"), by="weeks"))

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
tmax=nrow(dates)
wa_times <- seq(1, tmax, by =1) 

#upload parameters fit by MLE 
fitLL = readRDS("parameters_27Oct2023.rds")

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
source("Interventions/intervention_models_SMH.R")


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
  
  
  t0 <- 385 #just save results after burn-in period
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
  
  contact2 = npi[1933:2318]
  intro2 = introductions[1933:2318]
  
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
  
  H2 = cbind(rowSums(H1[,1:6]),rowSums(H1[,7:8]),rowSums(H1[,9:12]),H1[,13])
  
  results = as.data.frame(H2[356:385,1:4]) %>% mutate(date=dates$date[2289:2318])
  names(results) = c("<1y","1-4y","5-64y","65+y","date")
  
  
  
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
    
    
    t0 <- 385
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
    
    
    contact2 = npi[1933:2318]
    intro2 = introductions[1933:2318]
    
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
    
    
    H2 = cbind(rowSums(H1[,1:6]),rowSums(H1[,7:8]),rowSums(H1[,9:12]),H1[,13])
    
    projection_intervals = as.data.frame(H2[356:385,1:4]) %>% mutate(date=dates$date[2289:2318],
                                                                     lhs="point estimate")
    names(projection_intervals) = c("<1y","1-4y","5-64y","65+y","date","lhs")
   
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
  h <- 100        
  set.seed(123)
  lhs<-maximinLHS(h,5)
  
  new_parms <- cbind(
    baseline.txn.rate = lhs[,1]*(beta_upper-beta_lower)+beta_lower,
    b1 = lhs[,2]*(b1_upper-b1_lower)+b1_lower,
    phi = lhs[,3]*(phi_upper-phi_lower)+phi_lower,
    DurationMatImmunityDays = lhs[,4]*(duration_upper-duration_lower)+duration_lower,
    I_ex = lhs[,5]*(I_ex_upper-I_ex_lower)+I_ex_lower)
  
  #parameters for fitting the pandemic and rebound periods 
  for(l in 1:100){
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


t0 <- 385
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


contact2 = npi[1933:2318]
intro2 = introductions[1933:2318]

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

H2 = cbind(rowSums(H1[,1:6]),rowSums(H1[,7:8]),rowSums(H1[,9:12]),H1[,13])

estimates = as.data.frame(H2[356:385,1:4]) %>% mutate(date=dates$date[2289:2318],
                                                      lhs=l)
names(estimates) = c("<1y","1-4y","5-64y","65+y","date","lhs")

projection_intervals = rbind(projection_intervals, estimates)

  } 
  return(projection_intervals)
  }


#timing of senior vaccines 
#76% cumulative coverage
senior_vaccination = read_csv("https://raw.githubusercontent.com/midas-network/flu-scenario-modeling-hub_resources/main/Rd2_datasets/Age_Specific_Coverage_Flu_RD2_2022_23_Sc_A_B_C_D.csv") %>% 
  filter(Geography=="Maryland", Age=="65+ Years") %>% 
  select("date"=Week_Ending_Sat, "flu_coverage"=flu.coverage.rd2.sc_A_B_C_D) %>% 
  mutate(flu_cov_cum = flu_coverage/100,
         flu_cov_week = flu_cov_cum - lag(flu_cov_cum,1),
         flu_lag_cum = lag(flu_cov_cum,2), #start coverage first week of October 
         flu_lag_week = rollmean(lag(flu_cov_week,2),k=3,align="center",fill="extend")) %>% 
  replace(is.na(.),0)

#vaccination starts in August 
senior_opt = c(rep(0,2274),senior_vaccination$flu_lag_week*.4)
plot(cumsum(senior_opt[2274:2318]))
senior_pess = c(rep(0,2274),senior_vaccination$flu_lag_week*.2)
plot(cumsum(senior_pess[2274:2318]))

#timing of vaccination for catch-up vaccines (based on flu in 6mo-4yrs)
#75% cumulative coverage 
child_vaccination = read_csv("https://raw.githubusercontent.com/midas-network/flu-scenario-modeling-hub_resources/main/Rd2_datasets/Age_Specific_Coverage_Flu_RD2_2022_23_Sc_A_B_C_D.csv") %>% 
  filter(Geography=="Maryland", Age=="6 Months - 4 Years") %>% 
  select("date"=Week_Ending_Sat, "flu_coverage"=flu.coverage.rd2.sc_A_B_C_D) %>% 
  mutate(flu_cov_cum = flu_coverage/100,
         flu_cov_week = flu_cov_cum - lag(flu_cov_cum,1),
         flu_lag_cum = lag(flu_cov_cum,6), #start coverage first week of October 
         flu_lag_week = rollmean(lag(flu_cov_week,6),k=3,align="center",fill="extend")) %>% 
  replace(is.na(.),0)

infant_opt = c(rep(0,2274),child_vaccination$flu_lag_week[1:34]*.6,rep(0,10))
plot(cumsum(infant_opt[2274:2318]))
infant_pess = c(rep(0,2274),child_vaccination$flu_lag_week[1:34]*.2,rep(0,10))
plot(cumsum(infant_pess[2274:2318]))


cov_data = data.frame(cov = c(cumsum(senior_opt[2274:2318]),cumsum(senior_pess[2274:2318]),
                              cumsum(infant_opt[2274:2318]),cumsum(infant_pess[2274:2318])),
                      age = c(rep("Seniors",90),rep("Infants",90)),
                      scenario = c(rep("Optimistic",45),rep("Pessimistic",45),
                                   rep("Optimistic",45),rep("Pessimistic",45)),
                      date = dates$date[2274:2318])


plot1 = ggplot(data=cov_data)+
  geom_line(aes(x=date, y=cov*100,color=scenario),linewidth=1)+
  facet_grid(cols=vars(age))+
  labs(x=NULL, y="% Covered")
plot1
ggsave(plot=plot1, "coverage plots KC.png",height=6, width=13.5,units="in")

# Scenario A - Optimistic x Optimistic ------------------------------------
#Infant
#coverage = 60% of flu
#VE = 80% 
#Seniors
#coverage = 40% of flu
#VE = 90% 
scenario_A = interventions(birth_dose=infant_opt, 
                           cover_n=infant_opt, 
                           waningN=180,
                           RRIn=1, 
                           RRHn=.2, 
                           cover_s=senior_opt, 
                           waningS=730.5, 
                           RRIs=1, 
                           RRHs=.1)


scenarioPI_A = interventions_PI(birth_dose=infant_opt, 
                                cover_n=infant_opt, 
                                waningN=180,
                                RRIn=1, 
                                RRHn=.2, 
                                cover_s=senior_opt, 
                                waningS=730.5, 
                                RRIs=1, 
                                RRHs=.1) %>% 
  mutate(scenario = "A")

saveRDS(scenarioPI_A, "Scenario_A.rds")

SA_ts = scenarioPI_A %>% 
  pivot_longer(cols=`<1y`:`65+y`,names_to="age",values_to="value") %>% 
  group_by(date, age) %>% 
  summarize(median = median(value),
            lower = quantile(value, probs = 0.025),
            upper = quantile(value, probs = 0.975)) %>% 
  mutate(scenario = "A")

SA_cum = SA_ts %>% 
  ungroup() %>% 
  group_by(age) %>% 
  summarize(median = sum(median),
            lower = sum(lower),
            upper = sum(upper))%>% 
  mutate(scenario = "A")
sum(SA_cum$median)
sum(SA_cum$lower)
sum(SA_cum$upper)

 
# Scenario B - Optimistic x Pessimistic -----------------------------------
#Infant
#coverage = 60% of flu
#VE = 80% 
#Seniors
#coverage = 20% of flu
#VE = 70% 

scenario_B = interventions(birth_dose=infant_opt, 
                           cover_n=infant_opt, 
                           waningN=180,
                           RRIn=1, 
                           RRHn=.2, 
                           cover_s=senior_pess, 
                           waningS=730.5, 
                           RRIs=1, 
                           RRHs=.3)

scenario_B = interventions_PI(birth_dose=infant_opt, 
                           cover_n=infant_opt, 
                           waningN=180,
                           RRIn=1, 
                           RRHn=.2, 
                           cover_s=senior_pess, 
                           waningS=730.5, 
                           RRIs=1, 
                           RRHs=.3) %>% 
  mutate(scenario = "B")

saveRDS(scenario_B, "Scenario_B.rds")

SB_ts = scenario_B %>% 
  pivot_longer(cols=`<1y`:`65+y`,names_to="age",values_to="value") %>% 
  group_by(date, age) %>% 
  summarize(median = median(value),
            lower = quantile(value, probs = 0.025),
            upper = quantile(value, probs = 0.975)) %>% 
  mutate(scenario = "B")

SB_cum = SB_ts %>% 
  ungroup() %>% 
  group_by(age) %>% 
  summarize(median = sum(median),
            lower = sum(lower),
            upper = sum(upper))%>% 
  mutate(scenario = "B")
sum(SB_cum$median)
sum(SB_cum$lower)
sum(SB_cum$upper)


# Scenario C - Pessimistic x Optimistic -----------------------------------
#Infant
#coverage = 20% of flu
#VE = 60% 
#Seniors
#coverage = 40% of flu
#VE = 90% 

scenario_C = interventions(birth_dose=infant_pess, 
                           cover_n=infant_pess, 
                           waningN=180,
                           RRIn=1, 
                           RRHn=.4, 
                           cover_s=senior_opt, 
                           waningS=730.5, 
                           RRIs=1, 
                           RRHs=.1)

scenario_C = interventions_PI(birth_dose=infant_pess, 
                           cover_n=infant_pess, 
                           waningN=180,
                           RRIn=1, 
                           RRHn=.4, 
                           cover_s=senior_opt, 
                           waningS=730.5, 
                           RRIs=1, 
                           RRHs=.1) %>% 
  mutate(scenario = "C")

saveRDS(scenario_C, "Scenario_C.rds")

SC_ts = scenario_C %>% 
  pivot_longer(cols=`<1y`:`65+y`,names_to="age",values_to="value") %>% 
  group_by(date, age) %>% 
  summarize(median = median(value),
            lower = quantile(value, probs = 0.025),
            upper = quantile(value, probs = 0.975)) %>% 
  mutate(scenario = "C")

SC_cum = SC_ts %>% 
  ungroup() %>% 
  group_by(age) %>% 
  summarize(median = sum(median),
            lower = sum(lower),
            upper = sum(upper))%>% 
  mutate(scenario = "C")
sum(SC_cum$median)
sum(SC_cum$lower)
sum(SC_cum$upper)



# Scenario D - Pessimistic x Pessimistic  ---------------------------------
#Infant
#coverage = 20% of flu
#VE = 60% 
#Seniors
#coverage = 20% of flu
#VE = 70% 

scenario_D = interventions(birth_dose=infant_pess, 
                           cover_n=infant_pess, 
                           waningN=180,
                           RRIn=1, 
                           RRHn=.4, 
                           cover_s=senior_pess, 
                           waningS=730.5, 
                           RRIs=1, 
                           RRHs=.3) %>% 
  mutate(scenario=="D")

scenario_D = interventions_PI(birth_dose=infant_pess, 
                           cover_n=infant_pess, 
                           waningN=180,
                           RRIn=1, 
                           RRHn=.4, 
                           cover_s=senior_pess, 
                           waningS=730.5, 
                           RRIs=1, 
                           RRHs=.3) %>% 
  mutate(scenario="D")




saveRDS(scenario_D, "Scenario_D.rds")

SD_ts = scenario_D %>% 
  pivot_longer(cols=`<1y`:`65+y`,names_to="age",values_to="value") %>% 
  group_by(date, age) %>% 
  summarize(median = median(value),
            lower = quantile(value, probs = 0.025),
            upper = quantile(value, probs = 0.975)) %>% 
  mutate(scenario = "D")

SD_cum = SD_ts %>% 
  ungroup() %>% 
  group_by(age) %>% 
  summarize(median = sum(median),
            lower = sum(lower),
            upper = sum(upper))%>% 
  mutate(scenario = "D")
sum(SD_cum$median)
sum(SD_cum$lower)
sum(SD_cum$upper)





# Scenario E - Counter factual (no interventions) --------------------------
scenario_E = interventions(birth_dose=c(rep(0,2318)), 
                           cover_n=c(rep(0,2318)), 
                           waningN=36500,
                           RRIn=1, 
                           RRHn=1, 
                           cover_s=c(rep(0,2318)), 
                           waningS=36500, 
                           RRIs=1, 
                           RRHs=1)

scenario_E = interventions_PI(birth_dose=c(rep(0,2318)), 
                           cover_n=c(rep(0,2318)), 
                           waningN=36500,
                           RRIn=1, 
                           RRHn=1, 
                           cover_s=c(rep(0,2318)), 
                           waningS=36500, 
                           RRIs=1, 
                           RRHs=1) %>% 
  mutate(scenario = "E")

saveRDS(scenario_E, "Scenario_E.rds")

SE_ts = scenario_E %>% 
  pivot_longer(cols=`<1y`:`65+y`,names_to="age",values_to="value") %>% 
  group_by(date, age) %>% 
  summarize(median = median(value),
            lower = quantile(value, probs = 0.025),
            upper = quantile(value, probs = 0.975)) %>% 
  mutate(scenario = "E")

SE_cum = SE_ts %>% 
  ungroup(date, age) %>% 
  group_by(date) %>% 
  summarize(median = sum(median),
            lower = sum(lower),
            upper = sum(upper))%>% 
  mutate(scenario = "E")
sum(SE_cum$median)
sum(SE_cum$lower)
sum(SE_cum$upper)


all_ts = rbind(SA_ts, SB_ts, SC_ts, SD_ts, SE_ts)

plot2 = ggplot(data=all_ts)+
  theme_bw()+
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper,fill=scenario),alpha=.2)+
  geom_line(aes(x=date, y=median,color=scenario))+
  facet_grid(rows=vars(age),cols=vars(scenario))+
  labs(x=NULL, y="RSV Hospitalizations")
plot2
ggsave(plot=plot2,"scenario time series.png",height=6, width=13, units="in")


overall = all_ts %>% 
  group_by(date,scenario) %>% 
  summarize(median = sum(median),
            lower = sum(lower),
            upper = sum(upper))
plot3 = ggplot(data=overall)+
  theme_bw()+
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper,fill=scenario),alpha=.2)+
  geom_line(aes(x=date, y=median,color=scenario))+
  facet_grid(cols=vars(scenario))+
  labs(x=NULL, y="RSV Hospitalizations")
plot3
ggsave(plot=plot3,"scenario time series - overall.png",height=6, width=13, units="in")
