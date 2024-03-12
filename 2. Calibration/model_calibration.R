
rm(list=ls())
library(deSolve)
library(RColorBrewer)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(sjstats)
library(cdcfluview)
library(zoo)
library(imputeTS)
library(lhs)

"%notin%" = Negate('%in%')


# Step 1 - upload necessary data ------------------------------------------

#starting with demographic data 
#population is divided in to the M and S0 compartments 
#seed 1 infection in each compartment >6m, burn-in from 1995 to present 
#use the version you saved without the immunization compartments (Mn,Mv,N,Si,Vs1,Vs2)
yinit = readRDS("1. Data/Demographic Data/yinit.rds") 
yinit=as.matrix(yinit)
yinit.vector = readRDS("1. Data/Demographic Data/yinit.vector.rds")

#birth rates 
birth <-readRDS('1. Data/Demographic Data/births_kingcounty.rds') %>% select(-date)
birth = as.matrix(birth)

#vector of data from burn-in through projection period 
dates = data.frame(date=seq(from=as.Date("1995-01-07"), to=as.Date("2025-10-01"), by="weeks"))
tmax = 1317 #end of pre-pandemic period (2020-03-28) 
tmax2 = 1506 #end of rebound period (2023-11-11)
tmax3 = 1535#projection for the 2023-2024 season (2024-06-01)

#Upload the POLYMOD contact matrix, already aggregated to relevant age groups 
contact <- readRDS('1. Data/Demographic Data/contact_POLYMOD.rds')


#upload the other parameters - weekly seeding and death/migration rates 
seed = readRDS("1. Data/Demographic Data/other_parms.rds")[1] #seeding 1 infection per 100,000 population each week
introductions=c(rep(seed,tmax)) #vector through burn-in period 
um = readRDS("1. Data/Demographic Data/other_parms.rds")[2]*-1 #death death and migration rate 
npi = c(rep(1,tmax)) #vector through burn-in period - not NPIs during the -pre-pandemic period 
pop_2022 = readRDS("1. Data/Demographic Data/other_parms.rds")[3] #2022 population size (most recent available data)


#weekly time series of RSV hospitalizations (all ages)- before March 2020 
rsv_pre = readRDS('1. Data/RSV Data/rsv_ts.rds') %>% 
  select(date, rsvH=hrsv_smooth) %>% 
  filter(date<='2020-03-28') %>% 
  mutate(rsvH = round(rsvH))
plot(rsv_pre$rsvH)


burn_in_times <- seq(1, tmax, by = 1) # gives a sequence of weeks 

# list all fixed parameters  ------------------------------------------

parmset<-list(PerCapitaBirthsYear=birth,
              WidthAgeClassMonth=c(rep(2,times=6), 12,12*3, 60, 120, 240, 240, 240),#time spend in each age class (months)
              DurationMatImmunityDays = 112,#duration of maternal immunity (days)
              RRHm = 0.7,#relative risk of hospitalization given infection for those with maternal immunity
              recover1 = 365.25*.5, #days spent in R1 compartment
              recover2 = 365.25*.5, #days spent in R2 compartment
              recover3 = 365.25,#days spent in R3 compartment 
              recover4 = 365.25,#days spent in R4 compartment
              um=um, #net death and migration rate 
              rho1=0.75,#relative infectiousness following first infection 
              rho2=0.51,#relative infectiousness following 2+ infections 
              dur.days1=10, #duration of first infection (days)
              dur.days2=7, #duration of second infection (days)
              dur.days3=5, #duration of third+ infection (days)
              yinit.matrix=yinit, #initial compartments 
              q=1,
              c2=contact,
              sigma1=0.76,#reduced susceptibility following first infection 
              sigma2=0.6,#reduced susceptibility following second infection
              sigma3=0.4,#reduced susceptibility following third infection (and with maternal immunity)
              length.step = 7,
              time.step='week',
              npi=npi,
              seed=seed,#weekly seeding of infections 
              introductions=introductions)


# MSIRS Model - MLE calibration -------------------------------------------------------


## Read in transmission dynamic model
source("Calibration/model_dynamics.R")


fitmodel <-  function(parameters,dat) {  
  protrans <- parameters[1] # parameter for baseline transmission rate 
  baseline.txn.rate = 6+(3*(exp(protrans))) / (1+exp(protrans)) #transform between 6 and 9
  amp <- parameters[2] # parameter for seasonal amplitude
  b1 <- exp(amp) #ensure positive 
  trans <- parameters[3] # parameter for seasonal phase
  phi <-  (2*pi*(exp(trans))) / (1+exp(trans)) # transform to between 0 and 2pi 
  reporting_rate = 1/(1+exp(-parameters[4])) #reporting rate = how many infections are reported as hospitalizations in the data 

  # Simulate the model with initial conditions and timesteps defined above, and parameter values from function call
  results <- ode(y=yinit.vector, method = "ode45", t=burn_in_times,  
                 func=MSIRS_dynamics, 
                 parms=c(parmset,baseline.txn.rate=baseline.txn.rate,b1=b1,phi=phi))
  
  t0 <- nrow(rsv_pre)
  burnN <- tmax-t0
  St <- results[-c(1:burnN),-1]
  I1 <- St[,grep('I1', colnames(St))]
  I2 <- St[,grep('I2', colnames(St))]
  I3 <- St[,grep('I3', colnames(St))]
  I4 <- St[,grep('I4', colnames(St))]
  R1 <- St[,grep('R1', colnames(St))]
  R2 <- St[,grep('R2', colnames(St))]
  R3 <- St[,grep('R3', colnames(St))]
  R4 <- St[,grep('R4', colnames(St))]
  S1 <- St[,grep('S1', colnames(St))]
  S2 <- St[,grep('S2', colnames(St))]
  S3 <- St[,grep('S3', colnames(St))]
  S0 <- St[,grep('S0', colnames(St))]
  M <- St[,grep('M', colnames(St))]
  
  beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2# q depends on transmission type (whether depends on population density or not); c2 is the contact matrix
  
  lambda1=matrix(0,nrow=t0,ncol=13)#Force of infection
  for (t in 1:t0){
  lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+seed)%*%beta)/sum(St[t,]))}
  
   H1=matrix(0,nrow=t0,ncol=13)#Number of infections by age
  for (i in 1:13){
    H1[,i]=
      parmset$RRHm*parmset$sigma3*M[,i]*lambda1[,i]+ #assume infants with maternal immunity have protection from infection similar to adults 
      S0[,i]*lambda1[,i]+
      parmset$sigma1*S1[,i]*lambda1[,i]+
      parmset$sigma2*S2[,i]*lambda1[,i]+
      parmset$sigma3*S3[,i]*lambda1[,i]}
  

  H = rowSums(H1)*reporting_rate #combine into single time series 
  
  LL <- sum(dpois(x = dat, lambda =H, log = TRUE)) # fit to timeseries

  return(LL)
}


# Run optimization function  --------------------------------------

fitLL <- optim(par = c(0,-2,0,0),
               fn = fitmodel,        # the distance function to optimize
               dat = rsv_pre$rsvH,  # the dataset to fit to (dpois function)
               control = list(fnscale=-1))# the log likelihood is negative; here we maximize the log likelihood)
#save your parameters 
#saveRDS(fitLL, "Calibration/parameters_6Mar24.rds")
fitLL = readRDS("Calibration/parameters_6Mar24.rds")



baseline.txn.rate=6+(3*(exp(fitLL$par[1]))) / (1+exp(fitLL$par[1]))
b1=exp(fitLL$par[2])
phi=(2*pi*(exp(fitLL$par[3]))) / (1+exp(fitLL$par[3]))
reporting_rate = 1/(1+exp(-fitLL$par[4]))

# Plot results with fit parameters   ---------------------------------------------

output <- ode(y=yinit.vector, t=burn_in_times,method = "ode45",
              func=MSIRS_dynamics, 
              parms=c(parmset,
                      baseline.txn.rate=baseline.txn.rate,
                      b1=b1,
                      phi=phi))
                    

t0 <- nrow(rsv_pre)
al <- nrow(yinit)
output <- tail(output,t0)
St <- output[,-1]
I1 <- St[,grep('I1', colnames(St))]
I2 <- St[,grep('I2', colnames(St))]
I3 <- St[,grep('I3', colnames(St))]
I4 <- St[,grep('I4', colnames(St))]
S1 <- St[,grep('S1', colnames(St))]
S2 <- St[,grep('S2', colnames(St))]
S3 <- St[,grep('S3', colnames(St))]
S0 <- St[,grep('S0', colnames(St))]
R1 <- St[,grep('R1', colnames(St))]
R2 <- St[,grep('R2', colnames(St))]
R3 <- St[,grep('R3', colnames(St))]
R4 <- St[,grep('R4', colnames(St))]
M<- St[,grep('M', colnames(St))]

beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2


lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
for (t in 1:t0){
lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+parmset$seed)%*%beta)/sum(St[t,]))}

#Check infections and attack rates first 

weekly_infections=matrix(0,nrow=t0,ncol=al)#Number of incident infections by age
for (i in 1:al){
  weekly_infections[,i]=parmset$sigma3*M[,i]*lambda1[,i]+
    S0[,i]*lambda1[,i]+
    parmset$sigma1*S1[,i]*lambda1[,i]+
    parmset$sigma2*S2[,i]*lambda1[,i]+
    parmset$sigma3*S3[,i]*lambda1[,i]}

Infections = data.frame(weekly_infections)
names(Infections)=c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y","2-4Y","5-9Y","10-19Y","20-39Y","40-59Y","60Y+")
Infections$date = dates$date[(tmax-nrow(rsv_pre)+1):tmax]
Infections = Infections %>% 
  mutate(mmwr_week(date),
         season = case_when(mmwr_week<40 ~ paste0(mmwr_year-1,"-",mmwr_year),
                            mmwr_week>=40 ~ paste0(mmwr_year,"-",mmwr_year+1)))%>% 
  group_by(season) %>% 
  summarize_at(vars(`<2m`,`2-3m`,`4-5m`,`6-7m`,`8-9m`,`10-11m`,`1Y`,`2-4Y`,`5-9Y`,`10-19Y`,`20-39Y`,`40-59Y`,`60Y+`),sum)


pop=matrix(0,nrow=t0,ncol=al)#population - denominator for the attack rates 
for (i in 1:al){
  pop[,i]=M[,i]+
    S0[,i]+S1[,i]+S2[,i]+S3[,i]+
    I1[,i]+I2[,i]+I3[,i]+I4[,i]+
    R1[,i]+R2[,i]+R3[,i]+R4[,i]}

pop = data.frame(pop)
names(pop)=c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y","2-4Y","5-9Y","10-19Y","20-39Y","40-59Y","60Y+")
pop$date = dates$date[(tmax-nrow(rsv_pre)+1):tmax]
pop = pop %>% 
  mutate(mmwr_week(date),
         season = case_when(mmwr_week<40 ~ paste0(mmwr_year-1,"-",mmwr_year),
                            mmwr_week>=40 ~ paste0(mmwr_year,"-",mmwr_year+1)))%>% 
  group_by(season) %>% 
  summarize_at(vars(`<2m`,`2-3m`,`4-5m`,`6-7m`,`8-9m`,`10-11m`,`1Y`,`2-4Y`,`5-9Y`,`10-19Y`,`20-39Y`,`40-59Y`,`60Y+`),mean)

sum(pop[4,2:14])#check to make sure population size is approximately equal to the population size in 2020 (end if burn-in period)
#Infants should have highest attack rate ~between 40-70% 
#ideally children should have higher attack rates than adults 
attack_rates = Infections[2:4,2:14]/pop[2:4,2:14]*100
attack_rates

#get hospitalizations and compare with target data

#scale down all infections to hospitalizations 
H=rowSums(weekly_infections)*reporting_rate
plot(H)
H <- data.frame(H)
H$date <- dates$date[(tmax-nrow(rsv_pre)+1):tmax]


plot1 = ggplot()+
  theme_bw()+
  geom_area(data=rsv_pre, aes(x=date, y=rsvH),fill="seashell3")+
  geom_line(data=H, aes(x=date, y=H),color="navy",linewidth=1.5)+
  labs(x=NULL, y="RSV Hospitalizations")
plot1



# Fit NPI periods  --------------------------------------------------------
#time series of RSV hospitalizations for the full time period 
rsv_pand = readRDS('Data/RSV Data/rsv_ts.rds') %>% 
  select(date, rsvH=hrsv_smooth) %>% 
  filter(date<'2023-11-12') %>% 
  mutate(rsvH = round(rsvH))
plot(rsv_pand$rsvH)


rebound_times <- seq(1,tmax2 , by = 1) 

fit_rebound <-  function(parameters,dat) { 
  
  #parameters for duration of NPI (in terms of weeks, treat as whole numbers)
  time1 = round(exp(parameters[1]))
  time2 = round(exp(parameters[2]))
  time3 = round(exp(parameters[3]))
  time4 = round(exp(parameters[4]))
  time5 = round(exp(parameters[5]))
  
  #parameters for magnitude or NPI (proportion decrease in beta - must be <1)
  npi1 = 1/(1+exp(-parameters[6]))
  npi2 = 1/(1+exp(-parameters[7]))
  npi3 = 1/(1+exp(-parameters[8]))
  npi4 = 1/(1+exp(-parameters[9]))
  
  #fit a new reporting rate during the pandemic 
  reporting_rate2 = reporting_rate+(reporting_rate*(exp(parameters[10]))) / (1+exp(parameters[10])) 
 
  
introductions = data.frame(intros=c(rep(seed,tmax),rep(0,44),rep(NA,24),rep(seed,179))) %>% 
  mutate(intros = na_interpolation(intros, method="linear"))
introductions = introductions$intros

npi = data.frame(npis=c(rep(1,tmax),rep(npi1,13),rep(npi2,time1),rep(NA,time2),rep(npi3,time3),rep(npi4,time4),rep(NA,time5),rep(1,160)))%>% 
  mutate(npis= na_interpolation(npis, method="linear"))
npi = npi$npis

parmset_rebound<-list(
              baseline.txn.rate=baseline.txn.rate,
              phi=phi,
              b1=b1,
              PerCapitaBirthsYear=birth,
              WidthAgeClassMonth=c(rep(2,times=6), 12,12*3, 60, 120, 240, 240, 240),
              DurationMatImmunityDays = 112,
              RRHm = 0.7,
              recover1 = 365.25*.5,
              recover2 = 365.25*.5,
              recover3 = 365.25,
              recover4 = 365.25,
              um=um,
              rho1=0.75,
              rho2=0.51,
              dur.days1=10,
              dur.days2=7,
              dur.days3=5,
              yinit.matrix=yinit,
              q=1,
              c2=contact,
              sigma1=0.76,
              sigma2=0.6,
              sigma3=0.4,
              length.step = 7,
              time.step='week',
              npi=npi,
              seed=seed,
              introductions=introductions)


output <- ode(y=yinit.vector, t=rebound_times,method = "ode45",
              func=MSIRS_dynamics, 
              parms=c(parmset_rebound))


t0 <- nrow(rsv_pand)
al <- nrow(yinit)
output <- tail(output,t0)
St <- output[,-1]
I1 <- St[,grep('I1', colnames(St))]
I2 <- St[,grep('I2', colnames(St))]
I3 <- St[,grep('I3', colnames(St))]
I4 <- St[,grep('I4', colnames(St))]
S1 <- St[,grep('S1', colnames(St))]
S2 <- St[,grep('S2', colnames(St))]
S3 <- St[,grep('S3', colnames(St))]
S0 <- St[,grep('S0', colnames(St))]
M<- St[,grep('M', colnames(St))]

contact2 = npi[(tmax2-t0+1):tmax2]
intro2 = introductions[(tmax2-t0+1):tmax2]

lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
for (t in 1:t0){
  beta <-  baseline.txn.rate/(parmset_rebound$dur.days1/7)/(sum(yinit)^(1-parmset_rebound$q))*parmset_rebound$c2*contact2[t]
  lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset_rebound$rho1*I2[t,]+parmset_rebound$rho2*I3[t,]+parmset_rebound$rho2*I4[t,]+intro2[t])%*%beta)/sum(St[t,]))}



H1=matrix(0,nrow=t0,ncol=al)#Number of hospitalizations by age
for (i in 1:al){
  H1[,i]=
    parmset_rebound$RRHm*parmset_rebound$sigma3*M[,i]*lambda1[,i]+
    S0[,i]*lambda1[,i]+
    parmset_rebound$sigma1*S1[,i]*lambda1[,i]+
    parmset_rebound$sigma2*S2[,i]*lambda1[,i]+
    parmset_rebound$sigma3*S3[,i]*lambda1[,i]}

H=rowSums(H1)
H2 = H[1:(nrow(rsv_pre)+10)]*reporting_rate #use the pre-pandemic reporting rate 
H3 = H[(nrow(rsv_pre)+11):nrow(rsv_pand)]*reporting_rate2 #use the new reporting_rate
H4 = c(H2,H3)

LL <- sum(dpois(x = dat, lambda =H4, log = TRUE)) # fit to timeseries

}

fitLL <- optim(par = c(3.9, 3.4, 2.3, 2.8, 3.2, .41, .41, 1.1, .62,0.5),
               fn = fit_rebound,        # the distance function to optimize
               dat = rsv_pand$rsvH,  # the dataset to fit to (dpois function)
               control = list(fnscale=-1)) # the log likelihood is negative; here we maximize the log likelihood

#saveRDS(fitLL,"Calibration/NPI_6Mar24.RDS")            
fitLL = readRDS("Calibration/NPI_6Mar24.RDS")

time1 = round(exp(fitLL$par[1]))
time2 = round(exp(fitLL$par[2]))
time3 = round(exp(fitLL$par[3]))
time4 = round(exp(fitLL$par[4]))
time5 = round(exp(fitLL$par[5]))

npi1 = 1/(1+exp(-fitLL$par[6]))
npi2 = 1/(1+exp(-fitLL$par[7]))
npi3 = 1/(1+exp(-fitLL$par[8]))
npi4 = 1/(1+exp(-fitLL$par[9]))

reporting_rate2 = reporting_rate+(reporting_rate*(exp(fitLL$par[10]))) / (1+exp(fitLL$par[10]))


introductions = data.frame(intros=c(rep(seed,tmax),rep(0,44),rep(NA,24),rep(seed,179))) %>% 
  mutate(intros = na_interpolation(intros, method="linear"))
introductions = introductions$intros

npi = data.frame(npis=c(rep(1,tmax),rep(npi1,13),rep(npi2,time1),rep(NA,time2),rep(npi3,time3),rep(npi4,time4),rep(NA,time5),rep(1,160)))%>% 
  mutate(npis= na_interpolation(npis, method="linear"))
npi = npi$npis


parmset_rebound<-list(
  baseline.txn.rate=baseline.txn.rate,
  phi=phi,
  b1=b1,
  PerCapitaBirthsYear=birth,
  WidthAgeClassMonth=c(rep(2,times=6), 12,12*3, 60, 120, 240, 240, 240),
  DurationMatImmunityDays = 112,
  RRHm = 0.7,
  recover1 = 365.25*.5,
  recover2 = 365.25*.5,
  recover3 = 365.25,
  recover4 = 365.25,
  um=um,
  rho1=0.75,
  rho2=0.51,
  dur.days1=10,
  dur.days2=7,
  dur.days3=5,
  yinit.matrix=yinit,
  q=1,
  c2=contact,
  sigma1=0.76,
  sigma2=0.6,
  sigma3=0.4,
  length.step = 7,
  time.step='week',
  npi=npi,
  seed=seed,
  introductions=introductions)



projection_times = seq(from=1,to=tmax3,by=1)
output <- ode(y=yinit.vector, t=projection_times,method = "ode45",
              func=MSIRS_dynamics, 
              parms=c(parmset_rebound))

projection_length = tmax3-tmax2
t0 <- nrow(rsv_pand)+projection_length
al <- nrow(yinit)
output <- tail(output,t0)
St <- output[,-1]
I1 <- St[,grep('I1', colnames(St))]
I2 <- St[,grep('I2', colnames(St))]
I3 <- St[,grep('I3', colnames(St))]
I4 <- St[,grep('I4', colnames(St))]

S1 <- St[,grep('S1', colnames(St))]
S2 <- St[,grep('S2', colnames(St))]
S3 <- St[,grep('S3', colnames(St))]
S0 <- St[,grep('S0', colnames(St))]
M<- St[,grep('M', colnames(St))]

contact2 = npi[(tmax3-t0+1):tmax3]
intro2 = introductions[(tmax3-t0+1):tmax3]

lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
for (t in 1:t0){
  beta <-  baseline.txn.rate/(parmset_rebound$dur.days1/7)/(sum(yinit)^(1-parmset_rebound$q))*parmset_rebound$c2*contact2[t]
lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset_rebound$rho1*I2[t,]+parmset_rebound$rho2*I3[t,]+parmset_rebound$rho2*I4[t,]+intro2[t])%*%beta)/sum(St[t,]))}


H1=matrix(0,nrow=t0,ncol=al)#Number of hospitalizations by age
for (i in 1:al){
  H1[,i]=
    parmset_rebound$RRHm*parmset_rebound$sigma3*M[,i]*lambda1[,i]+
    S0[,i]*lambda1[,i]+
    parmset_rebound$sigma1*S1[,i]*lambda1[,i]+
    parmset_rebound$sigma2*S2[,i]*lambda1[,i]+
    parmset_rebound$sigma3*S3[,i]*lambda1[,i]}

H=rowSums(H1)
H2 = H[1:(nrow(rsv_pre)+10)]*reporting_rate #use the pre-pandemic reporting rate 
H3 = H[(nrow(rsv_pre)+11):t0]*reporting_rate2 #use the new reporting_rate
H4 = c(H2,H3)
plot(H4)
H4=data.frame(H4)
H4$date = dates$date[(tmax3-t0+1):tmax3]

#all available data
rsv_all = readRDS('Data/RSV Data/rsv_ts.rds') %>% 
  select(date, rsvH=hrsv_smooth) %>% 
  mutate(rsvH = round(rsvH))

plot_rebound = ggplot()+
  theme_bw()+
  geom_area(data=rsv_all, aes(x=date, y=rsvH),fill="seashell3")+
  geom_line(data=H4, aes(x=date, y=H4),color="navy",linewidth=1.5)+
  theme(legend.position="bottom")+
  labs(x=NULL, y="RSV Hospitalizations")
plot_rebound

#apply agr distributions to the all ages curve 
agedist = readRDS("Data/RSV Data/age_distributions.rds") %>% 
  mutate(age=factor(age,levels=c("<6m",">6m","1-4yrs","5-59yrs","60+yrs"))) %>% 
  arrange(age)
agedist1= agedist %>% filter(period=="pre-pandemic") %>% select(prop)
agedist2= agedist %>% filter(period=="post-pandemic") %>% select(prop)

#we don't know what the age distribution will be for the projection period 
# we assume it will be between the pre-pandemic period and the rebound period 

#first apply age distributions for the pre-pandemic period to the rebound period
h1 = H3*as.numeric(agedist1[1,])
h2 = H3*as.numeric(agedist1[2,])
h3 = H3*as.numeric(agedist1[3,])
h4 = H3*as.numeric(agedist1[4,])
h5 = H3*as.numeric(agedist1[5,])

#apply age distributions for rebound period 
h1.2 = H3*as.numeric(agedist2[1,])
h2.2 = H3*as.numeric(agedist2[2,])
h3.2 = H3*as.numeric(agedist2[3,])
h4.2 = H3*as.numeric(agedist2[4,])
h5.2 = H3*as.numeric(agedist2[5,])

#age specific reporting rates 
Infections= cbind(rowSums(H1[,1:3]),rowSums(H1[,4:6]),rowSums(H1[,7:8]),rowSums(H1[,9:12]),H1[,13])
Infections=Infections[(nrow(rsv_pre)+11):t0,]

prop1 = sum(h1)/sum(Infections[,1])
prop2 = sum(h2)/sum(Infections[,2])
prop3 = sum(h3)/sum(Infections[,3])
prop4 = sum(h4)/sum(Infections[,4])
prop5 = sum(h5)/sum(Infections[,5])
age_rates1 = c(prop1,prop2,prop3,prop4,prop5)
age_rates1

prop1.2 = sum(h1.2)/sum(Infections[,1])
prop2.2 = sum(h2.2)/sum(Infections[,2])
prop3.2 = sum(h3.2)/sum(Infections[,3])
prop4.2 = sum(h4.2)/sum(Infections[,4])
prop5.2 = sum(h5.2)/sum(Infections[,5])
age_rates2 = c(prop1.2,prop2.2,prop3.2,prop4.2,prop5.2)
age_rates2


#when calculating confidence intervals we will sample between the range of these estimates 
#for point estimates we will use the mean between the two estimates 
age_rates_mean = (age_rates1+age_rates2)/2
age_rates = rbind(age_rates1, age_rates2, age_rates_mean)
saveRDS(age_rates, "age_specific_reporting_rates.rds")

# Latin Hypercube Sampling of fitted parameters ---------------------------
#LHS bounds - adding and subtracting 3% from the phase and amplitude 
b1_lower = b1*.97
b1_upper = b1*1.03
phi_lower = phi*.97
phi_upper = phi*1.03
h1_lower = age_rates1[1]
h1_upper = age_rates2[1]
h2_lower = age_rates1[2]
h2_upper = age_rates2[2]
h3_lower = age_rates1[3]
h3_upper = age_rates2[3]
h4_lower = age_rates1[4]
h4_upper = age_rates2[4]
h5_lower = age_rates1[5]
h5_upper = age_rates2[5]


#set-up LHS
#draw 100 samples 
set.seed(123)
h=100
lhs<-maximinLHS(h,7)


new_parms <- cbind(
  b1 = lhs[,1]*(b1_upper-b1_lower)+b1_lower,
  phi = lhs[,2]*(phi_upper-phi_lower)+phi_lower,
  h1 = lhs[,3]*(h1_upper-h1_lower)+h1_lower,
  h2 = lhs[,4]*(h2_upper-h2_lower)+h2_lower,
  h3 = lhs[,5]*(h3_upper-h3_lower)+h3_lower,
  h4 = lhs[,6]*(h4_upper-h4_lower)+h4_lower,
  h5 = lhs[,7]*(h5_upper-h5_lower)+h5_lower)
saveRDS(new_parms,"lhs_resampling100.rds")



#Draw 1000 samples 
set.seed(123)
h=1000
lhs<-maximinLHS(h,7)


new_parms2 <- cbind(
  b1 = lhs[,1]*(b1_upper-b1_lower)+b1_lower,
  phi = lhs[,2]*(phi_upper-phi_lower)+phi_lower,
  h1 = lhs[,3]*(h1_upper-h1_lower)+h1_lower,
  h2 = lhs[,4]*(h2_upper-h2_lower)+h2_lower,
  h3 = lhs[,5]*(h3_upper-h3_lower)+h3_lower,
  h4 = lhs[,6]*(h4_upper-h4_lower)+h4_lower,
  h5 = lhs[,7]*(h5_upper-h5_lower)+h5_lower)
saveRDS(new_parms2,"lhs_resampling1000.rds")


