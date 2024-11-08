rm(list=ls())
library(deSolve)
library(tidyverse)
library(zoo)
library(imputeTS)
library(lhs)
library(tidycensus)
library(cdcfluview)
library(cowplot)
library(pracma)
library(lhs)
library(bbmle)

source("R/MSIRS_immunization_dynamics.R")

data = readRDS("DATA/fixed_parameters.rds")
parmset = data[[1]]
yinit=data[[2]]
yinit.vector=data[[3]]


timeseries = readRDS("DATA/time_series_public.rds") #data has been adjusted to account for testing changes during the pandemic 

time_series_full = timeseries %>% filter(date<'2023-10-01') %>% 
  mutate(hosp = round(hosp_rate/100000*2267000))
time_series_full = c(time_series_full$hosp)
plot(time_series_full)

time_series_pre = timeseries %>% filter(date<'2020-04-01') %>% 
  mutate(hosp = round(hosp_rate/100000*2267000))

time_series_pre = c(time_series_pre$hosp)
plot(time_series_pre)

age_dist = readRDS("DATA/age_distribution_public.rds")
age_dist = c(age_dist$proportion)
age_dist

#no immunizations for calibration 
parmset$monoclonal_birth =  rep(0, length(time_series_full) + 104+52)
parmset$monoclonal_catchup_01 =  rep(0, length(time_series_full) + 104+52)
parmset$monoclonal_catchup_23 =  rep(0, length(time_series_full) + 104+52)
parmset$monoclonal_catchup_45 =  rep(0, length(time_series_full) + 104+52)
parmset$monoclonal_catchup_67 =  rep(0, length(time_series_full) + 104+52)
parmset$maternal_vax <- rep(0, length(time_series_full) + 104+52)
parmset$senior_vax_75 <- rep(0, length(time_series_full) + 104+52)
parmset$senior_vax_60_74 <- rep(0, length(time_series_full) + 104+52)
parmset$introductions <- rep(parmset$seed, length(time_series_full) + 104+52)
parmset$npi <- rep(1, length(time_series_full) + 104+52)


# Fit to pre-pandemic years  ----------------------------------------------
fit_times <- seq(1, length(time_series_pre) + 104, by = 1)

fitmodel <-  function(parameters,dat) {
  protrans <- parameters[1] # parameter for baseline transmission rate
  baseline.txn.rate = exp(protrans) 
  amp <- parameters[2] # parameter for seasonal amplitude
  b1 <-  exp(amp)
  trans <- parameters[3] # parameter for seasonal phase
  phi <-  (2*pi*(exp(trans))) / (1+exp(trans)) # transform to between 0 and 2pi
  #Age-specific reporting fractions 
  report_infants <- 1 / (1 + exp(-parameters[4])) 
  report_children <- 1 / (1 + exp(-parameters[5]))
  report_adults <- 1 / (1 + exp(-parameters[6]))
  report_seniors60 <- 1 / (1 + exp(-parameters[7]))
  report_seniors75 <- 1 / (1 + exp(-parameters[8]))

  
  results <- ode(y=yinit.vector, method = "ode45", times=fit_times,
                 func=MSIRS_immunization_dynamics,
                 parms=c(parmset,baseline.txn.rate=baseline.txn.rate,b1=b1,phi=phi))
                 
  
  t0 <- length(time_series_pre)
  al <- nrow(yinit)
  output <- tail(results,t0)
  St <- output[,-1]
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
  M0 <- St[,grep('M0', colnames(St))]
  Si<- St[,grep('Si', colnames(St))]
  Mn1<- St[,grep('Mn1', colnames(St))]
  Mn2<- St[,grep('Mn2', colnames(St))]
  Mv1<- St[,grep('Mv1', colnames(St))]
  Mv2<- St[,grep('Mv2', colnames(St))]
  N1<- St[,grep('N1', colnames(St))]
  N2<- St[,grep('N2', colnames(St))]
  Vs1<- St[,grep('Vs1', colnames(St))]
  Vs2<- St[,grep('Vs2', colnames(St))]
  
  beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2# q depends on transmission type (whether depends on population density or not); c2 is the contact matrix
  
  lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
  for (t in 1:t0){
    lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+parmset$seed)%*%beta)/sum(St[t,]))}
  
  hosp1 <- c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2, report_infants*0.15, report_infants*0.15, rep(report_children, 2), rep(.001, 6))
  hosp2 <- hosp1 * 0.4
  hosp3 <- c(rep(0.001, 8), rep(report_adults, 4), report_seniors60, report_seniors75)
  
  H1 <- matrix(0, nrow = t0, ncol = al)
  for (i in 1:al) {
    H1[, i] <-
      hosp1[i] * parmset$RRHm * parmset$sigmaM * M0[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHn1 * parmset$sigmaM * Mn1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHn2 *  Mn2[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHv1 * parmset$sigmaM * Mv1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHv2 * Mv1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHn1 * N1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHn2 * N2[, i] * lambda1[, i] +
      hosp1[i] * S0[, i] * lambda1[, i] +
      hosp1[i] * Si[, i] * lambda1[, i] +
      hosp2[i] * parmset$sigma1 * S1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma2 * S2[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * S3[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs2[, i] * lambda1[, i]
  }
  
  H <- rowSums(H1)
  
  H2 <- cbind(H1[,1],H1[,2],H1[,3],H1[,4],H1[,5],H1[,6], rowSums(H1[, 7:8]), rowSums(H1[, 9:12]), H1[, 13],H1[,14])
  age_dist2 <- colSums(H2)
  
  
  LLall <- sum(dpois(x = dat, lambda = H, log = TRUE))
  LLmulti <- dmultinom(x = age_dist2, prob = age_dist, log = TRUE)
  
  LL <- LLall + LLmulti
  
  
  return(LL)
}

fitLL <- optim(par = c(2.3,-1.8,2,-2,-4,-8,-6,-4),
                fn = fitmodel, 
                dat = time_series_pre,  
                control = list(fnscale=-1, maxit=5000))

saveRDS(fitLL,"DATA/fitted_parameters_part1.rds")

baseline.txn.rate = exp(fitLL$par[1])
b1 <-  exp(fitLL$par[2])
phi=(2*pi*(exp(fitLL$par[3]))) / (1+exp(fitLL$par[3]))
report_infants <- 1 / (1 + exp(-fitLL$par[4]))
report_children <- 1 / (1 + exp(-fitLL$par[5]))
report_adults <- 1 / (1 + exp(-fitLL$par[6]))
report_seniors60 <- 1 / (1 + exp(-fitLL$par[7]))
report_seniors75 <- 1 / (1 + exp(-fitLL$par[8]))

#run and plot fitted parameters
results <- ode(y=yinit.vector, method = "ode45", times=fit_times,
               func=MSIRS_immunization_dynamics,
               parms=c(parmset,baseline.txn.rate=baseline.txn.rate,b1=b1,phi=phi))

t0 <- length(time_series_pre)
al <- nrow(yinit)
output <- tail(results,t0)
St <- output[,-1]
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
M0 <- St[,grep('M0', colnames(St))]
Si<- St[,grep('Si', colnames(St))]
Mn1<- St[,grep('Mn1', colnames(St))]
Mn2<- St[,grep('Mn2', colnames(St))]
Mv1<- St[,grep('Mv1', colnames(St))]
Mv2<- St[,grep('Mv2', colnames(St))]
N1<- St[,grep('N1', colnames(St))]
N2<- St[,grep('N2', colnames(St))]
Vs1<- St[,grep('Vs1', colnames(St))]
Vs2<- St[,grep('Vs2', colnames(St))]


beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2# q depends on transmission type (whether depends on population density or not); c2 is the contact matrix

lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
for (t in 1:t0){
  lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+parmset$seed)%*%beta)/sum(St[t,]))}

hosp1 <- c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2, report_infants*0.15, report_infants*0.15, rep(report_children, 2), rep(.001, 6))
hosp2 <- hosp1 * 0.4
hosp3 <- c(rep(0.001, 8), rep(report_adults, 4), report_seniors60, report_seniors75)

H1 <- matrix(0, nrow = t0, ncol = al)
for (i in 1:al) {
  H1[, i] <-
    hosp1[i] * parmset$RRHm * parmset$sigma3 * M0[, i] * lambda1[, i] +
    hosp1[i] * parmset$RRHm * parmset$RRHn1 * parmset$sigma3 * Mn1[, i] * lambda1[, i] +
    hosp1[i] * parmset$RRHm * parmset$RRHn2 *  Mn2[, i] * lambda1[, i] +
    hosp1[i] * parmset$RRHm * parmset$RRHv1 * parmset$sigma3 * Mv1[, i] * lambda1[, i] +
    hosp1[i] * parmset$RRHm * parmset$RRHv2 * Mv1[, i] * lambda1[, i] +
    hosp1[i] * parmset$RRHn1 * N1[, i] * lambda1[, i] +
    hosp1[i] * parmset$RRHn2 * N2[, i] * lambda1[, i] +
    hosp1[i] * S0[, i] * lambda1[, i] +
    hosp1[i] * Si[, i] * lambda1[, i] +
    hosp2[i] * parmset$sigma1 * S1[, i] * lambda1[, i] +
    hosp3[i] * parmset$sigma2 * S2[, i] * lambda1[, i] +
    hosp3[i] * parmset$sigma3 * S3[, i] * lambda1[, i] +
    hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs1[, i] * lambda1[, i] +
    hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs2[, i] * lambda1[, i]
}

H <- rowSums(H1)
plot(H)

H2 <- cbind(H1[,1],H1[,2],H1[,3],H1[,4],H1[,5],H1[,6], rowSums(H1[, 7:8]), rowSums(H1[, 9:12]), H1[, 13],H1[,14])
age_dist2 <-round(colSums(H2)/sum(H2),2)
age_dist
age_dist2


dates1 = seq(from=as.Date('2017-07-02'),to=as.Date('2020-03-31'),by='week')

plotA = ggplot()+
  theme_bw()+
  geom_area(aes(x=dates1,y=time_series_pre/22.65,fill='Data'))+
  geom_line(aes(x=dates1,y=H/22.65,color="Model Fit"),linewidth=1)+
  labs(x=NULL, y="RSV Hospitalization Rate")+
  scale_fill_manual(name=NULL,values='grey80')+
  scale_color_manual(name=NULL, values=c("blue"))+
  scale_y_continuous(name="RSV Hospitalization Rate",sec.axis = sec_axis( trans=~.*10, name="% of Contacts"))
plotA

#attack rates 
infections=matrix(0,nrow=t0,ncol=al)#Number of infections by age
for (i in 1:al){
  infections[,i]=
    parmset$sigma3*M0[,i]*lambda1[,i]+
    S0[,i]*lambda1[,i]+
    parmset$sigma1*S1[,i]*lambda1[,i]+
    parmset$sigma2*S2[,i]*lambda1[,i]+
    parmset$sigma3*S3[,i]*lambda1[,i]}
inf_dist = colSums(infections[118:144,])
inf_dist
pop = rowSums(parmset$yinit.matrix)
attack = inf_dist/pop
attack



# Fit pandemic period  ----------------------------------------------------
fit_times2 <- seq(1, length(time_series_full) + 104, by = 1)

fitpand <-  function(parameters,dat) {
  
  npi1 <- 1 / (1 + exp(-parameters[1]))  
  npi2 <- 1 / (1 + exp(-parameters[2])) 
  npi3 <- 1 / (1 + exp(-parameters[3])) 
  #npi4 <- 1 / (1 + exp(-parameters[4])) 
 
  introductions = data.frame(intros=c(rep(parmset$seed,248),rep(0,44),rep(NA,24),rep(parmset$seed,152))) %>% 
    mutate(intros = na_interpolation(intros, method="linear"))
  introductions = introductions$intros
  parmset$introductions = introductions
  
  #reductions in contacts during the pandemic 
  npi = data.frame(npis=c(rep(1,248),rep(npi1,52),rep(NA,22),rep(npi2,13),rep(npi3,18),rep(NA,13),rep(1,152)))%>% 
    mutate(npis= na_interpolation(npis, method="linear"))
  npi=npi$npis
  parmset$npi = npi
  
  
  # Simulate the model with initial conditions and timesteps defined above, and parameter values from function call
  results <- ode(y=yinit.vector, method = "ode45", times=fit_times2,
                 func=MSIRS_immunization_dynamics,
                 parms=c(parmset,baseline.txn.rate=baseline.txn.rate,b1=b1,phi=phi))
  
  t0 <- length(time_series_full)
  al <- nrow(yinit)
  output <- tail(results,t0)
  St <- output[,-1]
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
  M0 <- St[,grep('M0', colnames(St))]
  Si<- St[,grep('Si', colnames(St))]
  Mn1<- St[,grep('Mn1', colnames(St))]
  Mn2<- St[,grep('Mn2', colnames(St))]
  Mv1<- St[,grep('Mv1', colnames(St))]
  Mv2<- St[,grep('Mv2', colnames(St))]
  N1<- St[,grep('N1', colnames(St))]
  N2<- St[,grep('N2', colnames(St))]
  Vs1<- St[,grep('Vs1', colnames(St))]
  Vs2<- St[,grep('Vs2', colnames(St))]
  
  beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2# q depends on transmission type (whether depends on population density or not); c2 is the contact matrix
  
  lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
  for (t in 1:t0){
    lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+parmset$seed)%*%beta)/sum(St[t,]))}
  
  hosp1 <- c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2, report_infants*0.15, report_infants*0.15, rep(report_children, 2), rep(.001, 6))
  hosp2 <- hosp1 * 0.4
  hosp3 <- c(rep(0.001, 8), rep(report_adults, 4), report_seniors60, report_seniors75)
  
  H1 <- matrix(0, nrow = t0, ncol = al)
  for (i in 1:al) {
    H1[, i] <-
      hosp1[i] * parmset$RRHm * parmset$sigmaM * M0[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHn1 * parmset$sigmaM * Mn1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHn2 *  Mn2[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHv1 * parmset$sigmaM * Mv1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHv2 * Mv1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHn1 * N1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHn2 * N2[, i] * lambda1[, i] +
      hosp1[i] * S0[, i] * lambda1[, i] +
      hosp1[i] * Si[, i] * lambda1[, i] +
      hosp2[i] * parmset$sigma1 * S1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma2 * S2[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * S3[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs2[, i] * lambda1[, i]
  }
  
  H <- rowSums(H1)
  
  LLall <- sum(dpois(x = dat, lambda = H, log = TRUE))
 
  LL <- LLall 
  
  
  return(LL)
}

pandLL <- optim(par = c(1,1,1),
               fn = fitpand,  # the distance function to optimize
               dat = time_series_full,  # the dataset to fit to (dpois function)
               control = list(fnscale=-1, maxit=5000))#
saveRDS(pandLL,"DATA/fitted_parameters_part2.rds")

npi1 <- 1 / (1 + exp(-pandLL$par[1]))  
npi2 <- 1 / (1 + exp(-pandLL$par[2])) 
npi3 <- 1 / (1 + exp(-pandLL$par[3])) 


introductions = data.frame(intros=c(rep(parmset$seed,248),rep(0,44),rep(NA,24),rep(parmset$seed,152))) %>% 
  mutate(intros = na_interpolation(intros, method="linear"))
introductions = introductions$intros
parmset$introductions = introductions

npi = data.frame(npis=c(rep(1,248),rep(npi1,52),rep(NA,22),rep(npi2,13),rep(npi3,18),rep(NA,13),rep(1,152)))%>% 
  mutate(npis= na_interpolation(npis, method="linear"))
npi=npi$npis
parmset$npi = npi


# Simulate the model with initial conditions and timesteps defined above, and parameter values from function call
results <- ode(y=yinit.vector, method = "ode45", times=fit_times2,
               func=MSIRS_immunization_dynamics,
               parms=c(parmset,baseline.txn.rate=baseline.txn.rate,b1=b1,phi=phi))

t0 <- length(time_series_full)
al <- nrow(yinit)
output <- tail(results,t0)
St <- output[,-1]
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
M0 <- St[,grep('M0', colnames(St))]
Si<- St[,grep('Si', colnames(St))]
Mn1<- St[,grep('Mn1', colnames(St))]
Mn2<- St[,grep('Mn2', colnames(St))]
Mv1<- St[,grep('Mv1', colnames(St))]
Mv2<- St[,grep('Mv2', colnames(St))]
N1<- St[,grep('N1', colnames(St))]
N2<- St[,grep('N2', colnames(St))]
Vs1<- St[,grep('Vs1', colnames(St))]
Vs2<- St[,grep('Vs2', colnames(St))]

beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2# q depends on transmission type (whether depends on population density or not); c2 is the contact matrix

lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection

for (t in 1:t0){
  lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+parmset$seed)%*%beta)/sum(St[t,]))}

hosp1 <- c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2, report_infants*0.15, report_infants*0.15, rep(report_children, 2), rep(.001, 6))
hosp2 <- hosp1 * 0.4

hosp3 <- c(rep(0.001, 8), rep(report_adults, 4), report_seniors60, report_seniors75)

H1 <- matrix(0, nrow = t0, ncol = al)
for (i in 1:al) {
  H1[, i] <-
    hosp1[i] * parmset$RRHm * parmset$sigmaM * M0[, i] * lambda1[, i] +
    hosp1[i] * parmset$RRHm * parmset$RRHn1 * parmset$sigmaM * Mn1[, i] * lambda1[, i] +
    hosp1[i] * parmset$RRHm * parmset$RRHn2 *  Mn2[, i] * lambda1[, i] +
    hosp1[i] * parmset$RRHm * parmset$RRHv1 * parmset$sigmaM * Mv1[, i] * lambda1[, i] +
    hosp1[i] * parmset$RRHm * parmset$RRHv2 * Mv1[, i] * lambda1[, i] +
    hosp1[i] * parmset$RRHn1 * N1[, i] * lambda1[, i] +
    hosp1[i] * parmset$RRHn2 * N2[, i] * lambda1[, i] +
    hosp1[i] * S0[, i] * lambda1[, i] +
    hosp1[i] * Si[, i] * lambda1[, i] +
    hosp2[i] * parmset$sigma1 * S1[, i] * lambda1[, i] +
    hosp3[i] * parmset$sigma2 * S2[, i] * lambda1[, i] +
    hosp3[i] * parmset$sigma3 * S3[, i] * lambda1[, i] +
    hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs1[, i] * lambda1[, i] +
    hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs2[, i] * lambda1[, i]
}

H <- rowSums(H1)
plot(H)

dates1 = seq(from=as.Date('2017-07-02'),to=as.Date('2023-09-17'),by='week')


plotA = ggplot()+
  theme_bw()+
  #geom_area(aes(x=dates1,y=time_series_full/22.65,fill='Data'))+
  geom_point(aes(x=dates1,y=time_series_full/22.65,color='Data'),shape=1,size=3)+
  geom_line(aes(x=dates1,y=H/22.65,color="Model Fit"),linewidth=1)+
  geom_line(aes(x=dates1, y=npi[105:429]*10,color="% Contacts"))+
  labs(x=NULL, y="RSV Hospitalization Rate")+
  scale_fill_manual(name=NULL,values='grey80')+
  scale_color_manual(name=NULL, values=c("black","grey40","blue"))+
  scale_y_continuous(name="RSV Hospitalization Rate",sec.axis = sec_axis( trans=~.*10, name="% of Contacts"))
plotA
ggsave(plot=plotA,"FIGURES/eFigure6.png",height=5,width=9,units="in")




# Add noise to fitted parameters and resample -----------------------------

parms = c(baseline.txn.rate,b1,phi, report_seniors60,report_seniors75,report_infants,report_children,report_adults,npi1, npi2, npi3)

conf_int = cbind(parms*.9,parms*1.1)
#Less noise for Phi parameter which is more sensitive to small changes 
conf_int[3,2] = phi*1.025
conf_int[3,1] = phi*0.975
conf_int


h=100
lhs<-maximinLHS(h,11)
new_parms <- cbind(
  beta = lhs[,1]*(conf_int[1,2]-conf_int[1,1])+conf_int[1,1],
  b1 = lhs[,2]*(conf_int[2,2]-conf_int[2,1])+conf_int[2,1],
  phi = lhs[,3]*(conf_int[3,2]-conf_int[3,1])+conf_int[3,1],
  RS60 = lhs[,4]*(conf_int[4,2]-conf_int[4,1])+conf_int[4,1],
  RS75 = lhs[,5]*(conf_int[5,2]-conf_int[5,1])+conf_int[5,1],
  RI = lhs[,6]*(conf_int[6,2]-conf_int[6,1])+conf_int[6,1],
  RC = lhs[,7]*(conf_int[7,2]-conf_int[7,1])+conf_int[7,1],
  RA = lhs[,8]*(conf_int[8,2]-conf_int[8,1])+conf_int[8,1],
  npi1 = lhs[,9]*(conf_int[9,2]-conf_int[8,1])+conf_int[9,1],
  npi2 = lhs[,10]*(conf_int[10,2]-conf_int[9,1])+conf_int[10,1],
  npi3 = lhs[,11]*(conf_int[11,2]-conf_int[10,1])+conf_int[11,1])
 
saveRDS(new_parms,"DATA/fitted_parameters_100.rds")


h=1000
lhs<-maximinLHS(h,11)
new_parms <- cbind(
  beta = lhs[,1]*(conf_int[1,2]-conf_int[1,1])+conf_int[1,1],
  b1 = lhs[,2]*(conf_int[2,2]-conf_int[2,1])+conf_int[2,1],
  phi = lhs[,3]*(conf_int[3,2]-conf_int[3,1])+conf_int[3,1],
  RS60 = lhs[,4]*(conf_int[4,2]-conf_int[4,1])+conf_int[4,1],
  RS75 = lhs[,4]*(conf_int[5,2]-conf_int[5,1])+conf_int[5,1],
  RI = lhs[,4]*(conf_int[6,2]-conf_int[6,1])+conf_int[6,1],
  RC = lhs[,5]*(conf_int[7,2]-conf_int[7,1])+conf_int[7,1],
  RA = lhs[,6]*(conf_int[8,2]-conf_int[8,1])+conf_int[8,1],
  npi1 = lhs[,6]*(conf_int[9,2]-conf_int[8,1])+conf_int[9,1],
  npi2 = lhs[,6]*(conf_int[10,2]-conf_int[9,1])+conf_int[10,1],
  npi3 = lhs[,6]*(conf_int[11,2]-conf_int[10,1])+conf_int[11,1])

saveRDS(new_parms,"DATA/fitted_parameters_1000.rds")

