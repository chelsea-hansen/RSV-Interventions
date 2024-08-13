rm(list=ls())
library(deSolve)
library(tidyverse)
library(lhs)
library(cowplot)
library(imputeTS)

source("MSIRS_immunization_dynamics.R")


data = readRDS("Data/fixed_parameters.rds")
parmset = data[[1]]
yinit=data[[2]]
yinit.vector=data[[3]]


time_series = readRDS("Data/rsv time series.rds") %>% filter(date<'2023-10-01') %>% 
  mutate(rsv_rate = rsv_rate/100000*2267000)
time_series = round(time_series$rsv_rate)
age_dist = readRDS("Data/age_distributions.rds") %>% filter(season=="pre-pandemic")
age_dist = age_dist$prop


#setting all immunizations to 0 
parmset$monoclonal_01 =  rep(0, length(time_series) + 104)
parmset$monoclonal_23 =  rep(0, length(time_series) + 104)
parmset$monoclonal_45 =  rep(0, length(time_series) + 104)
parmset$monoclonal_67 =  rep(0, length(time_series) + 104)
parmset$maternal_vax <- rep(0, length(time_series) + 104)
parmset$senior_vax <- rep(0, length(time_series) + 104)

fit_times <- seq(1, length(time_series) + 104, by = 1)

# MLE Function ------------------------------------------------------------


fitmodel <- function(parameters) {
  
  
  baseline.txn.rate <- exp(parameters[1])
  b1 <- exp(parameters[2])
  phi <- (2 * pi * (exp(parameters[3]))) / (1 + exp(parameters[3]))
  
  report_seniors <- 1 / (1 + exp(-parameters[4]))
  report_infants <- 1 / (1 + exp(-parameters[5]))
  report_children <- 1 / (1 + exp(-parameters[6]))
  report_adults <- 1 / (1 + exp(-parameters[7]))
  
  npi1 = 1/(1+exp(-parameters[8]))
  npi2 = 1/(1+exp(-parameters[9]))
  npi3 = 1/(1+exp(-parameters[10]))
  npi4 = 1/(1+exp(-parameters[11]))
 
  #fit a new reporting rate during the pandemic 
  
  
  introductions = data.frame(intros=c(rep(parmset$seed,248),rep(0,44),rep(NA,24),rep(parmset$seed,114))) %>% 
    mutate(intros = na_interpolation(intros, method="linear"))
  introductions = introductions$intros
  parmset$introductions = introductions
  
  npi = data.frame(npis=c(rep(1,248),rep(npi1,13),rep(npi2,39),rep(NA,26),rep(npi3,13),rep(npi4,26),rep(NA,13),rep(1,52)))%>% 
    mutate(npis= na_interpolation(npis, method="linear"))
  npi=npi$npis
  parmset$npi = npi
  
  
  results <- ode(y = yinit.vector, method = "ode45", t = fit_times, func = MSIRS_immunization_dynamics, parms = c(parmset, baseline.txn.rate = baseline.txn.rate, b1 = b1, phi = phi))
  
  t0 <- length(time_series)
  al <- nrow(yinit)
  output <- tail(results, t0)
  St <- output[, -1]
  I1 <- St[, grep('I1', colnames(St))]
  I2 <- St[, grep('I2', colnames(St))]
  I3 <- St[, grep('I3', colnames(St))]
  I4 <- St[, grep('I4', colnames(St))]
  R1 <- St[, grep('R1', colnames(St))]
  R2 <- St[, grep('R2', colnames(St))]
  R3 <- St[, grep('R3', colnames(St))]
  R4 <- St[, grep('R4', colnames(St))]
  S1 <- St[, grep('S1', colnames(St))]
  S2 <- St[, grep('S2', colnames(St))]
  S3 <- St[, grep('S3', colnames(St))]
  S0 <- St[, grep('S0', colnames(St))]
  M0 <- St[, grep('M0', colnames(St))]
  Si <- St[, grep('Si', colnames(St))]
  Mn <- St[, grep('Mn', colnames(St))]
  Mv <- St[, grep('Mv', colnames(St))]
  N <- St[, grep('N', colnames(St))]
  Vs1 <- St[, grep('Vs1', colnames(St))]
  Vs2 <- St[, grep('Vs2', colnames(St))]
  
  beta <- baseline.txn.rate / (parmset$dur.days1 / 7) / (sum(yinit)^(1 - parmset$q)) * parmset$c2
  
  lambda1 <- matrix(0, nrow = t0, ncol = 13)
  for (t in 1:t0) {
    lambda1[t, ] <- as.vector((1 + b1 * cos(2 * pi * (t - phi * 52.1775) / 52.1775)) * ((I1[t, ] + parmset$rho1 * I2[t, ] + parmset$rho2 * I3[t, ] + parmset$rho2 * I4[t, ] + parmset$seed) %*% beta) / sum(St[t, ]))
  }
  
  hosp1 <- c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2, report_infants*0.15, report_infants*0.15, rep(report_children, 2), rep(.001, 5))
  hosp2 <- hosp1 * 0.4
  hosp3 <- c(rep(0.001, 8), rep(report_adults, 4), report_seniors)
  
  H1 <- matrix(0, nrow = t0, ncol = al)
  for (i in 1:al) {
    H1[, i] <-
      hosp1[i] * parmset$RRHm * parmset$sigma3 * M0[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$sigma3 * Mn[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$sigma3 * Mv[, i] * lambda1[, i] +
      hosp1[i] * N[, i] * lambda1[, i] +
      hosp1[i] * S0[, i] * lambda1[, i] +
      hosp1[i] * Si[, i] * lambda1[, i] +
      hosp2[i] * parmset$sigma1 * S1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma2 * S2[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * S3[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * Vs1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * Vs2[, i] * lambda1[, i]
  }
  
  H <- rowSums(H1)
  
  H2 <- cbind(H1[, 1], H1[, 2], H1[, 3], H1[, 4], H1[, 5], H1[, 6], rowSums(H1[, 7:8]), rowSums(H1[, 9:12]), H1[, 13])
  age_dist2 <- colSums(H2[1:144,])
  
 
  LLall <- sum(dpois(x = time_series, lambda = H, log = TRUE))
  LLmulti <- dmultinom(x = age_dist2, prob = age_dist, log = TRUE)
  
  LL <- LLall + LLmulti
  
  return(LL) # Return negative log-likelihood for mle2
}


start_params <-  c(2.12,-2.14,2.04,-5.93,-2.42,-4.23,-8.42,0.23,0.84,0.99,0.31)
lower_bound = c(2,-2.4,1.7,-10,-10,-10,-10,-10,-10,-10,-10)
upper_bound = c(2.7,-1.2,2.10,-2,-2,-2,-2,10,10,10,10)



fitLL <- optim(par = start_params,
               fn = fitmodel,        # the distance function to optimize
              # dat = round(time_series$rsv_rate/100000*2267000),  # the dataset to fit to (dpois function)
               lower = lower_bound,
               upper = upper_bound,
               method = "L-BFGS-B",
               control = list(fnscale=-1, maxit=6000,trace=1,REPORT=10))

saveRDS(fitLL,"Data/fit_13Aug2024.rds")
baseline.txn.rate = exp(fitLL$par[1])
b1 = exp(fitLL$par[2])
phi <- (2 * pi * (exp(fitLL$par[3]))) / (1 + exp(fitLL$par[3]))

report_seniors <- 1 / (1 + exp(-fitLL$par[4]))
report_infants <- 1 / (1 + exp(-fitLL$par[5]))
report_children <- 1 / (1 + exp(-fitLL$par[6]))
report_adults <- 1 / (1 + exp(-fitLL$par[7]))

npi1 = 1/(1+exp(-fitLL$par[8]))
npi2 = 1/(1+exp(-fitLL$par[9]))
npi3 = 1/(1+exp(-fitLL$par[10]))
npi4 = 1/(1+exp(-fitLL$par[11]))



introductions = data.frame(intros=c(rep(parmset$seed,248),rep(0,44),rep(NA,24),rep(parmset$seed,114))) %>% 
  mutate(intros = na_interpolation(intros, method="linear"))
introductions = introductions$intros
parmset$introductions = introductions

npi = data.frame(npis=c(rep(1,248),rep(npi1,13),rep(npi2,39),rep(NA,26),rep(npi3,13),rep(npi4,26),rep(NA,13),rep(1,52)))%>% 
  mutate(npis= na_interpolation(npis, method="linear"))
npi=npi$npis
parmset$npi = npi


results <- ode(y = yinit.vector, method = "ode45", t = fit_times, func = MSIRS_immunization_dynamics, parms = c(parmset, baseline.txn.rate = baseline.txn.rate, b1 = b1, phi = phi))

t0 <- length(time_series)
al <- nrow(yinit)
output <- tail(results, t0)
St <- output[, -1]
I1 <- St[, grep('I1', colnames(St))]
I2 <- St[, grep('I2', colnames(St))]
I3 <- St[, grep('I3', colnames(St))]
I4 <- St[, grep('I4', colnames(St))]
R1 <- St[, grep('R1', colnames(St))]
R2 <- St[, grep('R2', colnames(St))]
R3 <- St[, grep('R3', colnames(St))]
R4 <- St[, grep('R4', colnames(St))]
S1 <- St[, grep('S1', colnames(St))]
S2 <- St[, grep('S2', colnames(St))]
S3 <- St[, grep('S3', colnames(St))]
S0 <- St[, grep('S0', colnames(St))]
M0 <- St[, grep('M0', colnames(St))]
Si <- St[, grep('Si', colnames(St))]
Mn <- St[, grep('Mn', colnames(St))]
Mv <- St[, grep('Mv', colnames(St))]
N <- St[, grep('N', colnames(St))]
Vs1 <- St[, grep('Vs1', colnames(St))]
Vs2 <- St[, grep('Vs2', colnames(St))]

beta <- baseline.txn.rate / (parmset$dur.days1 / 7) / (sum(yinit)^(1 - parmset$q)) * parmset$c2

lambda1 <- matrix(0, nrow = t0, ncol = 13)

for (t in 1:t0) {
  lambda1[t, ] <- as.vector((1 + b1 * cos(2 * pi * (t - phi * 52.1775) / 52.1775)) * ((I1[t, ] + parmset$rho1 * I2[t, ] + parmset$rho2 * I3[t, ] + parmset$rho2 * I4[t, ] + parmset$seed) %*% beta) / sum(St[t, ]))
}

hosp1 <- c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2, report_infants*0.15, report_infants*0.15, rep(report_children, 2), rep(.001, 5))
hosp2 <- hosp1 * 0.4
hosp3 <- c(rep(0.001, 8), rep(report_adults, 4), report_seniors)

H1 <- matrix(0, nrow = t0, ncol = al)
for (i in 1:al) {
  H1[, i] <-
    hosp1[i] * parmset$RRHm * parmset$sigma3 * M0[, i] * lambda1[, i] +
    hosp1[i] * parmset$RRHm * parmset$sigma3 * Mn[, i] * lambda1[, i] +
    hosp1[i] * parmset$RRHm * parmset$sigma3 * Mv[, i] * lambda1[, i] +
    hosp1[i] * N[, i] * lambda1[, i] +
    hosp1[i] * S0[, i] * lambda1[, i] +
    hosp1[i] * Si[, i] * lambda1[, i] +
    hosp2[i] * parmset$sigma1 * S1[, i] * lambda1[, i] +
    hosp3[i] * parmset$sigma2 * S2[, i] * lambda1[, i] +
    hosp3[i] * parmset$sigma3 * S3[, i] * lambda1[, i] +
    hosp3[i] * parmset$sigma3 * Vs1[, i] * lambda1[, i] +
    hosp3[i] * parmset$sigma3 * Vs2[, i] * lambda1[, i]
}

H <- rowSums(H1)
plot(H)
H2 <- cbind(H1[, 1], H1[, 2], H1[, 3], H1[, 4], H1[, 5], H1[, 6], rowSums(H1[, 7:8]), rowSums(H1[, 9:12]), H1[, 13])
age_dist2 <- colSums(H2[1:144,])
age_dist2/sum(age_dist2)


plotA = ggplot()+
  theme_bw()+
  geom_area(aes(x=1:326,y=time_series/22.67),fill="grey")+
  geom_line(aes(x=1:326,y=H/22.67),color="blue",size=1)+
  geom_line(aes(x=1:326,y=npi[105:430]*10),color="black",size=1)+
  labs(x=NULL, y="RSV Hospitalization Rate")#+

plotA


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



# Latin Hypercube Sampling  -----------------------------------------------

parms = c(baseline.txn.rate,b1,phi, report_seniors,report_infants,report_children,report_adults,npi1, npi2, npi3, npi4)

conf_int = cbind(parms*.9,parms*1.1)
conf_int[3,2] = phi*1.03
conf_int[3,1] = phi*0.97
conf_int


h=100
lhs<-maximinLHS(h,11)
new_parms <- cbind(
  beta = lhs[,1]*(conf_int[1,2]-conf_int[1,1])+conf_int[1,1],
  b1 = lhs[,2]*(conf_int[2,2]-conf_int[2,1])+conf_int[2,1],
  phi = lhs[,3]*(conf_int[3,2]-conf_int[3,1])+conf_int[3,1],
  RS = lhs[,4]*(conf_int[4,2]-conf_int[4,1])+conf_int[4,1],
  RI = lhs[,4]*(conf_int[5,2]-conf_int[5,1])+conf_int[5,1],
  RC = lhs[,5]*(conf_int[6,2]-conf_int[6,1])+conf_int[6,1],
  RA = lhs[,6]*(conf_int[7,2]-conf_int[7,1])+conf_int[7,1],
  npi1 = lhs[,6]*(conf_int[8,2]-conf_int[8,1])+conf_int[8,1],
  npi2 = lhs[,6]*(conf_int[9,2]-conf_int[9,1])+conf_int[9,1],
  npi3 = lhs[,6]*(conf_int[10,2]-conf_int[10,1])+conf_int[10,1],
  npi4 = lhs[,6]*(conf_int[11,2]-conf_int[11,1])+conf_int[11,1])

saveRDS(new_parms,"Data/parameters_13Aug_100.rds")


h=1000
lhs<-maximinLHS(h,11)
new_parms <- cbind(
  beta = lhs[,1]*(conf_int[1,2]-conf_int[1,1])+conf_int[1,1],
  b1 = lhs[,2]*(conf_int[2,2]-conf_int[2,1])+conf_int[2,1],
  phi = lhs[,3]*(conf_int[3,2]-conf_int[3,1])+conf_int[3,1],
  RS = lhs[,4]*(conf_int[4,2]-conf_int[4,1])+conf_int[4,1],
  RI = lhs[,4]*(conf_int[5,2]-conf_int[5,1])+conf_int[5,1],
  RC = lhs[,5]*(conf_int[6,2]-conf_int[6,1])+conf_int[6,1],
  RA = lhs[,6]*(conf_int[7,2]-conf_int[7,1])+conf_int[7,1],
  npi1 = lhs[,6]*(conf_int[8,2]-conf_int[8,1])+conf_int[8,1],
  npi2 = lhs[,6]*(conf_int[9,2]-conf_int[9,1])+conf_int[9,1],
  npi3 = lhs[,6]*(conf_int[10,2]-conf_int[10,1])+conf_int[10,1],
  npi4 = lhs[,6]*(conf_int[11,2]-conf_int[11,1])+conf_int[11,1])

saveRDS(new_parms,"Data/parameters_13Aug_1000.rds")




