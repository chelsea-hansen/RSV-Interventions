rm(list=ls())
library(deSolve)
library(tidyverse)
library(imputeTS)
library(lhs)
library(bbmle)

source("MSIRS_immunization_dynamics.R")

data = readRDS("DATA/fixed_parameters.rds")
parmset = data[[1]]
yinit=data[[2]]
yinit.vector=data[[3]]


#RSV data not currently publicly available - will post when possible
#time_series = readRDS("DATA/rsv time series.rds") %>% filter(date<'2023-10-01') %>% 
  #mutate(hosp = round(rsv_rate/100000*2267000))

#time_series = c(time_series$hosp)
#plot(time_series)
#age_dist_all = readRDS("DATA/age_distributions.rds") %>% 
 # ungroup() %>% 
  #filter(season=="pre-pandemic") 
#age_dist = c(age_dist_all$prop)
#age_dist


parmset$monoclonal_01 =  rep(0, length(time_series) + 104+52)
parmset$monoclonal_23 =  rep(0, length(time_series) + 104+52)
parmset$monoclonal_45 =  rep(0, length(time_series) + 104+52)
parmset$monoclonal_67 =  rep(0, length(time_series) + 104+52)
parmset$maternal_vax <- rep(0, length(time_series) + 104+52)
parmset$senior_vax <- rep(0, length(time_series) + 104+52)


# MLE Function ------------------------------------------------------------


fit_bbmle <- function(protrans, amp, trans, RS, RI, RC, RA, npi1, npi2, npi3, npi4, parmset, yinit.vector, age_dist, time_series, yinit) {
  
  fit_times <- seq(1, length(time_series) + 104, by = 1)
  baseline.txn.rate <- exp(protrans)
  b1 <- exp(amp)
  phi <- (2 * pi * (exp(trans))) / (1 + exp(trans))
  
  report_seniors <- 1 / (1 + exp(-RS))
  report_infants <- 1 / (1 + exp(-RI))
  report_children <- 1 / (1 + exp(-RC))
  report_adults <- 1 / (1 + exp(-RA))
  
  npi1f = 1/(1+exp(-npi1))
  npi2f = 1/(1+exp(-npi2))
  
  npi3f = 1/(1+exp(-npi3))
  npi4f = 1/(1+exp(-npi4))
 
  #fit a new reporting rate during the pandemic 
  
  
  introductions = data.frame(intros=c(rep(parmset$seed,248),rep(0,44),rep(NA,24),rep(parmset$seed,114))) %>% 
    mutate(intros = na_interpolation(intros, method="linear"))
  introductions = introductions$intros
  parmset$introductions = introductions
  
  npi = data.frame(npis=c(rep(1,248),rep(npi1f,13),rep(npi2f,39),rep(NA,26),rep(npi3f,13),rep(npi4f,26),rep(NA,13),rep(1,52)))%>% 
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
  
  return(-LL) # Return negative log-likelihood for mle2
}

start_params <- list(protrans = 2.12, amp = -2.14, trans = 2.04, RS=-5.93,RI=-2.42,RC=-4.23,RA=-8.42,npi1=.23,npi2=.84,npi3=.99,npi4=.31)
lower_bound = list(protrans=2,amp=-2.4,trans=1.7,RS=-10,RI=-10,RC=-10,RA=-10,npi1=-10,npi2=-10,npi3=-10,npi4=-10)
upper_bound = list(protrans=2.9,amp=-1.2,trans=2.10,RS=-2,RI=-2,RC=-2,RA=-2,npi1=10,npi2=10,npi3=10,npi4=10)

# Fit the model using mle2
fit <- mle2(minuslogl = fit_bbmle, 
            start = start_params, 
            lower = lower_bound,
            upper=upper_bound,
            method="L-BFGS-B",
            control=list(maxit=5000),
            data = list(parmset = parmset, yinit.vector = yinit.vector, age_dist = age_dist, time_series = time_series, yinit = yinit))

saveRDS(fit,"DATA/fitted_parameters.rds")
estimates = coef(fit)
baseline.txn.rate <- exp(estimates[1])
b1 <- exp(estimates[2])
phi <- (2 * pi * (exp(estimates[3]))) / (1 + exp(estimates[3]))


report_seniors <- 1 / (1 + exp(-estimates[4]))
report_infants <- 1 / (1 + exp(-estimates[5]))
report_children <- 1 / (1 + exp(-estimates[6]))
report_adults <- 1 / (1 + exp(-estimates[7]))

npi1 <- 1 / (1 + exp(-estimates[8]))
npi2 <- 1 / (1 + exp(-estimates[9]))
npi3 <- 1 / (1 + exp(-estimates[10]))
npi4 <- 1 / (1 + exp(-estimates[11]))

introductions = data.frame(intros=c(rep(parmset$seed,248),rep(0,44),rep(NA,24),rep(parmset$seed,114))) %>% 
  mutate(intros = na_interpolation(intros, method="linear"))
introductions = introductions$intros
parmset$introductions = introductions

npi = data.frame(npis=c(rep(1,248),rep(npi1,13),rep(npi2,39),rep(NA,26),rep(npi3,13),rep(npi4,26),rep(NA,13),rep(1,52)))%>% 
  mutate(npis= na_interpolation(npis, method="linear"))
npi=npi$npis
parmset$npi = npi

fit_times <- seq(1, length(time_series) + 104, by = 1)
results <- ode(y = yinit.vector, method = "ode45", times = fit_times, 
               func = MSIRS_immunization_dynamics, 
               parms = c(parmset, baseline.txn.rate = baseline.txn.rate, b1 = b1, phi = phi))

t0 <- length(time_series)
al <- nrow(yinit)
output <- tail(results, t0)
St <- output[, -1]
I1 <- St[, grep('I1', colnames(St))]
I2 <- St[, grep('I2', colnames(St))]
I3 <- St[, grep('I3', colnames(St))]
I4 <- St[, grep('I4', colnames(St))]
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
hosp2 = hosp1*.4
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

Ha <- rowSums(H1)
plot(Ha)
dates1 = seq(from=as.Date('2017-07-02'),to=as.Date('2023-09-30'),by='week')
time_seriesa = time_series

plot(time_seriesa)
plotA = ggplot()+
  theme_bw()+
  geom_area(aes(x=dates1,y=time_seriesa/22.65),fill="grey")+
  geom_line(aes(x=dates1,y=Ha/22.65),color="blue",size=1)+
  geom_line(aes(x=dates1,y=npi[105:430]*10),color="black",size=1)+
  labs(x=NULL, y="RSV Hospitalization Rate")#+
#  scale_y_continuous(name="Hosp Rate",sec.axis = sec_axis( trans=~.*10, name="% of Contacts"))
plotA
#ggsave(plot=plotA,"eFigure6.png",height=5,width=9,units="in")
  
H2= cbind(H1[,1],H1[,2],H1[,3],H1[,4],H1[,5],H1[,6],rowSums(H1[,7:8]),rowSums(H1[,9:12]),H1[,13])
aged = colSums(H2[1:144,])/sum(H2[1:144,])
round(aged,2)
age_dist

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




# Add confidence intervals with Latin Hypercube Sampling --------------------------------------------
parms = c(baseline.txn.rate,b1,phi, report_seniors,report_infants,report_children,report_adults,npi1, npi2, npi3, npi4)

conf_int = cbind(parms*.9,parms*1.1)
#Less noise for Phi parameter which is more sensitive to small changes 
conf_int[3,2] = phi*1.02
conf_int[3,1] = phi*0.98
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

saveRDS(new_parms,"DATA/fitted_parameters_100.rds")


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

saveRDS(new_parms,"DATA/fitted_parameters_1000.rds")


