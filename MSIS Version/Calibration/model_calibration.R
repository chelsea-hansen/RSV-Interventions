
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

"%notin%" = Negate('%in%')


# Step 1 - upload necessary data ------------------------------------------


yinit <- as.matrix(readRDS('Data/yinit.rds')) # data for 2003 population size 
#start by seeding 1 infection in each age group 
yinit[,1:2]=yinit[,1:2]*0.7 # reset for 1980 population size instead for longer burn-in
p <- sum(yinit) 

#weekly birthrate per 1000 population starting from January 1980
#data was initially annual so used an interpolation step (not included here)
#you you could also just use the annual value for each week 
#columns represent age groups, only the first column has data (0 for all other columns)
birth <-readRDS('Data/birth.rds') %>% select(-date, -period)
birth = as.matrix(birth)

#get a vector of dates 
dates = readRDS('Data/birth.rds') %>% select(date) %>% 
  mutate(mmwr_week(date)) #pre-pandemic fitting seasons are line 1933-2100
dates = data.frame(date=seq(from=as.Date("1980-01-05"), to=as.Date("2025-10-01"), by="weeks"))

#Upload the POLYMOD contact matrix, already aggregated to relevant age groups 
contact <- readRDS('Data/contact_POLYMOD.rds')

#proportion of hospitalizations in each age group 
#data not available to be shared - see "dummy" version for general structure 
#wa_agedist = readRDS('Data/age_distribution_dummy.rds')
#wa_agedist <- readRDS('Data/age_distributions.rds') %>% arrange(age)
wa_agedist <- readRDS('update data/age_distributions.rds') %>% arrange(age)
wa_agedist <- as.matrix(subset(wa_agedist, period=="pre-pandemic",select=c("prop_hosp")))


#weekly time series of RSV hospitalizations (all ages) 
#data not publicly available, see dummy version for structure 
#values between 1 and 9 have been supressed and are interpolated 
#rsv_pre = readRDS('Data/rsv_ts.rds') %>% 
  #mutate(hrsv_smooth = round(rollmean(na_interpolation(hosp_rsv_adj), k=3, align="center",fill=NA)),
         #ersv_smooth = round(rollmean(na_interpolation(ed_rsv_adj), k=3, align="center",fill=NA))) %>% 
  #filter(!is.na(hrsv_smooth),date<='2020-03-28') %>% 
  #select(date, rsvH=hrsv_smooth, rsvED=ersv_smooth)

rsv_pre = readRDS('update data/rsv_ts.rds') %>% 
  select(date, rsvH=hrsv_smooth, rsvED=ersv_smooth) %>% 
  filter(date<='2020-03-28') 

#scale hosps to ED visits, pre-pandemic period 
#data not publicly available, see dummy version for structure 
#scaling = readRDS("Data/scale_ed_to_hosp.rds") %>%
scaling = readRDS("update data/scale_ed_to_hosp.rds") %>%
  filter(period=="pre-pandemic") %>% 
  select(scale) 
scaling = c(scaling$scale)


N_ages <- nrow(yinit) 
agenames <- c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y","2-4Y","5-9Y","10-19Y","20-39Y","40-64Y","65Y+")
#note that some of these age groups are collapsed later 
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

#Relative infectiousness for 2nd and subsequent infections
rho1 = 0.75
rho2 = 0.51

#Duration in days
dur.days1 <- 10 #days
dur.days2 <- 7 #days
dur.days3 <- 5 #days


WidthAgeClassMonth = c(rep(2,times=6), 12,12*3,  60, 120, 240, 300, 180)  #Aging rate=1/width age class (months) Vector of long N_age
um=-0.0002227 #death rate - use as starting point, might need to calibrate this 

PerCapitaBirthsYear=birth #birth rate 

#relative risk of infection following 1st, 2nd, infecitons 
sigma1=0.76
sigma2=0.6
sigma3=0.4

#length.step=30.44 #days per month
length.step = 7

#recovery rate 
gamma1= 1/(dur.days1/length.step)  #converts 1/days to 1/length.step
gamma2= 1/(dur.days2/length.step)  
gamma3= 1/(dur.days3/length.step)  


start_time = 1 
tmax = 2100 #end of pre-pandemic period (2020-03-28)
wa_times <- seq(start_time, tmax, by = 1) # gives a sequence of weeks 


## Read in transmission dynamic model
source("Calibration/rsv_dynamics.R")


# Step 2 - Use maximum likelihood estimation to calibrate parameters ------------------------------------------

parmset<-list(PerCapitaBirthsYear=PerCapitaBirthsYear,
              WidthAgeClassMonth=WidthAgeClassMonth,
              um=um,
              rho1=rho1,
              rho2=rho2,
              dur.days1=dur.days1,
              dur.days2=dur.days2,
              dur.days3=dur.days3,
              yinit.matrix=yinit,
              q=1,
              c2=contact,
              sigma1=sigma1,
              sigma2=sigma2,
              sigma3=sigma3,
              time.step='week'
)

fitmodel <-  function(parameters,dat) {  
  protrans <- parameters[1] # parameter for beta 
  amp <- parameters[2] # parameter for seasonal amplitude 
  trans <- parameters[3] # parameter for timing of peak
  DMD <- parameters[4] # parameter for duration of maternal immunity 
  report_seniors <- exp(parameters[5]) # proportion of infections leading to reported hospitalizations in seniors
  report_infants <- exp(parameters[6]) #proportion of infections leading to reported hospitalizations in infants <2m
  I_ex <- exp(parameters[7])#number of weekly imported infections 
  RRHm <- exp(parameters[8])#relative risk of hospitalization for M compartment 
  b1 <- exp(amp) #ensure positive 
  baseline.txn.rate <- exp(protrans) #ensure positive
  phi <-  (2*pi*(exp(trans))) / (1+exp(trans)) # transform to its scale in the model
  durx <- exp(DMD) #ensure positive
  
  # Simulate the model with initial conditions and timesteps defined above, and parameter values from function call
  results <- ode(y=yinit.vector, method = "ode45", t=wa_times,  
                 func=rsv_dynamics, 
                 parms=c(parmset,baseline.txn.rate=baseline.txn.rate,b1=b1,phi=phi, I_ex=I_ex, DurationMatImmunityDays=durx))
  
  t0 <- nrow(rsv_pre)
  burnN <- 1932
  St <- results[-c(1:burnN),-1]
  I1 <- St[,grep('I1', colnames(St))]
  I2 <- St[,grep('I2', colnames(St))]
  I3 <- St[,grep('I3', colnames(St))]
  I4 <- St[,grep('I4', colnames(St))]
  S1 <- St[,grep('S1', colnames(St))]
  S2 <- St[,grep('S2', colnames(St))]
  S3 <- St[,grep('S3', colnames(St))]
  S0 <- St[,grep('S0', colnames(St))]
  M <- St[,grep('M', colnames(St))]
  
  beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(parmset$yinit.matrix)^(1-parmset$q))*parmset$c2# q depends on transmission type (whether depends on population density or not); c2 is the contact matrix
  
  lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
  for (t in 1:t0)
  {lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+rho1*I2[t,]+rho2*I3[t,]+rho2*I4[t,]+I_ex)%*%beta)/sum(St[t,]))}
  
 #scaling reporting fractions to other age groups 
  hosp = c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2,report_infants*0.15, report_infants*0.15, report_seniors,report_seniors*0.06,report_seniors)
 
   H1=matrix(0,nrow=t0,ncol=al)#Number of hospitalizations by age
  for (i in 1:al){
    H1[,i]=
      RRHm*M[,i]*lambda1[,i]+
      S0[,i]*lambda1[,i]+
      sigma1*S1[,i]*lambda1[,i]+
      sigma2*S2[,i]*lambda1[,i]+
      sigma3*S3[,i]*lambda1[,i]}
  
   #combine middle age groups 
  H2= cbind(H1[,1],H1[,2],H1[,3],H1[,4],H1[,5],H1[,6],rowSums(H1[,7:8]),rowSums(H1[,9:12]),H1[,13])
  H2 = t(t(H2)*hosp) #scale infections to hospitalizations 
  H = rowSums(H2) #combine into single time series 
  
  age_dist <- colSums(H2) # age distribution 
  
  LLall <- sum(dpois(x = dat, lambda =H, log = TRUE)) # fit to timeseries
  LLmulti <- dmultinom(x= age_dist,prob = wa_agedist,log = T) # fit to age distribution
  
  #total log likelihood
  LL <- LLall+LLmulti
  return(LL)
}


# Run optimization function  --------------------------------------

fitLL <- optim(par = c(1.1,-1.6,-.1,4,-5,-2,4, -1.6),#starting values 
               fn = fitmodel,        # the distance function to optimize
               dat = rsv_pre$rsvH,  # the dataset to fit to (dpois function)
               control = list(fnscale=-1)) # the log likelihood is negative; here we maximize the log likelihood

#save your parameters 
saveRDS(fitLL, "parameters_27Oct2023.rds")
fitLL = readRDS("Calibration/parameters_27Oct2023.rds")

#list of fit parameters 
baseline.txn.rate=exp(fitLL$par[1])
b1=exp(fitLL$par[2])
phi=(2*pi*(exp(fitLL$par[3]))) / (1+exp(fitLL$par[3]))
DurationMatImmunityDays = exp(fitLL$par[4])
report_seniors <-  exp(fitLL$par[5])
report_infants <-  exp(fitLL$par[6])
I_ex = exp(fitLL$par[7])
RRHm = exp(fitLL$par[8])


# Plot results with fit parameters   ---------------------------------------------

output <- ode(y=yinit.vector, t=wa_times,method = "ode45",
              func=rsv_dynamics, 
              parms=c(parmset,
                      baseline.txn.rate=baseline.txn.rate,
                      b1=b1,
                      phi=phi,
                      DurationMatImmunityDays=DurationMatImmunityDays,
                      report_seniors=report_seniors,
                      report_infants = report_infants,
                      I_ex = I_ex,
                      RRHm=RRHm))

t0 <- 168
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

beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(parmset$yinit.matrix)^(1-parmset$q))*parmset$c2

lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
for (t in 1:t0){
  lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+rho1*I2[t,]+rho2*I3[t,]+rho2*I4[t,]+rho2*I_ex)%*%beta)/sum(St[t,]))}

hosp = c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2,report_infants*0.15, report_infants*0.15, report_seniors,report_seniors*0.06,report_seniors)
ed = hosp*scaling #scale from hospitalizations to ED visits 

H1=matrix(0,nrow=t0,ncol=al)#Number of hospitalizations by age
for (i in 1:al){
  H1[,i]=RRHm*M[,i]*lambda1[,i]+
    S0[,i]*lambda1[,i]+
    sigma1*S1[,i]*lambda1[,i]+
    sigma2*S2[,i]*lambda1[,i]+
    sigma3*S3[,i]*lambda1[,i]}

H2= cbind(H1[,1],H1[,2],H1[,3],H1[,4],H1[,5],H1[,6],rowSums(H1[,7:8]),rowSums(H1[,9:12]),H1[,13])

#scale to hospitalizations 
Hosp = t(t(H2)*hosp)

#scale to ED visits 
ED = t(t(H2)*ed)

age_labels=c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1-4Y","5-64Y","65Y+")
age_dist_hosp = data.frame(age_dist = colSums(Hosp)/sum(Hosp),age=age_labels)
age_dist_ed = data.frame(age_dist = colSums(ED)/sum(ED),age=age_labels)

H = rowSums(Hosp)
H <- data.frame(H)
H$date <- dates$date[1933:2100]

E = rowSums(ED)
E <- data.frame(E)
E$date <- dates$date[1933:2100]


# Plot results  -----------------------------------------------------------

# hospitalization time series 
plot1_h = ggplot()+
  theme_bw()+
  geom_area(data=rsv_pre, aes(x=date, y=rsvH),fill="lightsteelblue")+
  geom_line(data=H, aes(x=date, y=H),color="lightpink4",linewidth=2)+
  labs(x=NULL, y="RSV Hospitalizations")
plot1_h

#hospitalization age distribution 
age_dist_hosp2 = age_dist_hosp %>% select(age_dist, age) %>% mutate(data="model")
wa_agedist2 = as.data.frame(wa_agedist) %>% select("age_dist"=prop_hosp) %>% mutate(age=age_labels, data="data")
age_hosp_plot_dat = rbind(age_dist_hosp2,wa_agedist2) %>% mutate(age = factor(age, levels=age_labels))

plot2_h = ggplot(data=age_hosp_plot_dat)+
  theme_bw()+
  geom_bar(aes(x=age, y=age_dist,fill=data),stat="identity")+
  facet_grid(rows=vars(data))+
  labs(x=NULL, y="Age Distribution")+
  scale_fill_manual(name=NULL, values=c("lightsteelblue","lightpink4"))
plot2_h
h_plots = plot_grid(plot1_h, plot2_h, nrow=1)
h_plots
ggsave(plot=h_plots, "hospitalizations.png",height=6,width=12,units="in")
# ED time series 

plot1_ed = ggplot()+
  theme_bw()+
  geom_area(data=rsv_pre, aes(x=date, y=rsvED),fill="lightsteelblue")+
  geom_line(data=E, aes(x=date, y=E),color="lightpink4",linewidth=2)+
  labs(x=NULL, y="RSV ED Visits")
plot1_ed

#ED age distribution 
age_dist_ed2 = age_dist_ed %>% select(age_dist, age) %>% mutate(data="model")
wa_agedist3 = wa_agedist <- readRDS('Data/age_distributions.rds') %>% 
  arrange(age) %>% 
  filter(period=="pre-pandemic") %>% 
  select(age_dist = prop_ed) %>% 
  mutate(age=age_labels, data="data")
age_ed_plot_dat = rbind(age_dist_ed2,wa_agedist3) %>% mutate(age = factor(age, levels=age_labels))

plot2_ed = ggplot(data=age_ed_plot_dat)+
  theme_bw()+
  geom_bar(aes(x=age, y=age_dist,fill=data),stat="identity")+
  facet_grid(rows=vars(data))+
  labs(x=NULL, y="Age Distribution")+
  scale_fill_manual(name=NULL, values=c("lightsteelblue","lightpink4"))
plot2_ed
ed_plots = plot_grid(plot1_ed, plot2_ed, nrow=1)
ed_plots
ggsave(plot=ed_plots, "ED visits.png",height=6,width=12,units="in")
