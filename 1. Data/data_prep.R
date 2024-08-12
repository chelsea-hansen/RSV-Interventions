rm(list=ls())

library(tidyverse)
library(cdcfluview)
library(zoo)
library(imputeTS)
library(tidycensus)

"%notin%" = Negate('%in%')

#This code prepares the demographic data needed to calibrate the model


  #1. Population size 
  #2. birth rates
  #3. net mortality and migration rates to recreate population growth 



# 1. Population size ------------------------------------------------------

#Get population estimates for the most recent population size estimate 
#This will be scaled down to the 1995 population size for the burn-in period 
#you can fill in the initiate_compartments spreadsheet or use the code below 
pop = get_estimates(
  geography = "county",
  product = "characteristics",
  breakdown = c("AGEGROUP"),
  breakdown_labels = TRUE,
  state = "WA",
  county = "King",
  year = 2022)%>% 
  select(AGEGROUP,value) 
new_agegrps = c("all","0to4","5to9",rep("10to19",2),rep("20-39",4),rep("40to59",4),rep("60+",6))

pop = pop %>% 
  mutate(new_age = factor(new_agegrps,levels=unique(new_agegrps))) %>% 
  group_by(new_age) %>% 
  summarize(population = sum(value)) 

#total population in 2022 (most recent data available )
pop_2022 = pop %>% filter(new_age=="all") 
pop_2022=pop_2022$population

#total population in 1995 (start of burn-in period)
pop_1995 = readRDS("1. Data/Demographic Data/birth_rates_by_county.rds") %>% 
  filter(state=="WA",county=="King County",year==1995)
pop_1995=pop_1995$population



#This will look the same as the initiate_compartments spreadsheet 
M = c(rep(pop[2,2]*.034,3),rep(0,10))
S0 = c(rep(0,3),rep(pop[2,2]*.034,3),pop[2,2]*.2,pop[2,2]*.6,pop[3,2],pop[4,2],pop[5,2],pop[6,2],pop[7,2])
I1 = c(rep(0,3),rep(1,10))
R1 = rep(0,13)
S1 = rep(0,13)
I2 = rep(0,13)
R2 = rep(0,13)
S2 = rep(0,13)
I3 = rep(0,13)
R3 = rep(0,13)
S3 = rep(0,13)
I4 = rep(0,13)
R4 = rep(0,13)
Mn = rep(0,13)
Mv = rep(0,13)
N = rep(0,13)
Si = rep(0,13)
Vs1 = rep(0,13)
Vs2 = rep(0,13)
yinit=matrix(unlist(cbind(M,S0,I1,R1,S1,I2,R2,S2,I3,R3,S3,I4,R4,Mn,Mv,N,Si,Vs1,Vs2)),byrow=FALSE,ncol=19)
yinit[,1:2]= yinit[,1:2]*(pop_1995/pop_2022)#adjust population size to 1995 
colnames(yinit) = c("M","S0","I1","R1","S1","I2","R2","S2","I3","R3","S3","I4","R4","Mn","Mv","N","Si","Vs1","Vs2")
saveRDS(yinit,"1. Data/Demographic Data/yinit_interventions.rds") #save the full version 
saveRDS(yinit[,1:13],"1. Data/Demographic Data/yinit.rds")#save a version that does not have the intervention compartments 

#format the data for feeding into the model 
N_ages <- nrow(yinit) 
agenames <- c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y","2-4Y","5-9Y","10-19Y","20-39Y","40-59Y","60Y+")
al <- N_ages
rownames(yinit) <- agenames
yinit.vector <- as.vector(unlist(yinit))
name.array <- array(NA, dim=dim(yinit))
for(i in 1:dim(name.array)[1]){
  for(j in 1:dim(name.array)[2]){
    name.array[i,j] <- paste(dimnames(yinit)[[1]][i],dimnames(yinit)[[2]][j]  )
  }
}
name.vector <- as.vector(name.array)
names(yinit.vector) <- name.vector

saveRDS(yinit.vector, "1. Data/Demographic Data/yinit.vector_interventions.rds")#save the full version 
saveRDS(yinit.vector[1:169], "1. Data/Demographic Data/yinit.vector.rds")#save a version without the intervention compartments 

#seeding 1 infection per 100,000 (based on 2022 population)
seed = pop_2022/100000 #seeding 1 infection per 100,000 (based on 2022 population)


# 2. Birth Rates --------------------------------------------------------
births = readRDS("1. Data/Demographic Data/birth_rates_by_county.rds") %>% 
  filter(state=="WA",county=="King County") %>% 
  mutate(mmwr_year=as.numeric(year),
         mmwr_week=as.numeric(1)) %>% 
  select(mmwr_year,mmwr_week, birth_rate) 

fill_in = data.frame(date =seq(from=as.Date('1995-01-07'), to=as.Date("2025-10-01"), by=7)) %>% 
  mutate(mmwr_week(date))

birth_interp = births %>% 
  right_join(fill_in, by=c("mmwr_year","mmwr_week")) %>% 
  arrange(date) %>% 
  mutate(birth_rate = na_locf(birth_rate))
#for the 2023-2025 data we are just assuming the same as 2022 (most recent available data)

#matrix with column for each age group - for all age groups except the youngest birth rate = 0
birth_complete = data.frame(date = birth_interp$date,
                   V1 = birth_interp$birth_rate, 
                   V2 = rep(0,nrow(birth_interp)),
                   V3 = rep(0,nrow(birth_interp)),
                   V4 = rep(0,nrow(birth_interp)),
                   V5 = rep(0,nrow(birth_interp)),
                   V6 = rep(0,nrow(birth_interp)),
                   V7 = rep(0,nrow(birth_interp)),
                   V8 = rep(0,nrow(birth_interp)),
                   V9 = rep(0,nrow(birth_interp)),
                   V10 = rep(0,nrow(birth_interp)),
                   V11 = rep(0,nrow(birth_interp)),
                   V12 = rep(0,nrow(birth_interp)),
                   V13 = rep(0,nrow(birth_interp)))


saveRDS(birth_complete,"1. Data/Demographic Data/births_kingcounty.rds")



# 3. Net death and migration ----------------------------------------------
#estimated value to use to recreate the population growth
#will check later in the model, this value may need to be adjusted slightly 

weeks = 1461 #number of weeks from January 1995 to December 2022 (check in births data frame)
#weekly population growth divided by the average population size (from start of time series to end of time series)
um = ((pop_2022-pop_1995)/weeks)/((pop_2022+pop_1995)/2)

#save the parameters
other_parms = c("seed"=seed, "um"=um,"pop"=pop_2022)
saveRDS(other_parms,"1. Data/Demographic Data/other_parms.rds")

