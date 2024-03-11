intervention_models <- function(t,y,parms, time.step='week'){
  
  States<-array(y, dim=dim(parms$yinit.matrix))
  dimnames(States) <- dimnames(parms$yinit.matrix)
  um=parms$um
  
  if(parms$time.step=='month'){
    period=12
    length.step=30.44 #days
  }else if(parms$time.step=='week'){
    period=52.1775
    length.step=7 #days
  }
  
  omega = 1/(parms$DurationMatImmunityDays/length.step)
  
  #parameters for monoclonals 
  birth_n = parms$birth_dose[t]#infants who receive birth dose of nirsevimab
  V_n = parms$cover_n[t] #infants who receive nirsevimab catch-up dose  
  waningN=1/(parms$waningN/length.step) #duration of nirsevimab protection 
  RRIn=parms$RRIn #relative risk of infection for infants receiving nirsevimab (default is to set to 1)
 
 
  #parameters for seniors
  V_s=parms$cover_s[t] # indicator function for when the vaccine will be administered * vaccine coverage
  waningS=1/(parms$waningS/2/length.step) #duration of protection from vaccination (divided by 2 because 2 compartments)
  RRIs=parms$RRIs #relative risk of infection for vaccinated seniors
  
  mu= 1/parms$WidthAgeClassMonth
  if(parms$time.step=='week'){
    mu= 1/(WidthAgeClassMonth*4.345)
  }
  
  gamma1= 1/(parms$dur.days1/length.step)  #converts 1/days to 1/length.step
  gamma2= 1/(parms$dur.days2/length.step)  
  gamma3= 1/(parms$dur.days3/length.step)  
  gamma4= gamma3  
  
  #Pull out the states  for the model as vectors
  M <-  States[,'M'] # newborns who do not receive monoclonals 
  
 
  Mn <-  States[,'Mn'] # newborns who receive monoclonals at birth
  N <-  States[,'N'] #infants who receive nirsevimab ahead of the season but not at birth 
  
  S0 <-  States[,'S0']
  I1 <-  States[,'I1']
  
  S1 <-  States[,'S1']
  I2 <-  States[,'I2']
  
  S2 <-  States[,'S2']
  I3 <-  States[,'I3']
  
  S3 <-  States[,'S3']
  I4 <-  States[,'I4']
  
  Vs1 <-  States[,'Vs1'] # vaccination of seniors
  Vs2 <-  States[,'Vs2']
  
  N.ages <- length(M)
  
  ###Check the standardization of beta and overall structure of lambda here
  ####################
  seasonal.txn <- (1+parms$b1*cos(2*pi*(t-parms$phi*period)/period))# seasonality waves
  
  # baseline.txn.rate is the probability of transmission given contact per capita
  # (parms$dur.days1/30.44) is the duration of infectiousness of primary infection
  # q depends on transmission type (whether depends on population density or not)
  # density (q=0) vs frequency-dependent (q=1) transmission
  # c2 is the contact matrix 
  # beta is transimissibility per unit time
  transmission_unittime <-  parms$baseline.txn.rate/(parms$dur.days1/length.step)
  beta=transmission_unittime*parms$c2*parms$npi[t]
  
  beta_a_i <- seasonal.txn * beta/(sum(States)^parms$q)
  infectiousN <- I1 + parms$rho1*I2 + parms$rho2*I3 + parms$rho2*I4 + parms$introductions[t]
  
  lambda <- infectiousN %*% beta_a_i 
  lambda <- as.vector(lambda)
  ####################################  
  
  dy <- matrix(NA, nrow=N.ages, ncol=ncol(States))
  colnames(dy) <- colnames(States)
  
  period.birth.rate <- log(parms$PerCapitaBirthsYear[t,]+1)/period 
  
  #mu represents aging to the next class
  #um is death rate
  
  Aging.Prop <- c(0,mu[1:(N.ages-1)])
  
  # newborns who do not receive nirsevimab at birth (still protected by natural maternal immunity but can get infected)
  dy[,'M'] <- (1-birth_n)*period.birth.rate*sum(States) - 
    lambda*M - 
    V_n*sum(States[1:3,"M"]) - #infants who receive nirsevimab after birth (but still in M compartment)
    (omega+(mu+um))*M +
    Aging.Prop*c(0,M[1:(N.ages-1)]) 
  
  # newborns who receive monoclonals 
  dy[,'Mn'] <- birth_n*period.birth.rate*sum(States) - 
    RRIn*lambda*Mn - 
    (waningN+(mu+um))*Mn + #waning immunity, aging out, migration and death 
    Aging.Prop*c(0,Mn[1:(N.ages-1)]) #aging in 
  
  # infants <8 months who did not receive nirsevimab at birth but receive a catch-up dose ahead of the season
  # these infants can be in the M or S0 compartments 
  dy[,'N'] <- V_n*sum(States[1:3,"M"]) + V_n*sum(States[1:3,"S0"]) - 
    waningN*N-
    RRIn*lambda*N -
    (mu + um)*N + 
    Aging.Prop*c(0,N[1:(N.ages-1)]) 
  
  dy[,'S0'] <- omega*M + waningN*N + waningN*Mn-
    V_n*sum(States[1:3,"S0"]) -
    lambda*S0 -
    (mu + um)*S0 + 
    Aging.Prop*c(0,S0[1:(N.ages-1)]) 
  
  dy[,'I1'] <-   lambda*S0+ lambda*M + RRIn*lambda*N + RRIn*lambda*Mn - 
    (gamma1 + mu + um)*I1 + 
    Aging.Prop*c(0,I1[1:(N.ages-1)]) 
  
  dy[,'S1'] <- gamma1*I1 - 
    parms$sigma1*lambda*S1 - 
    (mu+um)*S1 + 
    Aging.Prop*c(0,S1[1:(N.ages-1)]) 
  
  dy[,'I2'] <- parms$sigma1*lambda*S1 - 
    gamma2*I2-
    (mu + um)*I2 + 
    Aging.Prop*c(0,I2[1:(N.ages-1)]) 
  
  dy[,'S2'] <- gamma2*I2 - 
    parms$sigma2*lambda*S2 -
    (mu+um)*S2 + 
    Aging.Prop*c(0,S2[1:(N.ages-1)]) 
  
  dy[,'I3'] <- parms$sigma2*lambda*S2 -
    (gamma3 + mu+um)*I3 +  
    Aging.Prop*c(0,I3[1:(N.ages-1)]) 
  
  dy[,'S3'] <- gamma3*I3 +  
    gamma4*I4 +
    waningS*Vs2-
    parms$sigma3*lambda*S3 -
    c(rep(0,11),.2*V_s*sum(States[12,"S3"]),V_s*sum(States[13,"S3"])) - #seniors who get vaccinated 
    (mu + um)*S3 + 
    Aging.Prop*c(0,S3[1:(N.ages-1)]) 
  
  
  #vaccination compartments for seniors - using "linear chain trick" to keep people in the vaccinated states longer
  #2 V compartments, 1 for each season that vaccination lasts 
    
    #make a vaccination compartment for the >65 
    dy[,'Vs1'] <- c(rep(0,11),.2*V_s*sum(States[12,"S3"]),V_s*sum(States[13,"S3"])) - 
      RRIs*parms$sigma3*lambda*Vs1 -
      (mu + um)*Vs1 -
      waningS*Vs1 +
      Aging.Prop*c(0,Vs1[1:(N.ages-1)])
    
    #add another vaccination compartment
    dy[,'Vs2'] <- waningS*Vs1 - 
      RRIs*parms$sigma3*lambda*Vs2 -
      (mu + um)*Vs2 -
      waningS*Vs2 +
      Aging.Prop*c(0,Vs2[1:(N.ages-1)])
    
  dy[,'I4'] <- parms$sigma3*lambda*S3 + 
    RRIs*parms$sigma3*lambda*Vs1+
    RRIs*parms$sigma3*lambda*Vs2-
    gamma4*I4 - 
    (mu + um)*I4 + 
    Aging.Prop*c(0,I4[1:(N.ages-1)]) 
  
  
  derivs <- as.vector(dy)
  
  res <- list(derivs)
  
  return(res)
}