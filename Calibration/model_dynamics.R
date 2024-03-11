MSIRS_dynamics <- function(t,y,parms){
  
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
  waning1 = 1/(parms$recover1/length.step)
  waning2 = 1/(parms$recover2/length.step)
  waning3 = 1/(parms$recover3/length.step)
  waning4 = 1/(parms$recover4/length.step)
  
  
  
  
  mu= 1/parms$WidthAgeClassMonth
  if(parms$time.step=='week'){
    mu= 1/(parms$WidthAgeClassMonth*4.345)
  }
  
  gamma1= 1/(parms$dur.days1/parms$length.step) 
  gamma2= 1/(parms$dur.days2/parms$length.step)  
  gamma3= 1/(parms$dur.days3/parms$length.step)  
  gamma4= gamma3  
  
  #Pull out the states  for the model as vectors
  M <-  States[,'M']
  S0 <-  States[,'S0']
  I1 <-  States[,'I1']
  R1 <-  States[,'R1']
  
  S1 <-  States[,'S1']
  I2 <-  States[,'I2']
  R2 <-  States[,'R2']
  
  S2 <-  States[,'S2']
  I3 <-  States[,'I3']
  R3 <-  States[,'R3']
  
  S3 <-  States[,'S3']
  I4 <-  States[,'I4']
  R4 <-  States[,'R4']
  
  N.ages <- length(M)
  
  ####################
  seasonal.txn <- (1+parms$b1*cos(2*pi*(t-parms$phi*period)/period))# seasonality waves
  
  # baseline.txn.rate is the probability of transmission given contact per capita
  # (parms$dur.days1/length.step) is the duration of infectiousness of primary infection
  # q depends on transmission type (whether depends on population density or not)
  # density (q=0) vs frequency-dependent (q=1) transmission
  # c2 is the contact matrix 
  # beta is transmisibility per unit time
  transmission_unittime <-  parms$baseline.txn.rate/(parms$dur.days1/length.step)
  beta=transmission_unittime*parms$c2*parms$npi[t]
  
  beta_a_i <- seasonal.txn * beta/(sum(States)^parms$q)

  infectiousN <- I1 + parms$rho1*I2 + parms$rho2*I3 + parms$rho2*I4 + parms$introductions[t]
  lambda <- infectiousN %*% beta_a_i 
  lambda <- as.vector(lambda)
  ##########transmission dynamics##########################  
  
  dy <- matrix(NA, nrow=N.ages, ncol=ncol(States))
  colnames(dy) <- colnames(States)
  
  period.birth.rate <- log(parms$PerCapitaBirthsYear[t,]+1)/period #adjust annual birthrates to weekly scale
  
  #mu represents aging to the next class
  #um is death rate/emigration; adjust it to reproduce population growth
  
  Aging.Prop <- c(0,mu[1:(N.ages-1)])
  
  dy[,'M'] <- period.birth.rate*sum(States) - 
   parms$sigma3*lambda*M - 
    (omega+(mu+um))*M +
    Aging.Prop*c(0,M[1:(N.ages-1)]) 
  
  dy[,'S0'] <- omega*M -
    lambda*S0 - 
    (mu + um)*S0 + 
    Aging.Prop*c(0,S0[1:(N.ages-1)])
  
  dy[,'I1'] <-   lambda*S0 + parms$sigma3*lambda*M - 
    (gamma1 + mu + um)*I1 + 
    Aging.Prop*c(0,I1[1:(N.ages-1)]) 
  
  dy[,'R1'] <-  gamma1*I1 -  
    (waning1 + mu + um)*R1 + 
    Aging.Prop*c(0,R1[1:(N.ages-1)]) 
  
  dy[,'S1'] <-waning1*R1 - 
    parms$sigma1*lambda*S1 - 
    (mu+um)*S1 + 
    Aging.Prop*c(0,S1[1:(N.ages-1)]) 
  
  dy[,'I2'] <- parms$sigma1*lambda*S1 - 
    gamma2*I2-(mu + um)*I2 + 
    Aging.Prop*c(0,I2[1:(N.ages-1)]) 
  
  dy[,'R2'] <-  gamma2*I2 -  
    (waning2 + mu + um)*R2 + 
    Aging.Prop*c(0,R2[1:(N.ages-1)]) 
  
  dy[,'S2'] <- waning2*R2  -
    parms$sigma2*lambda*S2 -
    (mu+um)*S2 + 
    Aging.Prop*c(0,S2[1:(N.ages-1)]) 
  
  dy[,'I3'] <- parms$sigma2*lambda*S2 -
    (gamma3 + mu+um)*I3 +  
    Aging.Prop*c(0,I3[1:(N.ages-1)]) 
  
  dy[,'R3'] <-  gamma3*I3 -  
    (waning3 + mu + um)*R3 + 
    Aging.Prop*c(0,R3[1:(N.ages-1)]) 
  
  dy[,'S3'] <- waning3*R3 + waning4*R4  - 
    parms$sigma3*lambda*S3 -
    (mu + um)*S3 + 
    Aging.Prop*c(0,S3[1:(N.ages-1)]) 
  
  dy[,'I4'] <- parms$sigma3*lambda*S3 - 
    gamma4*I4 - 
    (mu + um)*I4 + 
    Aging.Prop*c(0,I4[1:(N.ages-1)]) 
  
  dy[,'R4'] <-  gamma4*I4 -  
    (waning4 + mu + um)*R4 + 
    Aging.Prop*c(0,R4[1:(N.ages-1)]) 
  
  derivs <- as.vector(dy)
  
  res <- list(derivs)
  
  return(res)
}
