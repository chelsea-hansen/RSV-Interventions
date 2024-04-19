interventions = function(birth_dose, #timing and coverage of nirsevimab birth doses
                         cover_n, #timing and coverage of nirsevimab catch-up doses
                         waningN, #duration of nirsevimab protection (in days)
                         RRIn, #relative risk of infection while protected by nirsevimab (default = 1)
                         RRHn,#relative risk of hospitalization while protected by nirsevimab
                         maternal_dose, #timing and coverage of infants born following maternal doses
                         waningV, #duration of maternal vaccination protection (in days)
                         RRIv, #relative risk of infection while protected by maternal vaccination (default = 1)
                         RRHv,#relative risk of hospitalization while protected by maternal vaccination 
                         cover_s, #timing and coverage of senior vaccination 
                         waningS, #duration of protection from vaccines (in days) for seniors
                         RRIs, #relative risk of infection after vaccination in seniors (default = 1)
                         RRHs #relative risk of hospitalization after vaccination in seniors
){
  
  
  parmset<-list(PerCapitaBirthsYear=birth,
                WidthAgeClassMonth=c(rep(2,times=6), 12,12*3, 60, 120, 240, 240, 240),
                DurationMatImmunityDays=112,
                baseline.txn.rate=baseline.txn.rate,
                RRHm = 0.7,
                b1=b1,
                phi=phi,
                recover1 = 365.25*.5,
                recover2 = 365.25*.5,
                recover3 = 365.25,
                recover4 = 365.25, 
                age_reporting=age_reporting,
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
                npi=npi,
                introductions = introductions,
                time.step='week',
                #intervention parameters
                birth_dose = birth_dose, 
                cover_n = cover_n,
                waningN = waningN,
                RRIn = RRIn,
                RRHn = RRHn,
                maternal_dose = maternal_dose, 
                waningV = waningV,
                RRIv = RRIv,
                RRHv = RRHv,
                cover_s = cover_s, 
                waningS = waningS, 
                RRIs = RRIs, 
                RRHs =RRHs
  )
  
  
  output <- ode(y=yinit.vector, t=times,
                method = "ode45",
                func=MSIRS_interventions, 
                parms=c(parmset))
  
  
  t0 <- tmax3-tmax0+1 #just save results after burn-in period
  al <- nrow(yinit)
  output <- tail(output,t0)
  St <- output[,-1]
  
  I1 <- St[,grep('I1', colnames(St))]
  I2 <- St[,grep('I2', colnames(St))]
  I3 <- St[,grep('I3', colnames(St))]
  I4 <- St[,grep('I4', colnames(St))]
  R1 <- St[,grep('R1', colnames(St))]
  R2 <- St[,grep('R2', colnames(St))]
  R3 <- St[,grep('R3', colnames(St))]
  R4 <- St[,grep('R4', colnames(St))]
  S0 <- St[,grep('S0', colnames(St))]
  Si <- St[,grep('Si', colnames(St))]
  S1 <- St[,grep('S1', colnames(St))]
  S2 <- St[,grep('S2', colnames(St))]
  S3 <- St[,grep('S3', colnames(St))]
  
  M <- St[,grep('M', colnames(St))][,1:13]
  Mn <- St[,grep('Mn', colnames(St))]
  Mv <- St[,grep('Mv', colnames(St))]
  N <- St[,grep('N', colnames(St))]

  
  Vs1 <- St[,grep('Vs1', colnames(St))]
  Vs2 <- St[,grep('Vs2', colnames(St))]
  
  contact2 = npi[tmax0:tmax3]
  intro2 = introductions[tmax0:tmax3]
  
  lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
  for (t in 1:t0){
    beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(parmset$yinit.matrix)^(1-parmset$q))*parmset$c2*contact2[t]
    lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+intro2[t])%*%beta)/sum(St[t,]))}
  
 hosp = c(rep(parmset$age_reporting[1],3),rep(parmset$age_reporting[2],3),rep(parmset$age_reporting[3],2),rep(parmset$age_reporting[4],4),parmset$age_reporting[5])
  
  H1=matrix(0,nrow=t0,ncol=al)#Number of hospitalizations by age
  for (i in 1:al){
    H1[,i]=
      hosp[i]*parmset$RRHm*parmset$sigma3*M[,i]*lambda1[,i]+
      hosp[i]*parmset$RRHn*parmset$RRIn*parmset$sigma3*Mn[,i]*lambda1[,i]+
      hosp[i]*parmset$RRHm*parmset$RRHv*parmset$RRIv*parmset$sigma3*Mv[,i]*lambda1[,i]+
      hosp[i]*parmset$RRHn*parmset$RRIn*N[,i]*lambda1[,i]+
      hosp[i]*S0[,i]*lambda1[,i]+
      hosp[i]*Si[,i]*lambda1[,i]+
      hosp[i]*parmset$sigma1*S1[,i]*lambda1[,i]+
      hosp[i]*parmset$sigma2*S2[,i]*lambda1[,i]+
      hosp[i]*parmset$sigma3*S3[,i]*lambda1[,i]+
      hosp[i]*parmset$RRIs*parmset$RRHs*parmset$sigma3*Vs1[,i]*lambda1[,i]+
      hosp[i]*parmset$RRIs*parmset$RRHs*parmset$sigma3*Vs2[,i]*lambda1[,i]}
  
  
  H2= cbind(rowSums(H1[,1:3]),rowSums(H1[,4:6]),rowSums(H1[,7:8]),rowSums(H1[,9:12]),H1[,13]) 
  results = data.frame(cbind(H2[(t0-(tmax3-tmax2)):t0,],rowSums(H2[(t0-(tmax3-tmax2)):t0,])))
  names(results)=c("<6m",">6m","1-4yr","5-59yrs","60+yrs","All ages")
  results$date = dates$date[tmax2:tmax3]
  
  return(results)
}






interventions_PI = function(birth_dose, #timing and coverage of nirsevimab birth doses
                             cover_n, #timing and coverage of nirsevimab catch-up doses
                             waningN, #duration of nirsevimab protection (in days)
                             RRIn, #relative risk of infection while protected by nirsevimab (default = 1)
                             RRHn,#relative risk of hospitalization while protected by nirsevimab
                             maternal_dose, #timing and coverage of infants born following maternal doses
                             waningV, #duration of maternal vaccination protection (in days)
                             RRIv, #relative risk of infection while protected by maternal vaccination (default = 1)
                             RRHv,#relative risk of hospitalization while protected by maternal vaccination 
                             cover_s, #timing and coverage of senior vaccination 
                             waningS, #duration of protection from vaccines (in days) for seniors
                             RRIs, #relative risk of infection after vaccination in seniors (default = 1)
                             RRHs #relative risk of hospitalization after vaccination in seniors
){
  results = data.frame(`<6m`=NA,`>6m`=NA,`1-4yr`=NA,`5-59yrs`=NA,`60+yrs`=NA,`All ages`=NA,date=NA,sample=NA)
  names(results)=c("<6m",">6m","1-4yr","5-59yrs","60+yrs","All ages","date","sample") 
  
  for(l in 1:rep_num){
    parmset<-list(PerCapitaBirthsYear=birth,
                  WidthAgeClassMonth=c(rep(2,times=6), 12,12*3, 60, 120, 240, 240, 240),
                  DurationMatImmunityDays=112,
                  baseline.txn.rate=new_parms[l,"beta"],
                  RRHm = 0.7,
                  b1=new_parms[l,"b1"],
                  phi=new_parms[l,"phi"],
                  recover1 = 365.25*.5,
                  recover2 = 365.25*.5,
                  recover3 = 365.25,
                  recover4 = 365.25, 
                  h1=new_parms[l,"h1"],
                  h2=new_parms[l,"h2"],
                  h3=new_parms[l,"h3"],
                  h4=new_parms[l,"h4"],
                  h5=new_parms[l,"h5"],
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
                  npi=npi,
                  introductions = introductions,
                  time.step='week',
                  #intervention parameters
                  birth_dose = birth_dose, 
                  cover_n = cover_n,
                  waningN = waningN,
                  RRIn = RRIn,
                  RRHn = RRHn,
                  maternal_dose = maternal_dose, 
                  waningV = waningV,
                  RRIv = RRIv,
                  RRHv = RRHv,
                  cover_s = cover_s, 
                  waningS = waningS, 
                  RRIs = RRIs, 
                  RRHs =RRHs
                  )

    
    output <- ode(y=yinit.vector, t=times,
                  method = "ode45",
                  func=MSIRS_interventions, 
                  parms=c(parmset))
    
    
    t0 <- tmax3-tmax0+1 #just save results after burn-in period
    al <- nrow(yinit)
    output <- tail(output,t0)
    St <- output[,-1]
    
    I1 <- St[,grep('I1', colnames(St))]
    I2 <- St[,grep('I2', colnames(St))]
    I3 <- St[,grep('I3', colnames(St))]
    I4 <- St[,grep('I4', colnames(St))]
    R1 <- St[,grep('R1', colnames(St))]
    R2 <- St[,grep('R2', colnames(St))]
    R3 <- St[,grep('R3', colnames(St))]
    R4 <- St[,grep('R4', colnames(St))]
    S0 <- St[,grep('S0', colnames(St))]
    Si <- St[,grep('Si', colnames(St))]
    S1 <- St[,grep('S1', colnames(St))]
    S2 <- St[,grep('S2', colnames(St))]
    S3 <- St[,grep('S3', colnames(St))]
    
    M <- St[,grep('M', colnames(St))][,1:13]
    Mn <- St[,grep('Mn', colnames(St))]
    Mv <- St[,grep('Mv', colnames(St))]
    N <- St[,grep('N', colnames(St))]
    
    Vs1 <- St[,grep('Vs1', colnames(St))]
    Vs2 <- St[,grep('Vs2', colnames(St))]
    
    contact2 = npi[tmax0:tmax3]
    intro2 = introductions[tmax0:tmax3]
    
    lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
    for (t in 1:t0){
      beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(parmset$yinit.matrix)^(1-parmset$q))*parmset$c2*contact2[t]
      lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+intro2[t])%*%beta)/sum(St[t,]))}
    
    
    hosp = c(rep(parmset$h1,3),rep(parmset$h2,3),rep(parmset$h3,2),rep(parmset$h4,4),parmset$h5)
    
    H1=matrix(0,nrow=t0,ncol=al)#Number of hospitalizations by age
    for (i in 1:al){
      H1[,i]=
        hosp[i]*parmset$RRHm*parmset$sigma3*M[,i]*lambda1[,i]+
        hosp[i]*parmset$RRHn*parmset$RRIn*parmset$sigma3*Mn[,i]*lambda1[,i]+
       hosp[i]*parmset$RRHm*parmset$RRHv*parmset$RRIv*parmset$sigma3*Mv[,i]*lambda1[,i]+
       hosp[i]*parmset$RRHn*parmset$RRIn*parmset$RRHm*N[,i]*lambda1[,i]+
       hosp[i]*S0[,i]*lambda1[,i]+
       hosp[i]*Si[,i]*lambda1[,i]+
      hosp[i]*parmset$sigma1*S1[,i]*lambda1[,i]+
       hosp[i]*parmset$sigma2*S2[,i]*lambda1[,i]+
       hosp[i]*parmset$sigma3*S3[,i]*lambda1[,i]+
       hosp[i]*parmset$RRIs*parmset$RRHs*parmset$sigma3*Vs1[,i]*lambda1[,i]+
       hosp[i]*parmset$RRIs*parmset$RRHs*parmset$sigma3*Vs2[,i]*lambda1[,i]}
    
    
    
    H2= cbind(rowSums(H1[,1:3]),rowSums(H1[,4:6]),rowSums(H1[,7:8]),rowSums(H1[,9:12]),H1[,13]) 
    new_results = data.frame(cbind(H2[(t0-(tmax3-tmax2)):t0,],rowSums(H2[(t0-(tmax3-tmax2)):t0,])))
    names(new_results)=c("<6m",">6m","1-4yr","5-59yrs","60+yrs","All ages")
    new_results$date = dates$date[tmax2:tmax3]
    new_results$sample = l
    
    results = rbind(results,new_results)
    
  }
  return(results)
}