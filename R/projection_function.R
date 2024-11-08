projection_function = function(RRHn1,
                               RRHn2,
                               RRHv1,
                               RRHv2,
                               RRHs,
                               RRIn,
                               RRIv,
                               RRIs,
                               waningN1,
                               waningN2,
                               waningV1,
                               waningV2,
                               waningS,
                               monoclonal_birth,
                               monoclonal_01,
                               monoclonal_23,
                               monoclonal_45,
                               monoclonal_67,
                               maternal_vax, 
                               senior_vax_75,
                               senior_vax_60_74,
                               parmset,
                               lhs_parms,
                               min_date,
                               max_date,
                               fit_times){

  parmset$RRHn1 = RRHn1
  parmset$RRHn2 = RRHn2
  parmset$RRHv1 = RRHv1
  parmset$RRHv2 = RRHv2
  parmset$RRHs = RRHs
  parmset$RRIn = RRIn
  parmset$RRIv = RRIv
  parmset$RRIs = RRIs
  parmset$waningN1 = waningN1
  parmset$waningN2 = waningN2
  parmset$waningV1 = waningV1
  parmset$waningV2 = waningV2
  parmset$waningS = waningS 
  parmset$monoclonal_birth = monoclonal_birth
  parmset$monoclonal_catchup_01 = monoclonal_01
  parmset$monoclonal_catchup_23 = monoclonal_23
  parmset$monoclonal_catchup_45 = monoclonal_45
  parmset$monoclonal_catchup_67 = monoclonal_67
  parmset$maternal_vax = maternal_vax
  parmset$senior_vax_75 = senior_vax_75
  parmset$senior_vax_60_74 = senior_vax_60_74
  

  newH=data.frame(date=NA,sample=NA,Age=NA,value=NA)
  dates = seq(from=as.Date(min_date),to=as.Date(max_date),by='week')
  for(l in 1:nrow(lhs_parms)){
    parmset$baseline.txn.rate=lhs_parms[l,"beta"]
    parmset$b1=lhs_parms[l,"b1"]
    parmset$phi=lhs_parms[l,"phi"]
    report_seniors75 <- lhs_parms[l,"RS75"]
    report_seniors60 <- lhs_parms[l,"RS60"]
    report_infants <- lhs_parms[l,"RI"]
    report_children <- lhs_parms[l,"RC"]
    report_adults <- lhs_parms[l,"RA"]
    npi1 <- lhs_parms[l,"npi1"]
    npi2 <- lhs_parms[l,"npi2"]
    npi3 <- lhs_parms[l,"npi3"]
   # npi4 <- lhs_parms[l,"npi4"]
    
    introductions = data.frame(intros=c(rep(parmset$seed,248),rep(0,44),rep(NA,24),rep(parmset$seed,250))) %>% 
      mutate(intros = na_interpolation(intros, method="linear"))
    parmset$introductions = introductions$intros
    
    npi = data.frame(npis=c(rep(1,248),rep(npi1,52),rep(NA,22),rep(npi2,13),rep(npi3,18),rep(NA,13),rep(1,188)))%>% 
      mutate(npis= na_interpolation(npis, method="linear"))
    parmset$npi = npi$npis
    
    
    
    
    output <- ode(y=yinit.vector, times=fit_times,method = "ode45",
                  func=MSIRS_immunization_dynamics, 
                  parms=parmset)
    
    t0=length(fit_times)-104
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
    M0<- St[,grep('M0', colnames(St))]
    Si<- St[,grep('Si', colnames(St))]
    Mn1<- St[,grep('Mn1', colnames(St))]
    Mn2<- St[,grep('Mn2', colnames(St))]
    Mv1<- St[,grep('Mv1', colnames(St))]
    Mv2<- St[,grep('Mv2', colnames(St))]
    N1<- St[,grep('N1', colnames(St))]
    N2<- St[,grep('N2', colnames(St))]
    Vs1<- St[,grep('Vs1', colnames(St))]
    Vs2<- St[,grep('Vs2', colnames(St))]
    
    beta <-  parmset$baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2
    
    
    lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
    for (t in 1:t0){
      lambda1[t,] <- as.vector((1+parmset$b1*cos(2*pi*(t-parmset$phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+parmset$seed)%*%beta)/sum(St[t,]))}
    
   
    hosp1 <- c(report_infants,report_infants*.59,report_infants*.33,report_infants*.2,report_infants*.15,report_infants*.15,rep(report_children,2),rep(.001,6))
    hosp2 <- hosp1 * 0.4
    hosp3 <- c(rep(0.001, 8), rep(report_adults, 4),report_seniors60, report_seniors75)
    
  H1=matrix(0,nrow=t0,ncol=al)#Number of hospitalizations by age
    for (i in 1:al){
      H1[,i]=
        hosp1[i]*parmset$sigmaM*M0[,i]*lambda1[,i]+
        hosp1[i]*parmset$RRHv1*parmset$RRIv*parmset$sigmaM*Mv1[,i]*lambda1[,i]+
        hosp1[i]*parmset$RRHv2*parmset$RRIv*Mv2[,i]*lambda1[,i]+
        hosp1[i]*parmset$RRHn1*parmset$RRIn*parmset$sigmaM*Mn1[,i]*lambda1[,i]+
        hosp1[i]*parmset$RRHn2*parmset$RRIn*Mn2[,i]*lambda1[,i]+
        hosp1[i]*parmset$RRHn1*parmset$RRIn*N1[,i]*lambda1[,i]+
        hosp1[i]*parmset$RRHn2*parmset$RRIn*N2[,i]*lambda1[,i]+
        hosp1[i]*S0[,i]*lambda1[,i]+
        hosp1[i]*Si[,i]*lambda1[,i]+
        hosp2[i]*parmset$sigma1*S1[,i]*lambda1[,i]+
        hosp3[i]*parmset$sigma2*S2[,i]*lambda1[,i]+
        hosp3[i]*parmset$sigma3*S3[,i]*lambda1[,i]+
        hosp3[i]*parmset$RRIs*parmset$RRHs*parmset$sigma3*Vs1[,i]*lambda1[,i]+
        hosp3[i]*parmset$RRIs*parmset$RRHs*parmset$sigma3*Vs2[,i]*lambda1[,i]}
    
    
    H2 = cbind(rowSums(H1[,1:3]),
               rowSums(H1[,4:6]),
               rowSums(H1[,7:8]),
               rowSums(H1[,9:12]),
               H1[,13], H1[,14])
    
    
    
    
    # H3 = t(t(H2) * hosp_prop)
 
  
    
    H = data.frame(H2)
    names(H)=c("<6m","6-11m","1-4yrs","5-59yrs","60-74yrs","75+yrs")
    H$All = rowSums(H2)
    H$date = dates
    
    H = H %>% 
      mutate(sample=l) %>% 
      pivot_longer(cols=c("<6m":"All"),names_to="Age",values_to="value")
    
    
    newH=rbind(newH,H) %>% filter(!is.na(date))
    
  }
  return(newH)
}


