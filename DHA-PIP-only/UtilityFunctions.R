# Ref. White NJ, Chapman D, Watt G (1992) The effects of multiplication and synchronicity 
# on the vascular distribution of parasites in falciparum malaria. Trans R Soc Trop Med Hyg 86:590-597.


h<-1 #what is h, used in function 5, drugaction
initconc<-1 #initial drug concentration

#1. age distribution of parasites####
InitAgeDistribution <- function(initN,lifecycle,mu,sigma)
{
  #initN: initial number of total parasites (across all stages)
  #lifecycle: number of hours parasites survive at birth. also serves as the length of storage at each stage
  #mu and sigma: parameters for standard density distribution 
  #centered at mu with SD sigma
  distr<-dnorm(1:lifecycle,mu,sigma)
  tot<-sum(distr[1:lifecycle])
  distr*initN/tot  
}
#2. aging####
ShiftOneHour<-function(lst,pmf)
{
  #lst is the density distribution of all parasites at a single timepoint
  #i.e. initial distribution is created by InitAgeDistribution function
  #parasites are assumed to live up to 48 hours
  #after which they multiply by pmf (parasite multiplication factor) and die (therefore shift in the array)
  n<-length(lst) #n is 48
  c(lst[n]*pmf,lst[1:(n-1)])
}
#3. for comparing with the observed####
countrings<-function(lst){
  MaxAge<-26   #maximum age that can be observed
  sum(lst[1:MaxAge]) #only up to hour 26
}


#4. drug concentration####
drugconcentration<-function(t,initconc,drugloss,halflife){
  #t timestep
  #initconc: initial concentration
  #drugloss: rate of drug loss?
  #halflife: time at which only 50% of original concentration will remain
  firstDose <- initconc*exp(-drugloss*(t/halflife))
  secondDose <- initconc*exp(-drugloss*(t/halflife))*1.24
  thirdDose <- initconc*exp(-drugloss*(t/halflife))*1.24^2
  
  (firstDose) + ((t>=24)*secondDose) + ((t>=48)*thirdDose)
}

#4.1. DHA or piperaquine concentration####
# source("DHApip.R")
# DHAconcentration <- DHAconc()
# PIPconcentration <- PIPconc()
#sourcing the code takes a long time, therefore 
#concentration data was saved in advance and loaded here
DHAconcentration <- readRDS("data/DHAconc.RDS")
PIPconcentration <- readRDS("data/PIPconc.RDS")
SpecificDrugConc <- function(t, drugname){
  if(drugname=="DHA"){
    out <- DHAconcentration[t,2]
  }
  if(drugname=="pip"){
    out <- PIPconcentration[t,2]
  }
  out
}

#4.2: concentration of third drug: PK function scaled from Piperaquine model####

library(Matrix)
library(pracma)

ThirdDrug <- function(dose=960, scaling){
  #for piperaquine
  #%solver for age and time model%
  #% uses analytical solution of pdes for constant source equation at each time step
  #clear
  #% global mu tau1 tau2 sigma beta n J BC alpha lam k r phi amp import
  
  # % Names of variables
  # % y(1) - dose compartment
  # % y(2) - central compartment
  # % y(3) - peripheral 1
  # % y(4) - peripheral 2
  # % y(5) - transit 1
  # % y(6) - transit 2
  # 
  # library(Matrix)
  # library(pracma)
  
  npeople = 1
  dosing = dose #scaling
  
  # %dimension of system
  n=6
  
  # %time steps
  da=1 #0.001
  maxa=(100*24) #10*24
  asteps=maxa/da
  nasteps=asteps
  
  # % dosing times
  dt=da
  mda1=(24/dt)
  mda2=(48/dt)
  mda3=(72/dt)
  
  
  # %%PK piperaquine parameters
  ACL = 0.75 
  AV = 1 
  M_WE = 54
  
  # THETA=[55.4;2910;310;4910;105;30900;2.11;1;1.06;0.236;0.575;5.51];
  # % 1.CL/F L/h
  # % 2.V2/F  L
  # % 3.Q1/F L/h
  # % 4.V3/F  L
  # % 5.Q2/F L/h
  # % 6.V4/F  L
  # % 7. MT
  # % 8. F
  # % 9. Scale
  # % 10.F dos occ
  # % 11.ET50
  # % 12.Hill
  
  Theta_1_CL_F = 55.4*scaling
  Theta_2_V2_F = 2910*scaling
  Theta_3_Q1_F = 310*scaling
  Theta_4_V3_F = 4910*scaling
  Theta_5_Q2_F = 105*scaling
  Theta_6_V4_F = 30900*scaling
  Theta_7_MT = 2.11*scaling
  #scaling to be done for above 7 parameters
  
  Theta_8_F = 1
  Theta_9_Scale = 1.06
  Theta_10_F_dos = 0.236
  Theta_11_ET50 = 0.575
  Theta_12_Hill = 5.51
  
  
  AGE=25
  WT=54
  EM50 = Theta_11_ET50
  HILL = Theta_12_Hill
  MF = ((AGE)^HILL)/(((EM50)^HILL)+((AGE)^HILL))
  
  OCC=1
  F1D = Theta_10_F_dos
  F1COVD = (1 + F1D*(OCC-1))
  
  
  TVCL = Theta_1_CL_F*MF*(WT/M_WE)^ACL
  TVV2 = Theta_2_V2_F*(WT/M_WE)^AV
  TVQ1 = Theta_3_Q1_F*(WT/M_WE)^ACL
  TVV3 = Theta_4_V3_F*(WT/M_WE)^AV
  TVQ2 = Theta_5_Q2_F*(WT/M_WE)^ACL
  TVV4 = Theta_6_V4_F*(WT/M_WE)^AV
  TVMT = Theta_7_MT
  TVF1 = Theta_8_F*F1COVD
  
  
  # % initialize solution matrix - rows are people , columns are timesteps
  # vv=zeros(npeople,asteps+1);
  #vv = matrix(0, npeople, asteps+1)
  #tic
  # % run system for each person
  # for(jj in 1:npeople){
  
  # % assign parameters to each person - mean estimate times randomly
  # % distributed number with mean 0 and variance omega
  CL = TVCL
  V2 = TVV2
  Q1 = TVQ1
  V3 = TVV3
  Q2 = TVQ2
  V4 = TVV4
  MT = TVMT
  F10 = TVF1
  F1 = F10
  
  NN = 2
  KTR = (NN+1)/MT
  
  K15 = KTR
  K56 = KTR
  K62 = KTR
  K23 = Q1/V2
  K32 = Q1/V3
  K24 = Q2/V2
  K42 = Q2/V4
  K20 = CL/V2
  
  S2 = V2*1000
  VD = V2+V3+V4
  
  
  # %set up V = solution matrix - rows is number of compartments or variable,
  # %columsn is number timesteps
  v=matrix(0,n,asteps+1)
  
  # %initialize jacobian like matrix
  J = matrix(0, n, n)
  
  for(m in 2:(asteps+1)){
    # %  state transition matrix
    diag(J) <- c(-K15,(-K23-K20-K24),-K32,-K42,-K56,-K62)
    J[2,3] <- K32
    J[2,4] <- K42
    J[2,6] <- K62
    J[3,2] <- K23
    J[4,2] <- K24
    J[5,1] <- K15
    J[6,5] <- K56
    
    
    # % calc new V
    ej=expm(J*dt)
    v[,m]= ej %*% v[,m-1]
    
    
    # % assign different dosing absorption rates per dose.
    
    if (m==mda1){
      v[1,m]=dosing*F1*0.577+v[1,m]
    }
    
    if (m==mda2){
      v[1,m]=dosing*1.24*F1*0.577+v[1,m]
    }
    
    if (m==mda3){
      v[1,m]=dosing*(1.24^2)*F1*0.577+v[1,m]
    }
  }
  vv2 <- v[2,]
  data.frame(time=seq(0,maxa),normal=(vv2))
}



drugconcentration_obsolete<-function(t,initconc,drugloss,halflife){
  #t timestep
  #initconc: initial concentration
  #drugloss: rate of drug loss?
  #halflife: time at which only 50% of original concentration will remain
  initconc*exp(-drugloss*(t/halflife))
}
#so far only one dose

#5. drug action: killing the parasites depending on the drug concentration and ce50####
drugaction<-function(t,killrate,drugconcentration,ce50,h){
  killrate/(1+(ce50/drugconcentration)^h) #ce50 is the tipping point of the sigmoidal curve, h is curvature of the sigmoidal curve
}
#overall killrate would be the sum of 2: artemisinin + partner drug


##6. NJW model main function + Immunity####
K<-function(a,t,delay,tpar,parasites,k0){
  k0*(parasites/tpar)+a*(t/(delay))
}
#adaptive immunity, tpar is the threshold for adaptive immunity, which grows linearly
#a innate immune response (lower than the adaptive immuntiy) which is delayed
#innate and 

#7. simulate parasite age distribution over time after taking a single dose of artemisinin####
#MIC detection is implemented here (not in NJWIm2, where 2 drugs are acting together) as it is separate for each drug#
NJWIm<-function(initn,lc,mu,sig,pmf,k0,a,tpar,delay,runtime,initconc,drugloss,halflife,killrate,ce50,h){
  biglst<-matrix(0,nrow=runtime,ncol=lc) #store parasite age distribution
  druglst<-matrix(0,nrow=runtime,ncol=lc) #store drug effect or conc? not used in this function!!
  
  lst <- InitAgeDistribution(initn,lc,mu,sig)
  biglst[1,]<-lst
  
  #to identify where is MIC
  growing <- rep(NA,runtime)
  growing[1] <- FALSE
  
  i=2
  while(i <= runtime)
  {
    priorLst <- sum(lst)
    #sai: not required anymore
    #if(i==0){parasites=initn}
    #else{parasites<-sum(biglst[i-1,])}
    
    #if timestep is = dose 1, add this much concentration
    #if dose 2, then add to the previous concentration
    #first dose is at time 0;
    
    #output should be like figure 6
    #integral of the drug concentration curve (this could pop up in the graph)
    #ACT failure is the result of the resistance to partner drug (if integral of parasites under regiem with partner drug is high)
    
    #user may change the killrates, halflife and dosages
    #default: everything sensitive
    #user: with resistance
    #toggle to do that
    
    #k<-K(a,i,delay,tpar,parasites,k0)
    drugconc<-drugconcentration(i,initconc,drugloss,halflife)
    drugeffect<-drugaction(i,killrate,drugconc,ce50,h)
    #lst<-ShiftOneHour(lst,pmf)*exp(-drugeffect)
    lst<-((lst<1)*0)+((lst>=1)*ShiftOneHour(lst,pmf)*exp(-drugeffect))
    #lst <- ((lst<=0)*0)+((lst>0)*lst)
    afterLst <- sum(lst)
    biglst[i,]<-lst
    
    #TODO
    if((afterLst-priorLst)<0){
      growing[i] <- FALSE
    } else {growing[i] <- TRUE}
    
    i<-i+1   
  }
  
  data.frame(time=seq(1,runtime),log10=log10(apply(biglst,1,countrings)), MIC=growing)
  #output the log of total observable parasites
  
}



#8. drug concentration over time (drug effect on the next function drugeff)####
drugf<-function(runtime,initconc,drugloss,halflife,killrate,ce50,h){
  druglst<-matrix(0,nrow=runtime,ncol=2)
  
  drugc <- drugconcentration(1,initconc,drugloss,halflife)
  druglst[1,1]<-drugc
  
  druge<-drugaction(1,killrate,drugc,ce50,h)
  druglst[1,2]<-druge
  
  i=2
  while(i <= runtime)
  {
    drugconc<-drugconcentration(i,initconc,drugloss,halflife)
    drugeffect<-drugaction(i,killrate,drugconc,ce50,h)
    druglst[i,1]<-drugconc
    druglst[i,2]<-drugeffect
    i<-i+1   
  }
  
  data.frame(time=seq(1,runtime),log10=log10(druglst[,1])) #<- only difference is here (column subset:1)
  
}



#9. drug effect over time (drug concentration is on the previous function drugf)####
drugeff<-function(runtime,initconc,drugloss,halflife,killrate,ce50,h){
  druglst<-matrix(0,nrow=runtime,ncol=2)
  
  drugc <- drugconcentration(1,initconc,drugloss,halflife)
  druglst[1,1]<-drugc
  
  druge<-drugaction(1,killrate,drugc,ce50,h)
  druglst[1,2]<-druge
  
  i=2
  while(i <= runtime)
  {
    drugconc<-drugconcentration(i,initconc,drugloss,halflife)
    drugeffect<-drugaction(i,killrate,drugconc,ce50,h)
    druglst[i,1]<-drugconc
    druglst[i,2]<-drugeffect
    i<-i+1   
  }
  
  data.frame(time=seq(1,runtime),log10=druglst[,2]) #<- only difference is here (column subset:2)
  
}

#9.1 DHA and piperaquine drug effect over time (drug concentration is on the previous function drugf)####
SpecificDrugEffect<-function(runtime,drugname, killrate,ce50,h){
  druglst<-matrix(0,nrow=runtime,ncol=2)
  
  if(drugname!="thirdDrug"){
    drugc <- SpecificDrugConc(1, drugname) #drugconcentration(1,initconc,drugloss,halflife)
    druglst[1,1]<-drugc
  }
  
  if(drugname=="thirdDrug"){}
  
  druge<-drugaction(1,killrate,drugc,ce50,h)
  druglst[1,2]<-druge
  
  i=2
  while(i <= runtime)
  {
    if(drugname!="thirdDrug"){drugconc<-SpecificDrugConc(i, drugname)} #drugconcentration(i,initconc,drugloss,halflife)
    if(drugname=="thirdDrug"){}
    drugeffect<-drugaction(i,killrate,drugconc,ce50,h)
    druglst[i,1]<-drugconc
    druglst[i,2]<-drugeffect
    i<-i+1   
  }
  
  data.frame(time=seq(1,runtime),log10=druglst[,2]) #<- only difference is here (column subset:2)
  
}

#10. takes 2 drugs: simulate parasite age distribution over time after taking a single dose of artemisinin####
NJWIm_2<-function(initn,lc,mu,sig,pmf,k0,a,tpar,delay,runtime,initconc,drugloss,halflife,killrate,ce50,h, initconc_2,drugloss_2,halflife_2,killrate_2,ce50_2, h_2){
  biglst<-matrix(0,nrow=runtime,ncol=lc) #store parasite age distribution
  druglst<-matrix(0,nrow=runtime,ncol=lc) #store drug effect or conc? not used in this function!!
  
  lst <- InitAgeDistribution(initn,lc,mu,sig)
  biglst[1,]<-lst
  i=2
  while(i <= runtime)
  {
    drugconc_1<-drugconcentration(i,initconc,drugloss,halflife)
    drugeffect_1<-drugaction(i,killrate,drugconc_1,ce50,h)
    
    drugconc_2<-drugconcentration(i,initconc_2,drugloss_2,halflife_2)
    drugeffect_2<-drugaction(i,killrate_2,drugconc_2,ce50_2, h_2)
    #lst <-ShiftOneHour(lst,pmf)*exp(-(drugeffect_1+drugeffect_2))
    lst<-((lst<1)*0)+((lst>=1)*ShiftOneHour(lst,pmf)*exp(-(drugeffect_1+drugeffect_2)))
    #lst<-((lst<=0)*0)+((lst>0)*ShiftOneHour(lst,pmf)*exp(-(drugeffect_1+drugeffect_2)))
    #lst <- ((lst<=0)*0)+((lst>0)*lst)
    #tmp <- c(1,1,2,0,0,3,-10,-2, -5,-9)
    #((tmp<=0)*0)+((tmp>0)*tmp)
    
    biglst[i,]<-lst
    i<-i+1   
  }
  
  data.frame(time=seq(1,runtime),log10=log10(apply(biglst,1,countrings)), normal=apply(biglst,1,countrings))
  #output the log of total observable parasites
}

#11. takes 3 drugs: simulate parasite age distribution over time after taking triple doses of respective drugs####
NJWIm_3<-function(initn,lc,mu,sig,pmf,k0,a,tpar,delay,runtime,initconc,drugloss,halflife,killrate,ce50,h, 
                  initconc_2,drugloss_2,halflife_2,killrate_2,ce50_2, h_2,
                  initconc_3,drugloss_3,halflife_3,killrate_3,ce50_3, h_3){
  biglst<-matrix(0,nrow=runtime,ncol=lc) #store parasite age distribution
  druglst<-matrix(0,nrow=runtime,ncol=lc) #store drug effect or conc? not used in this function!!
  
  lst <- InitAgeDistribution(initn,lc,mu,sig)
  biglst[1,]<-lst
  i=2
  while(i <= runtime)
  {
    drugconc_1<-drugconcentration(i,initconc,drugloss,halflife)
    drugeffect_1<-drugaction(i,killrate,drugconc_1,ce50,h)
    
    drugconc_2<-drugconcentration(i,initconc_2,drugloss_2,halflife_2)
    drugeffect_2<-drugaction(i,killrate_2,drugconc_2,ce50_2, h_2)
    
    drugconc_3<-drugconcentration(i,initconc_3,drugloss_3,halflife_3)
    drugeffect_3<-drugaction(i,killrate_3,drugconc_3,ce50_3, h_3)
    #lst <-ShiftOneHour(lst,pmf)*exp(-(drugeffect_1+drugeffect_2))
    lst<-((lst<1)*0)+((lst>=1)*ShiftOneHour(lst,pmf)*exp(-(drugeffect_1+drugeffect_2+drugeffect_3)))
    #lst<-((lst<=0)*0)+((lst>0)*ShiftOneHour(lst,pmf)*exp(-(drugeffect_1+drugeffect_2)))
    #lst <- ((lst<=0)*0)+((lst>0)*lst)
    #tmp <- c(1,1,2,0,0,3,-10,-2, -5,-9)
    #((tmp<=0)*0)+((tmp>0)*tmp)
    
    biglst[i,]<-lst
    i<-i+1   
  }
  
  data.frame(time=seq(1,runtime),log10=log10(apply(biglst,1,countrings)), normal=apply(biglst,1,countrings))
  #output the log of total observable parasites
}

#12. finding MIC ####
#TODO
TrueMIC <- function(MICvector){
  logicalSwitch <- FALSE
  trueMIC <- 0
  for(i in 1:length(MICvector)){
    if(!isTRUE(MICvector[i])){
      logicalSwitch <- FALSE
      trueMIC <- 0
    }
    if(isTRUE(MICvector[i]) && !logicalSwitch){
      logicalSwitch <- TRUE
      trueMIC <- i
    }
  }
  trueMIC
}

# the reason MIC is reached but the total observable parasites keeps falling is because
# 1. the way the drug effec is modelled (affecting all ages of parasites)
#     is different from the way the paraistes are multiplied (only those reaching the oldest age)
# 2. MIC is calulated based on all parasites. plot is only on observable parasites

#13. NJWIm_DHApip() ####
#copied from 10. takes 2 drugs: simulate parasite age distribution over time after taking a single dose of artemisinin
NJWIm_DHApip<-function(initn,lc,mu,sig,pmf,k0,a,tpar,delay,runtime,killrate,ce50,h, killrate_2,ce50_2, h_2){
  biglst<-matrix(0,nrow=runtime,ncol=lc) #store parasite age distribution
  druglst<-matrix(0,nrow=runtime,ncol=lc) #store drug effect or conc? not used in this function!!
  #gametocytes<-matrix(0,nrow=runtime,ncol=lc) #store gametocytes counts
  gametocytes <- NA
  MatureGametocytes <- NA
  
  lst <- InitAgeDistribution(initn,lc,mu,sig)
  biglst[1,]<-lst
  gametocytes[1] <- sum(lst)*0.001*pmf #kgam*pmf
  #MatureGametocytes[1] <- 
  i=2
  while(i <= runtime)
  {
    drugconc_1<-SpecificDrugConc(i,drugname="DHA")
    drugeffect_1<-drugaction(i,killrate,drugconc_1,ce50,h)
    
    drugconc_2<-SpecificDrugConc(i,drugname="pip")
    drugeffect_2<-drugaction(i,killrate_2,drugconc_2,ce50_2, h_2)
    lst<-((lst<1)*0)+((lst>=1)*ShiftOneHour(lst,pmf)*exp(-(drugeffect_1+drugeffect_2)))

    biglst[i,]<-lst
    if((i%%2)==0){
      gametocytes[i] <- gametocytes[i-1]+(sum(lst)*0.001*pmf) #kgam*pmf
      gametocytes[i] <- gametocytes[i]*exp(-1/9) #exp(-kg)
    } else{
      gametocytes[i] <- gametocytes[i-1]*exp(-1/9) #exp(-kg)
    }
    
    i<-i+1   
  }
  
  data.frame(time=seq(1,runtime),log10=log10(apply(biglst,1,countrings)), normal=apply(biglst,1,countrings), gam=log10(gametocytes)) #apply(gametocytes,1,sum)))
  #output the log of total observable parasites
}


#14. NJWIm_DHApip_C() ####
#copied from 13
NJWIm_DHApip_C<-function(initn,lc,mu,sig,pmf,k0,a,tpar,delay,runtime,killrate,ce50,h, killrate_2,ce50_2, h_2,
                         initconc_3,drugloss_3,halflife_3,killrate_3,ce50_3, h_3){
  biglst<-matrix(0,nrow=runtime,ncol=lc) #store parasite age distribution
  druglst<-matrix(0,nrow=runtime,ncol=lc) #store drug effect or conc? not used in this function!!
  
  lst <- InitAgeDistribution(initn,lc,mu,sig)
  biglst[1,]<-lst
  i=2
  while(i <= runtime)
  {
    drugconc_1<-SpecificDrugConc(i,drugname="DHA")
    drugeffect_1<-drugaction(i,killrate,drugconc_1,ce50,h)
    
    drugconc_2<-SpecificDrugConc(i,drugname="pip")
    drugeffect_2<-drugaction(i,killrate_2,drugconc_2,ce50_2, h_2)
    
    drugconc_3<-drugconcentration(i,initconc_3,drugloss_3,halflife_3)
    drugeffect_3<-drugaction(i,killrate_3,drugconc_3,ce50_3, h_3)
    
    lst<-((lst<1)*0)+((lst>=1)*ShiftOneHour(lst,pmf)*exp(-(drugeffect_1+drugeffect_2+drugeffect_3)))
    
    biglst[i,]<-lst
    i<-i+1   
  }
  
  data.frame(time=seq(1,runtime),log10=log10(apply(biglst,1,countrings)), normal=apply(biglst,1,countrings))
  #output the log of total observable parasites
}

