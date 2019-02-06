# Ref. White NJ, Chapman D, Watt G (1992) The effects of multiplication and synchronicity 
# on the vascular distribution of parasites in falciparum malaria. Trans R Soc Trop Med Hyg 86:590-597.

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

initconc<-1
#4. drug concentration####
drugconcentration<-function(t,initconc,drugloss,halflife){
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

h<-1
##6. NJW model main function + Immunity####
K<-function(a,t,delay,tpar,parasites,k0){
  k0*(parasites/tpar)+a*(t/(delay))
}
#adaptive immunity, tpar is the threshold for adaptive immunity, which grows linearly
#a innate immune response (lower than the adaptive immuntiy) which is delayed
#innate and 

#7. simulate parasite age distribution over time after taking a single dose of artemisinin####
NJWIm<-function(initn,lc,mu,sig,pmf,k0,a,tpar,delay,runtime,initconc,drugloss,halflife,killrate,ce50,h){
  biglst<-matrix(0,nrow=runtime,ncol=lc) #store parasite age distribution
  druglst<-matrix(0,nrow=runtime,ncol=lc) #store drug effect or conc? not used in this function!!

  lst <- InitAgeDistribution(initn,lc,mu,sig)
  biglst[1,]<-lst
  i=2
  while(i <= runtime)
  {
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
    lst<-ShiftOneHour(lst,pmf)*exp(-drugeffect)
    lst <- ((lst<=0)*0)+((lst>0)*lst) #test
    biglst[i,]<-lst
    i<-i+1   
  }
  
  data.frame(time=seq(1,runtime),log10=log10(apply(biglst,1,countrings)) ) #,normal=apply(biglst,1,countrings))
  #output the log of total observable parasites

}

#testing
outd <- NJWIm(10000,48,28,7,8,0,0,0,0,2400,72,0.693,54,.2,15,4)
head(outd)
min(outd[,2]) #-3.99315
#min(outd[,3]) #0.0001015897

library(manipulate)

# NJW + Immunity
manipulate(plot(NJWIm(initn,48,mu,sig,pmf,k0,a,tpar,delay,2400,initconc,0.693,halflife,killrate,ce50,h),xlim=c(0,245),ylim=c(0,8),type=show),
           initn=slider(10,1000),
           mu=slider(1,48),
           sig=slider(1,48),
           pmf=slider(1,30),
           # k0=slider(0,1,step = 0.01),
           # a=slider(0.0001,0.01,step = 0.001),
           # tpar=slider(1000,100000,step = 10),
           # delay=slider(0,4,step=1),
           initconc=slider(25,150,step=5),
           halflife=slider(5,240,step=1),
           killrate=slider(0.001,1,step=0.05),
           ce50=slider(10,100,step=5),
           h=slider(1,10,step=1),
           show=picker("line"="l", "points"="p")
)

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
  
  data.frame(time=seq(1,runtime),log10=druglst[,1]) #<- only difference is here (column subset:1)
  
}


manipulate(plot(drugf(2400,initconc,0.693,halflife,killrate,ce50,h),xlim=c(0,245),ylim=c(0,80),type=show),
           initconc=slider(25,150,step=5),
           halflife=slider(5,240,step=1),
           killrate=slider(0.001,1,step=0.05),
           ce50=slider(10,100,step=5),
           h=slider(1,10,step=1),
           show=picker("line"="l", "points"="p")
)

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
manipulate(plot(drugeff(2400,initconc,0.693,halflife,killrate,ce50,h),xlim=c(0,125),ylim=c(0,100),type=show),
      initconc=slider(25,150,step=5),
      halflife=slider(5,240,step=1),
      killrate=slider(1,100,step=0.5),
      ce50=slider(10,100,step=5),
      h=slider(1,10,step=1),
      show=picker("line"="l", "points"="p")
)



# ### plot model output and data
# PlotDatNJW<-function(data,initn,mu,sig,pmf,k0,a,tpar,delay){
#   
#   tmp<-NJWIm(initn,48,mu,sig,pmf,k0,a,tpar,delay,2400)[seq(1,7)*24,]
#   NJWout<-data.frame(day=seq(1,7),log10=tmp$log10)
#   
#   plot(data[,c(1,3)],type="b",xlim=c(0,8),ylim=c(0,8),ylab="",xlab="")
#   par(new=TRUE)
#   plot(NJWout,type="l",xlim=c(0,8),ylim=c(0,8),col = "red",xlab="day",
#        ylab="log10(parasitaemia)")
#   
#   legend(6,7,c("data","model"),lty=c(1,1), lwd=c(1,1), col=c("black","red") )
#   
# }
# 
# data1<-data.frame(day=c(1,2,3,4,5,6,7),
#                   asex=c(273,31,12160,60,24300,386,16040),
#                   log=log10(c(273,31,12160,60,24300,386,16040)));
# 
# manipulate(  
#   PlotDatNJW(data1,initn,mu,sig,pmf,k0,a,tpar,delay),
#   initn=slider(10,1000),
#   mu=slider(1,48),
#   sig=slider(1,48),
#   pmf=slider(1,30),
#   k0=slider(0,0.1),
#   a=slider(0.0001,0.01),
#   delay=slider(1,4),
#   tpar=slider(1000,10000)
# )
# 
