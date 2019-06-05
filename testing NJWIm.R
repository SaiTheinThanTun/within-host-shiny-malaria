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
    #lst <- ((lst<=0)*0)+((lst>0)*lst) #test
    biglst[i,]<-lst
    i<-i+1   
  }
  
  data.frame(time=seq(1,runtime),log10=log10(apply(biglst,1,countrings)), normal=apply(biglst,1,countrings))
  #output the log of total observable parasites
  
}

#testing
outd <- NJWIm(10000,48,28,7,8,0,0,0,0,2400,72,0.693,54,.2,15,4)
head(outd)
min(outd[,2]) #-3.99315
min(outd[,3])
