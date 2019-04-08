#1 person version 20190405

# % Names of variables
# % COMP = (1)  ;(ABSORB, DEFDOSE)
# % COMP = (2)  ;(CENTRAL)
# % COMP = (3)  ;(TRANSIT 1)
# % COMP = (4)  ;(TRANSIT 2)
# % COMP = (5)  ;(TRANSIT 3)
# % COMP = (6)  ;(TRANSIT 4)
# % COMP = (7)  ;(TRANSIT 5)
# % COMP = (8)  ;(TRANSIT 6)
# % COMP = (9)  ;(TRANSIT 7)

library(Matrix)
library(pracma)

DHAconc <- function(){
  
npeople = 1
dosing = 40
# % flag = pregnant
# % para is initial parasitaemia
# % scale is the natural log of patient parasitaemia


# %dimension of system
n=9

# %time steps
# %size of time step
da=1 #.003
days=10
# %maximum age in hours for 10 days period
maxa=days*24
# %no. of timesteps
asteps=maxa/da
nasteps=asteps
# %no. of timesteps in a day
oneday=24*(1/da)

# % dosing times, 
# % calculating index for mda
dt=da
mda1=(24/dt)
mda2=(48/dt)
mda3=(72/dt)


# %%PK piperaquine parameters
ACL = 0.75
AV = 1
M_WE = 48.5

# THETA=[78.0;129;0.982;1;0.58;-0.375;0.278];
# % 1.CL/F L/h
# % 2.V2/F  L
# % 3. MT
# % 4. F1
# % 5. w
# % 6. F1 Flag
# % 7.F1 para
Theta_CL_F <- 78.0
Theta_V2_F <- 129
Theta_MT <- 0.982
Theta_F1 <- 1
Theta_w <- 0.58
Theta_F1_flag <- -0.375
Theta_F1_para <- 0.278

FLAG=1

# OMEGA=[1 ; 0.0162; 1; 0.088; 1; 0.23 ];
Omega_CL <- 1
Omega_V2 <- .0162
Omega_MT <-1
Omega_F1 <- .088
Omega_W <- 1
Omega_MTT_IOV <- 0.23

WT=48.5
PARA=3.98

# TVCL = THETA(1)*(WT/M_WE)^ACL;
# TVV2 = THETA(2)*(WT/M_WE)^AV;
# TVMT = THETA(3);
# TVF1 = THETA(4)*(1+THETA(6)*FLAG)*(1+THETA(7)*(PARA-3.98));
TVCL = Theta_CL_F*(WT/M_WE)^ACL
TVV2 = Theta_V2_F*(WT/M_WE)^AV
TVMT = Theta_MT
TVF1 = Theta_F1*(1+Theta_F1_flag*FLAG)*(1+Theta_F1_para*(PARA-3.98))


# % initialize solution matrix - rows are people , columns are timesteps
vv=matrix(0,npeople,asteps+1)

# tic
  
  # % assign parameters to each person - mean estimate times randomly
  # % distributed number with mean 0 and variance omega
  CL = TVCL
  V2 = TVV2
  MT = TVMT
  # %     F10 = TVF1;  % CHECK THIS
  F10 = 1
  F1 = F10
  
  NN = 7
  KTR = (NN+1)/MT
  
  K13 = KTR
  K34 = KTR
  K45 = KTR
  K56 = KTR
  K67 = KTR
  K78 = KTR
  K89 = KTR
  K92 = KTR
  K20 = CL/V2
  
  S2 = V2/1000
  VD = V2
  
  W= Theta_w #THETA(5);
  
  
  # %set up V = solution matrix - rows is number of compartments or variable,
  # %colums is number timesteps
  v=matrix(0,n,asteps+1)
  
  # %initialize jacobian like matrix
  J=matrix(0,n,n)
  
  # for m=2:asteps+1
  for(m in 2:(asteps+1))
  {
    
    # %  state transition matrix
    diag(J) <- c(-K13,-K20,-K34,-K45,-K56,-K67,-K78,-K89,-K92)
    J[2,9] <- K92
    J[3,1] <- K13
    J[4,3] <- K34
    J[5,4] <- K45
    J[6,5] <- K56
    J[7,6] <- K67
    J[8,7] <- K78
    J[9,8] <- K89
    
    # % calc new V
    ej=expm(J*dt)
    v[,m]= ej %*% v[,m-1]
    #v[,m]= as.matrix(ej %*% v[,m-1])
    
    # % assign different dosing absorption rates per dose.
    if (m==mda1){
      v[1,m]=dosing*F1+v[1,m]
    }
    
    if (m==mda2){
      v[1,m]=dosing*F1+v[1,m]
    }
    
    if (m==mda3){
      v[1,m]=dosing*F1+v[1,m]
    }
    
  }
  
  # % keep solution for compartment 2 for each individual
  vv[1,]=v[2,] #1 person version, matrix
  vv2 <- v[2,]
  vv2
}
  # xticktimes=seq(from=24,by=24,to=maxa)
  # xticklabels=1:10
  # v[2,ncol(v)]=NaN

#toc

# # millitres of blood in adult
# bl=4.5*1000
# # nanograms in miligrams
# ng=1*10^6
# 
# mvv=vv[1,]
# semilogy(seq(0,maxa, by=da), mvv*ng/bl, type='l', xlim=c(0,71), ylim=c(0.5,10000), xlab="Time in hours")
# 
# #
# semilogy(seq(0,maxa, by=da), mvv, type='l', xlim=c(0,71), ylim=c(0.5,10000), xlab="Time in hours")
# 
# 
# #checking out if the scale is correct
# vv[1,25]
# mvv[25]
# mvv[25]*ng/bl
# 
# plot(seq(0,maxa, by=da), mvv, type='l')


PIPconc <- function(){
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
  dosing = 960
  
  # %dimension of system
  n=6
  
  # %time steps
  da=1 #0.001
  maxa=10*24
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
  Theta_1_CL_F = 55.4
  Theta_2_V2_F = 2910
  Theta_3_Q1_F = 310
  Theta_4_V3_F = 4910
  Theta_5_Q2_F = 105
  Theta_6_V4_F = 30900
  Theta_7_MT = 2.11
  Theta_8_F = 1
  Theta_9_Scale = 1.06
  Theta_10_F_dos = 0.236
  Theta_11_ET50 = 0.575
  Theta_12_Hill = 5.51
  
  
  
  # OMEGA=[0.0752 ; 0.371; 0; 0.0558; 0.0541; 0.114 ; 0.135; 0.158; 0; 0.252; 0.252;0.252;0.195;0.195;0.195];
  # %1. CL
  # %2. V2
  # %3. Q1
  # %4. V3
  # %5. Q2
  # %6. V4
  # %7. MT
  # %8. F
  # 
  # Omega_1_CL = 0.0752
  # Omega_2_V2 = 0.371
  # Omega_3_Q1 = 0
  # Omega_4_V3 = 0.0558
  # Omega_5_Q2 = 0.0541
  # Omega_6_V4 = #NO USAGE OF OMEGA FOUND!!!
  
  
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
    # v=zeros(n,asteps+1);
    v=matrix(0,n,asteps+1)
    
    # %initialize jacobian like matrix
    # J=zeros(n,n);
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
      
      # J=[
      #   -K15 0 0 0 0 0;
      #   0 (-K23-K20-K24) K32 K42 0 K62 ;
      #   0 K23 -K32 0 0 0;
      #   0 K24 0 -K42 0 0;
      #   K15 0 0 0 -K56 0;
      #   0 0 0 0 K56 -K62;
      #   ];
      
      
      # % calc new V
      ej=expm(J*dt)
      v[,m]= ej %*% v[,m-1]
      
      # ej=expm(J(:,:)*dt);
      # v(:,m)=ej*v(:,m-1);
      
      # % assign different dosing absorption rates per dose.
      
      if (m==mda1){
        v[1,m]=dosing*F1*0.577+v[1,m]
        # v(1,m)=dosing*F1*0.577+v(1,m);
      }
      
      if (m==mda2){
        v[1,m]=dosing*1.24*F1*0.577+v[1,m]
        # v(1,m)=dosing*1.24*F1*0.577+v(1,m);
      }
      
      if (m==mda3){
        v[1,m]=dosing*(1.24^2)*F1*0.577+v[1,m]
        # v(1,m)=dosing*F1*(1.24^2)*0.577+v(1,m);
      }
      
    }
    
    # % keep solution for compartment 2 for each individual
    # vv(jj,:)=v(2,:);
    #vv[jj,]=v[2,]#matrix output
    vv2 <- v[2,]
    vv2
    # xticktimes=(24:24*5:maxa);
    # xticklabels=(1:100);
    # v(2,end)=NaN;
    # patch((0:maxa/asteps:maxa),v(2,:),'black','EdgeColor','black','FaceAlpha',0.2);
    # set(gca,'FontSize',14,'XTick',xticktimes,'XTickLabel',1:2:100)
    # hold on
  #}
  
  # toc
  
  #mvv=vv[1,] 
  
  
}
