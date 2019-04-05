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
# % run system for each person
# for jj = 1: npeople
for(jj in 1:npeople){

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
# J=[
#   -K13 0 0 0 0 0 0 0 0;
#   0 -K20 0 0 0 0 0 0 K92;
#   K13 0 -K34 0 0 0 0 0 0;
#   0 0 K34 -K45 0 0 0 0 0;
#   0 0 0 K45 -K56 0 0 0 0;
#   0 0 0 0 K56 -K67 0 0 0;
#   0 0 0 0 0 K67 -K78 0 0;
#   0 0 0 0 0 0 K78 -K89 0;
#   0 0 0 0 0 0 0 K89 -K92;
#   ];

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
vv[jj,]=v[2,]

xticktimes=seq(from=24,by=24,to=maxa)
xticklabels=1:10
v[2,ncol(v)]=NaN
# %     patch((0:da:maxa),v(2,:),'black','EdgeColor','black','FaceAlpha',0.2);
# %     set(gca,'FontSize',14,'XTick',xticktimes,'XTickLabel',1:1:days)
# %     hold on
}
#end

#toc

# %millitres of blood in adult
bl=4.5*1000
# %nanograms in miligrams
ng=1*10^6

mvv=vv[1,]
semilogy(seq(0,maxa, by=da), mvv*ng/bl, type='l', xlim=c(0,71), ylim=c(0.5,10000), xlab="Time in hours")

# figure(2)
# semilogy((0:da:maxa),mvv*ng/bl,'k','LineWidth',3)
# set(gca,'FontSize',14,'XTick',xticktimes,'XTickLabel',xticklabels)
# ylim([0.5 10000])
# xlim([0 72])
# xlabel('Time in days')
# hold on

# % figure(3)
# % semilogy((0:da:maxa),mvv,'k','LineWidth',3)
# % set(gca,'FontSize',14,'XTick',xticktimes,'XTickLabel',xticklabels)
# % ylim([0.5 100])
# % xlim([0 48])
# % xlabel('Time in days')


