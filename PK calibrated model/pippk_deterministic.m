%solver for age and time model%
% uses analytical solution of pdes for constant source equation at each time step
clear
close all
% global mu tau1 tau2 sigma beta n J BC alpha lam k r phi amp import

% Names of variables
% y(1) - dose compartment
% y(2) - central compartment
% y(3) - peripheral 1
% y(4) - peripheral 2
% y(5) - transit 1
% y(6) - transit 2

npeople = 1;
dosing = 960;

%dimension of system
n=6;

%time steps
%da=0.001;
da=1;
maxa=10*24;
asteps=maxa/da;
nasteps=asteps;

% dosing times
dt=da;
mda1=int32(24/dt);
mda2=int32(48/dt);
mda3=int32(72/dt);


%%PK piperaquine parameters
ACL = 0.75;
AV = 1;
M_WE = 54;

THETA=[55.4;2910;310;4910;105;30900;2.11;1;1.06;0.236;0.575;5.51];
% 1.CL/F L/h
% 2.V2/F  L
% 3.Q1/F L/h
% 4.V3/F  L
% 5.Q2/F L/h
% 6.V4/F  L
% 7. MT
% 8. F
% 9. Scale
% 10.F dos occ
% 11.ET50
% 12.Hill

OMEGA=[0.0752 ; 0.371; 0; 0.0558; 0.0541; 0.114 ; 0.135; 0.158; 0; 0.252; 0.252;0.252;0.195;0.195;0.195];
%CL
%V2
%Q1
%V3
%Q2
%V4
%MT
%F
AGE=25;
WT=54;
EM50 = THETA(11);
HILL = THETA(12);
MF = ((AGE)^HILL)/(((EM50)^HILL)+((AGE)^HILL));

OCC=1;
F1D = THETA(10);
F1COVD = (1 + F1D*(OCC-1));


TVCL = THETA(1)*MF*(WT/M_WE)^ACL;
TVV2 = THETA(2)*(WT/M_WE)^AV;
TVQ1 = THETA(3)*(WT/M_WE)^ACL;
TVV3 = THETA(4)*(WT/M_WE)^AV;
TVQ2 = THETA(5)*(WT/M_WE)^ACL;
TVV4 = THETA(6)*(WT/M_WE)^AV;
TVMT = THETA(7);
TVF1 = THETA(8)*F1COVD;


% initialize solution matrix - rows are people , columns are timesteps
vv=zeros(npeople,asteps+1);

tic
% run system for each person
for jj = 1: npeople
    
    % assign parameters to each person - mean estimate times randomly
    % distributed number with mean 0 and variance omega
    CL = TVCL;
    V2 = TVV2;
    Q1 = TVQ1;
    V3 = TVV3;
    Q2 = TVQ2;
    V4 = TVV4;
    MT = TVMT;
    F10 = TVF1;
    F1 = F10;
    
    NN = 2;
    KTR = (NN+1)/MT;
    
    K15 = KTR;
    K56 = KTR;
    K62 = KTR;
    K23 = Q1/V2;
    K32 = Q1/V3;
    K24 = Q2/V2;
    K42 = Q2/V4;
    K20 = CL/V2;
    
    S2 = V2*1000;
    VD = V2+V3+V4;
    
    
    %set up V = solution matrix - rows is number of compartments or variable,
    %columsn is number timesteps
    v=zeros(n,asteps+1);
    
    %initialize jacobian like matrix
    J=zeros(n,n);
    
    for m=2:asteps+1
        
        %  state transition matrix
        J=[
            -K15 0 0 0 0 0;
            0 (-K23-K20-K24) K32 K42 0 K62 ;
            0 K23 -K32 0 0 0;
            0 K24 0 -K42 0 0;
            K15 0 0 0 -K56 0;
            0 0 0 0 K56 -K62;
            ];
        
        
        % calc new V
        ej=expm(J(:,:)*dt);
        v(:,m)=ej*v(:,m-1);
        
        % assign different dosing absorption rates per dose.
        
        if (m==mda1)
            v(1,m)=dosing*F1*0.577+v(1,m);
        end
        
        if (m==mda2)
            
            v(1,m)=dosing*1.24*F1*0.577+v(1,m);
        end
        
        if (m==mda3)
            
            v(1,m)=dosing*F1*(1.24^2)*0.577+v(1,m);
        end
    end
    
    % keep solution for compartment 2 for each individual
    vv(jj,:)=v(2,:);
    
    %xticktimes=(24:24*5:maxa);
    xticktimes=(10:10:maxa);
    xticklabels=(1:100);
    v(2,end)=NaN;
    patch((0:maxa/asteps:maxa),v(2,:),'black','EdgeColor','black','FaceAlpha',0.2);
    %set(gca,'FontSize',14,'XTick',xticktimes,'XTickLabel',1:2:100)
    set(gca,'FontSize',14,'XTick',xticktimes,'XTickLabel',10:10:1200)
    hold on
end

toc

mvv=vv(1,:);

figure(2)
semilogy((0:maxa/asteps:maxa),mvv,'k','LineWidth',3)
%set(gca,'FontSize',14,'XTick',xticktimes,'XTickLabel',1:5:100)
set(gca,'FontSize',14,'XTick',xticktimes,'XTickLabel',10:10:1200)
ylim([0.1 1000])
xlabel('Time in days')



% dlmwrite('pippk_deterministic.txt',vv)
