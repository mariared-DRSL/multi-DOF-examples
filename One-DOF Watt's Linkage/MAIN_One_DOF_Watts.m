%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%------------ Watt's Linkage: Continuous Equilibrium ------------%%%
%%%------------------ Internal Torsional Springs ------------------%%%
%  "Programming Continuous Equilibrium Motions in Multi-DOF Systems" %
%%% ------------------ Redoutey, M., Filipov, E.T. ----------------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; 
clear; clc 

%% Load one-DOF Watt's linkage kinemtics
load('GetGeo_One_DOF_Watts','Phi1','ThetaA','ThetaB','ThetaC','ThetaD','PE_G');

%% Compute rest angles and stiffnesses that optimize for continuous equilibrium
[as,ks,history,exitflag,lambda] = FindRest_One_DOF_Watts(ThetaA,ThetaB,ThetaC,ThetaD,PE_G); 

%% Calculate PE from torsional spring j: 1/2*k_j*(theta_j-alpha_j)^2
PE_A = (1/2)*ks(1)*(ThetaA-(as(1))).^2;
PE_B = (1/2)*ks(2)*(ThetaB-(as(2))).^2;
PE_C = (1/2)*ks(3)*(ThetaC-(as(3))).^2;
PE_D = (1/2)*ks(4)*(ThetaD-(as(4))).^2;
    
% Sum all PE contributions
PE_T = PE_A + PE_B + PE_C + PE_D + PE_G;
    
%% Plot potential energy breakdown 
figure()
hold on
area(Phi1*180/pi,PE_A+PE_B+PE_C+PE_D+(PE_G-min(PE_G)*0),'FaceColor',[0.5 0.5 0.5])
area(Phi1*180/pi,PE_A+PE_B+PE_C+PE_D,'FaceColor',[0 0.4470 0.7410])
area(Phi1*180/pi,PE_A+PE_B+PE_C,'FaceColor',[0.4660 0.6740 0.1880])
area(Phi1*180/pi,PE_A+PE_B,'FaceColor',[0.9290 0.6940 0.1250])
area(Phi1*180/pi,PE_A,'FaceColor','r')
xlim([min(Phi1)*180/pi max(Phi1)*180/pi])
xlabel('\phi')
ylabel('PE')

%% Calculate RMSE
PEavg = mean(PE_T);  
N = length(PE_T);
RMSE = sqrt(sum((PE_T-PEavg).^2)/N);

PE_Gavg = mean(PE_G);
RMSEgrav = sqrt(sum((PE_G-PE_Gavg).^2)/N);

normRMSE = RMSE/RMSEgrav;
