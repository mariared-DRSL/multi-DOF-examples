%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%-------- Three-DOF Watt's Linkage: Continuous Equilibrium ------%%%
%%%------------------ Internal Torsional Springs ------------------%%%
%  "Programming Continuous Equilibrium Motions in Multi-DOF Systems" %
%%% ------------------ Redoutey, M., Filipov, E.T. ----------------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear 
clc

%% Load one-DOF Watt's linkage kinemtics
load('GetGeo_Three_DOF_Watts','Phi1','Phi2','Phi3','ThetaA','ThetaB','ThetaC','ThetaD','ThetaE','ThetaF','PE_G')

%% Compute rest angles and stiffnesses that optimize for continuous equilibrium
[aA,aB,aC,aD,aE,aF,kA,kB,kC,kD,kE,kF,history,exitflag,output,lambda] = FindRest_Three_DOF_Watts(ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,ThetaF,PE_G);

%% Calculate PE from torsional spring j: 1/2*k_j*(theta_j-alpha_j)^2
PE_A = (1/2)*kA*(ThetaA-aA).^2;
PE_B = (1/2)*kB*(ThetaB-aB).^2;
PE_C = (1/2)*kC*(ThetaC-aC).^2;
PE_D = (1/2)*kD*(ThetaD-aD).^2;
PE_E = (1/2)*kE*(ThetaE-aE).^2;
PE_F = (1/2)*kF*(ThetaF-aF).^2;
PE_S = PE_A+PE_B+PE_C+PE_D+PE_E+PE_F;

PE_T = PE_S + PE_G;

%% Plot PE_G volume
figure()
phi3slice = 100; phi2slice = 130; phi1slice = 210;
slice(Phi3*180/pi,Phi2*180/pi,Phi1*180/pi,PE_G,phi3slice,phi2slice,phi1slice)
shading interp
s = contourslice(Phi3*180/pi,Phi2*180/pi,Phi1*180/pi,PE_G,phi3slice,phi2slice,phi1slice,0:2:500);
set(s,'EdgeColor','k')
colorbar
colormap jet
lighting none
ylim([80 130])
xlim([100 150])
zlim([170 210])
view([-120 30])
caxis([200 450])
xlabel('\phi_3')
ylabel('\phi_2')
zlabel('\phi_1')

%% Plot PE_T volume
figure()
phi3slice = 100; phi2slice = 130; phi1slice = 210;
slice(Phi3*180/pi,Phi2*180/pi,Phi1*180/pi,PE_T,phi3slice,phi2slice,phi1slice)
shading interp
s = contourslice(Phi3*180/pi,Phi2*180/pi,Phi1*180/pi,PE_T,phi3slice,phi2slice,phi1slice,0:2:500);
set(s,'EdgeColor','k')
colorbar
colormap jet
lighting none
ylim([80 130])
xlim([100 150])
zlim([170 210])
view([-120 30])
caxis([200 450])
xlabel('\phi_3')
ylabel('\phi_2')
zlabel('\phi_1')

%% Plot potential energy surface for Phi1 = 180
figure()
idx1 = 1:6; idx2 = 1:6; idx3 = 3;
surf(Phi3*180/pi,Phi2*180/pi,PE_T(idx1,idx2,idx3))
hold on
contour3(Phi3*180/pi,Phi2*180/pi,PE_T(idx1,idx2,idx3),0:2:500,'k')
surf(Phi3*180/pi,Phi2*180/pi,PE_G(idx1,idx2,idx3))
contour3(Phi3*180/pi,Phi2*180/pi,PE_G(idx1,idx2,idx3),0:2:500,'k')
xlabel('\phi_3')
ylabel('\phi_2')
zlabel('PE')
xlim([100 150])
ylim([80 130])
zlim([200 450])
caxis([200 450])
shading interp
colormap jet
colorbar
view([-120 30])

%% Plot potential energy breakdown for Phi1 = 180, Phi3 = 110
figure()
idx1 = 1:6; idx2 = 2; idx3 = 3;
area(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)+PE_C(idx1,idx2,idx3)+PE_D(idx1,idx2,idx3)+PE_E(idx1,idx2,idx3)+PE_F(idx1,idx2,idx3)+PE_G(idx1,idx2,idx3),'FaceColor',[0.5 0.5 0.5])
hold on
area(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)+PE_C(idx1,idx2,idx3)+PE_D(idx1,idx2,idx3)+PE_E(idx1,idx2,idx3)+PE_F(idx1,idx2,idx3),'FaceColor','m')
area(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)+PE_C(idx1,idx2,idx3)+PE_D(idx1,idx2,idx3)+PE_E(idx1,idx2,idx3),'FaceColor',[102/255 51/255 153/255])
area(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)+PE_C(idx1,idx2,idx3)+PE_D(idx1,idx2,idx3),'FaceColor','b')
area(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)+PE_C(idx1,idx2,idx3),'FaceColor','g')
area(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3),'FaceColor','y')
area(PE_A(idx1,idx2,idx3),'FaceColor','r')
xlabel('\phi_2')
ylabel('PE')
ylim([0 450])

%% Plot potential energy surface for Phi2 = 90
figure()
idx1 = 2; idx2 = 1:6; idx3 = 1:11;
surf(Phi1*180/pi,Phi3*180/pi,squeeze(PE_T(idx1,idx2,idx3)))
hold on
contour3(Phi1*180/pi,Phi3*180/pi,squeeze(PE_T(idx1,idx2,idx3)),0:2:500,'k')
surf(Phi1*180/pi,Phi3*180/pi,squeeze(PE_G(idx1,idx2,idx3)))
contour3(Phi1*180/pi,Phi3*180/pi,squeeze(PE_G(idx1,idx2,idx3)),0:2:500,'k')
xlabel('\phi_1')
ylabel('\phi_3')
zlabel('PE')
xlim([170 210])
ylim([100 150])
zlim([200 450])
caxis([200 450])
shading interp
colormap jet
colorbar
view([-120 30])

%% Plot potential energy breakdown for Phi2 = 90, Phi1 = 200
figure()
idx1 = 2; idx2 = 1:6; idx3 = 7;
area(squeeze(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)+PE_C(idx1,idx2,idx3)+PE_D(idx1,idx2,idx3)+PE_E(idx1,idx2,idx3)+PE_F(idx1,idx2,idx3)+PE_G(idx1,idx2,idx3)),'FaceColor',[0.5 0.5 0.5])
hold on
area(squeeze(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)+PE_C(idx1,idx2,idx3)+PE_D(idx1,idx2,idx3)+PE_E(idx1,idx2,idx3)+PE_F(idx1,idx2,idx3)),'FaceColor','m')
area(squeeze(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)+PE_C(idx1,idx2,idx3)+PE_D(idx1,idx2,idx3)+PE_E(idx1,idx2,idx3)),'FaceColor',[102/255 51/255 153/255])
area(squeeze(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)+PE_C(idx1,idx2,idx3)+PE_D(idx1,idx2,idx3)),'FaceColor','b')
area(squeeze(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)+PE_C(idx1,idx2,idx3)),'FaceColor','g')
area(squeeze(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)),'FaceColor','y')
area(squeeze(PE_A(idx1,idx2,idx3)),'FaceColor','r')
xlabel('\phi_3')
ylabel('PE')
ylim([0 450])

%% Plot potential energy surface for Phi3 = 120
figure()
idx1 = 1:6; idx2 = 3; idx3 = 1:11;
surf(Phi1*180/pi,Phi2*180/pi,squeeze(PE_T(idx1,idx2,idx3)))
hold on
contour3(Phi1*180/pi,Phi2*180/pi,squeeze(PE_T(idx1,idx2,idx3)),0:2:500,'k')
surf(Phi1*180/pi,Phi2*180/pi,squeeze(PE_G(idx1,idx2,idx3)))
contour3(Phi1*180/pi,Phi2*180/pi,squeeze(PE_G(idx1,idx2,idx3)),0:2:500,'k')
xlabel('\phi_1')
ylabel('\phi_2')
zlabel('PE')
xlim([170 210])
ylim([80 130])
zlim([200 450])
caxis([200 450])
shading interp
colormap jet
colorbar
view([-120 30])

%% Plot potential energy breakdown for Phi3 = 120, Phi2 = 100
figure()
idx1 = 3; idx2 = 3; idx3 = 1:11;
area(squeeze(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)+PE_C(idx1,idx2,idx3)+PE_D(idx1,idx2,idx3)+PE_E(idx1,idx2,idx3)+PE_F(idx1,idx2,idx3)+PE_G(idx1,idx2,idx3)),'FaceColor',[0.5 0.5 0.5])
hold on
area(squeeze(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)+PE_C(idx1,idx2,idx3)+PE_D(idx1,idx2,idx3)+PE_E(idx1,idx2,idx3)+PE_F(idx1,idx2,idx3)),'FaceColor','m')
area(squeeze(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)+PE_C(idx1,idx2,idx3)+PE_D(idx1,idx2,idx3)+PE_E(idx1,idx2,idx3)),'FaceColor',[102/255 51/255 153/255])
area(squeeze(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)+PE_C(idx1,idx2,idx3)+PE_D(idx1,idx2,idx3)),'FaceColor','b')
area(squeeze(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)+PE_C(idx1,idx2,idx3)),'FaceColor','g')
area(squeeze(PE_A(idx1,idx2,idx3)+PE_B(idx1,idx2,idx3)),'FaceColor','y')
area(squeeze(PE_A(idx1,idx2,idx3)),'FaceColor','r')
xlabel('\phi_1')
ylabel('PE')
ylim([0 450])

%% Calculate RMSE
PEavg = mean(mean(mean(PE_T)));
RMSE = sqrt(mean(mean(mean((PE_T-PEavg).^2))));

PE_Gavg = mean(mean(mean(PE_G)));
RMSEgrav = sqrt(mean(mean(mean((PE_G-PE_Gavg).^2))));

normRMSE = RMSE/RMSEgrav;
