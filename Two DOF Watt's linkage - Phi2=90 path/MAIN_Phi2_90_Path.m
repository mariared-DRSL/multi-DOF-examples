%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%------------ Watt's Linkage: Continuous Equilibrium ------------%%%
%%%------------------ Internal Torsional Springs ------------------%%%
%  "Programming Continuous Equilibrium Motions in Multi-DOF Systems" %
%%% ------------------ Redoutey, M., Filipov, E.T. ----------------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear 
clc

%% Load two-DOF Watt's linkage kinemtics
load('GetGeo_Phi2_90_Path','Phi1','Phi2','PHI1','PHI2','ThetaA','ThetaB','ThetaC','ThetaD','ThetaE','PE_G');

%% Define path
pathPhi1 = Phi1;
pathPhi2 = repmat(1*pi/2,size(pathPhi1));

% Get PE_G along path
for i = 1:length(pathPhi1)
    j = find(round(Phi2,4)==round(pathPhi2(i),4));
    pathPE_G(i) = PE_G(j,i);
    pathThetaA(i) = ThetaA(j,i);
    pathThetaB(i) = ThetaB(j,i);
    pathThetaC(i) = ThetaC(j,i);
    pathThetaD(i) = ThetaD(j,i);
    pathThetaE(i) = ThetaE(j,i);
end 

%% Optimization
[rest,history,exitflag,output,lambda] = FindRest_Phi2_90_Path(PHI1,PHI2,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,PE_G);

aA = rest(1); aE = rest(2); aB = rest(3); aC = rest(4); aD = rest(5);
kA = rest(6); kE = rest(7); kB = rest(8); kC = rest(9); kD = rest(10);

%% Calculate PE from torsional spring j: 1/2*k_j*(theta_j-alpha_j)^2
PE_A = (1/2)*kA*(ThetaA-aA).^2;
PE_E = (1/2)*kE*(ThetaE-aE).^2;
PE_B = (1/2)*kB*(ThetaB-aB).^2;
PE_C = (1/2)*kC*(ThetaC-aC).^2;
PE_D = (1/2)*kD*(ThetaD-aD).^2;
PE_S = PE_A+PE_E+PE_B+PE_C+PE_D;

PE_T = PE_S + PE_G;

%% Calculate PE along path
pathPE_A = (1/2)*kA*(pathThetaA-aA).^2;
pathPE_E = (1/2)*kE*(pathThetaE-aE).^2;
pathPE_B = (1/2)*kB*(pathThetaB-aB).^2;
pathPE_C = (1/2)*kC*(pathThetaC-aC).^2;
pathPE_D = (1/2)*kD*(pathThetaD-aD).^2;

pathPE_T = pathPE_A + pathPE_E + pathPE_B + pathPE_C + pathPE_D + pathPE_G;

%% Plot spring kinematics
figure()
plot(Phi1*180/pi,pathThetaA*180/pi,'r','LineWidth',2)
hold on
plot(Phi1*180/pi,pathThetaE*180/pi,'m','LineWidth',2)
plot(Phi1*180/pi,pathThetaB*180/pi,'y','LineWidth',2)
plot(Phi1*180/pi,pathThetaC*180/pi,'g','LineWidth',2)
plot(Phi1*180/pi,pathThetaD*180/pi,'b','LineWidth',2)
xlabel('\phi_1')
xlim([min(Phi1)*180/pi max(Phi1)*180/pi])
legend('\phi_1','\phi_2','\phi_3','\phi_4','\phi_5')

%% Plot potential energy surface 
figure()
surf(Phi1*180/pi,Phi2*180/pi,PE_T,'FaceAlpha',0.75)
hold on
contour3(Phi1*180/pi,Phi2*180/pi,PE_T,100:5:300,'k')
colormap('jet')
colorbar;
surf(Phi1*180/pi,Phi2*180/pi,PE_G,'FaceAlpha',0.75)
contour3(Phi1*180/pi,Phi2*180/pi,PE_G,100:5:300,'k')
plot3(pathPhi1*180/pi,pathPhi2*180/pi,pathPE_T,'k','LineWidth',2)
shading interp
xlabel('\phi_1')
ylabel('\phi_2')
zlim([100 300])
caxis([100 300])
view([-125 41])

%% Plot gradient of PE_T
GetQuiver(Phi1,pathPhi1,Phi2,pathPhi2,PE_T)

%% Plot PE breakdown along path
figure()
area(1:51,[(pathPE_A + pathPE_E + pathPE_B + pathPE_C + pathPE_D + pathPE_G)],'FaceColor',[0.5 0.5 0.5])
hold on
area(1:51,[(pathPE_A + pathPE_E + pathPE_B + pathPE_C + pathPE_D)],'FaceColor','b')
area(1:51,[(pathPE_A + pathPE_E + pathPE_B + pathPE_C)],'FaceColor','g')
area(1:51,[(pathPE_A + pathPE_E + pathPE_B)],'FaceColor','y')
area(1:51,[(pathPE_A + pathPE_E)],'FaceColor','m')
area(1:51,[(pathPE_A)],'FaceColor','r')
xlim([1 51])

%% Plot PE breakdown perpendicular to the path
p1Phi1 = repmat(pi,1,51); p1Phi2 = [130:-1:80]*pi/180;
p2Phi1 = [180:220]*pi/180; p2Phi2 = repmat(pi/2,1,41);

for k = 1:length(p1Phi1)
    i = find(round(Phi1,4)==round(p1Phi1(k),4));
    j = find(round(Phi2,4)==round(p1Phi2(k),4));
    p1PE_G(k) = PE_G(j,i);
    p1PE_A(k) = PE_A(j,i);
    p1PE_E(k) = PE_E(j,i);
    p1PE_B(k) = PE_B(j,i);
    p1PE_C(k) = PE_C(j,i);
    p1PE_D(k) = PE_D(j,i);
end 

figure()
area(1:51,[(p1PE_A + p1PE_E + p1PE_B + p1PE_C + p1PE_D + p1PE_G)],'FaceColor',[0.5 0.5 0.5])
hold on
area(1:51,[(p1PE_A + p1PE_E + p1PE_B + p1PE_C + p1PE_D)],'FaceColor','b')
area(1:51,[(p1PE_A + p1PE_E + p1PE_B + p1PE_C)],'FaceColor','g')
area(1:51,[(p1PE_A + p1PE_E + p1PE_B)],'FaceColor','y')
area(1:51,[(p1PE_A + p1PE_E)],'FaceColor','m')
area(1:51,[(p1PE_A)],'FaceColor','r')
xlim([1 51])

%% Calculate RMSE
PEavg_path = mean(pathPE_T);
PE_G_avg_path = mean(pathPE_G);

N = numel(pathPE_T);

RMSE = sqrt(sum((pathPE_T-PEavg_path).^2)/N);
RMSEgrav = sqrt(sum((pathPE_G-PE_G_avg_path).^2)/N);

normRMSE = RMSE/RMSEgrav;
