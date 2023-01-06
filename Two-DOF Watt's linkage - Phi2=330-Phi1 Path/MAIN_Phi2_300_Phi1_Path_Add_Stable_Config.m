%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%--------- Two-DOF Watt's Linkage: Continuous Equilibrium -------%%%
%  "Programming Continuous Equilibrium Motions in Multi-DOF Systems" %
%%% ------------------ Redoutey, M., Filipov, E.T. ----------------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear 
clc

%% Load two-DOF Watt's linkage kinemtics
load('GetGeo_Linear_Phi1_Phi2_Path','Phi1','Phi2','PHI1','PHI2','ThetaA','ThetaB','ThetaC','ThetaD','ThetaE','PE_G')

%% Define path
pathPhi1 = Phi1;
pathPhi2 = fliplr(Phi2);

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

%% Load continuous equilibrium spring properties (from MAIN_Phi2_300-Phi1_Path)
aA = 0.1303; aE = 2.1695; aB = 2.6842; aC = 3.1058; aD = 1.8887;
kA = 8.1398; kE = 17.0791; kB = 9.0621; kC = 7.7594; kD = 14.3379;

% Calculate potential energy in continuous equilibrium spring set
PE_A = (1/2)*kA*(ThetaA-aA).^2;
PE_E = (1/2)*kE*(ThetaE-aE).^2;
PE_B = (1/2)*kB*(ThetaB-aB).^2;
PE_C = (1/2)*kC*(ThetaC-aC).^2;
PE_D = (1/2)*kD*(ThetaD-aD).^2;
PE_CE = PE_A+PE_E+PE_B+PE_C+PE_D;

%% Define stable configuration
SS_Phi1 = pathPhi1(21);
SS_Phi2 = pathPhi2(21);

SS_ThetaA = ThetaA(find(Phi2==SS_Phi2),find(Phi1==SS_Phi1));
SS_ThetaE = ThetaE(find(Phi2==SS_Phi2),find(Phi1==SS_Phi1));
SS_ThetaB = ThetaB(find(Phi2==SS_Phi2),find(Phi1==SS_Phi1));
SS_ThetaC = ThetaC(find(Phi2==SS_Phi2),find(Phi1==SS_Phi1));
SS_ThetaD = ThetaD(find(Phi2==SS_Phi2),find(Phi1==SS_Phi1));

SS_aA = SS_ThetaA; SS_aE = SS_ThetaE; SS_aB = SS_ThetaB; SS_aC = SS_ThetaC; SS_aD = SS_ThetaD;
SS_kA = 20; SS_kE = 20; SS_kB = 20; SS_kC = 20; SS_kD = 20;

% Compute potential energy in stable configuration springs
SS_PE_A = (1/2)*SS_kA*(ThetaA-SS_aA).^2;
SS_PE_E = (1/2)*SS_kE*(ThetaE-SS_aE).^2;
SS_PE_B = (1/2)*SS_kB*(ThetaB-SS_aB).^2;
SS_PE_C = (1/2)*SS_kC*(ThetaC-SS_aC).^2;
SS_PE_D = (1/2)*SS_kD*(ThetaD-SS_aD).^2;
PE_SS = SS_PE_A+SS_PE_E+SS_PE_B+SS_PE_C+SS_PE_D;

%% Total potential energy 
PE_T = PE_CE + PE_SS + PE_G;

%% Potential energy along path 
pathPE_A = (1/2)*kA*(pathThetaA-aA).^2;
pathPE_E = (1/2)*kE*(pathThetaE-aE).^2;
pathPE_B = (1/2)*kB*(pathThetaB-aB).^2;
pathPE_C = (1/2)*kC*(pathThetaC-aC).^2;
pathPE_D = (1/2)*kD*(pathThetaD-aD).^2;

SS_pathPE_A = (1/2)*SS_kA*(pathThetaA-SS_aA).^2;
SS_pathPE_E = (1/2)*SS_kE*(pathThetaE-SS_aE).^2;
SS_pathPE_B = (1/2)*SS_kB*(pathThetaB-SS_aB).^2;
SS_pathPE_C = (1/2)*SS_kC*(pathThetaC-SS_aC).^2;
SS_pathPE_D = (1/2)*SS_kD*(pathThetaD-SS_aD).^2;
SS_pathPE = SS_pathPE_A+SS_pathPE_E+SS_pathPE_B+SS_pathPE_C+SS_pathPE_D;

pathPE_T = pathPE_A + pathPE_E + pathPE_B + pathPE_C + pathPE_D + SS_pathPE + pathPE_G;

%% Plot spring kinematics
figure('Position',[-836,380,560,420])
plot(Phi1*180/pi,pathThetaA*180/pi,'r','LineWidth',2)
hold on
plot(Phi1*180/pi,pathThetaE*180/pi,'m','LineWidth',2)
plot(Phi1*180/pi,pathThetaB*180/pi,'y','LineWidth',2)
plot(Phi1*180/pi,pathThetaC*180/pi,'g','LineWidth',2)
plot(Phi1*180/pi,pathThetaD*180/pi,'b','LineWidth',2)
xlabel('\phi_1')
xlim([min(Phi1)*180/pi max(Phi1)*180/pi])
legend('\phi_1','\phi_2','\phi_3','\phi_4','\phi_5')

%% Plot potential energy along the path
figure('Position',[-836,380,560,420])
area(pathPE_A + pathPE_E + pathPE_B + pathPE_C + pathPE_D + SS_pathPE + pathPE_G,'FaceColor',[0.5 0.5 0.5])
hold on
area(pathPE_A + pathPE_E + pathPE_B + pathPE_C + pathPE_D + SS_pathPE,'FaceColor',[0.5 0.5 1])
area(pathPE_A + pathPE_E + pathPE_B + pathPE_C + pathPE_D,'FaceColor','b')
area(pathPE_A + pathPE_E + pathPE_B + pathPE_C,'FaceColor','g')
area(pathPE_A + pathPE_E + pathPE_B,'FaceColor','y')
area(pathPE_A + pathPE_E,'FaceColor','m')
area(pathPE_A,'FaceColor','r')
xlabel('\phi_1')
ylabel('PE_T')
xlim([1 51])
ylim([0 300])

%% Plot potential energy perpendicular to the path 
p1Phi1 = [170:220]*pi/180; p1Phi2 = [80:130]*pi/180;
p2Phi1 = [196:220]*pi/180; p2Phi2 = [106:130]*pi/180;

jj = 11:51; ii = 1:41;

for m = 1:length(jj)
    PE_perp(m) = PE_T(jj(m),ii(m));
    perpPE_G(m) = PE_G(jj(m),ii(m));
    perpPE_A(m) = PE_A(jj(m),ii(m));
    perpPE_E(m) = PE_E(jj(m),ii(m));
    perpPE_B(m) = PE_B(jj(m),ii(m));
    perpPE_C(m) = PE_C(jj(m),ii(m));
    perpPE_D(m) = PE_D(jj(m),ii(m));
    perpPE_SS(m) = PE_SS(jj(m),ii(m));
end

figure('Position',[-836,380,560,420])
area(1:41,[(perpPE_A + perpPE_E + perpPE_B + perpPE_C + perpPE_D + perpPE_SS + perpPE_G)],'FaceColor',[0.5 0.5 0.5])
hold on
area(1:41,[(perpPE_A + perpPE_E + perpPE_B + perpPE_C + perpPE_D + perpPE_SS)],'FaceColor',[0.5 0.5 1])
area(1:41,[(perpPE_A + perpPE_E + perpPE_B + perpPE_C + perpPE_D)],'FaceColor','b')
area(1:41,[(perpPE_A + perpPE_E + perpPE_B + perpPE_C) ],'FaceColor','g')
area(1:41,[(perpPE_A + perpPE_E + perpPE_B)],'FaceColor','y')
area(1:41,[(perpPE_A + perpPE_E)],'FaceColor','m')
area(1:41,[(perpPE_A) ],'FaceColor','r')
xline(21)
xlim([1 51])
ylim([0 300])

%% Plot potential energy surface 
figure('Position',[-836,380,560,420])
surf(Phi1*180/pi,Phi2*180/pi,PE_T,'FaceAlpha',0.75)
hold on
contour3(Phi1*180/pi,Phi2*180/pi,PE_T,100:1:1000,'k')
colormap('jet')
colorbar;
surf(Phi1*180/pi,Phi2*180/pi,PE_G,'FaceAlpha',0.75)
contour3(Phi1*180/pi,Phi2*180/pi,PE_G,100:1:1000,'k')
plot3(pathPhi1*180/pi,pathPhi2*180/pi,pathPE_T,'k','LineWidth',2)
scatter3(SS_Phi1*180/pi,SS_Phi2*180/pi,PE_T(21,21),'kx')
shading interp
xlabel('\phi_1')
ylabel('\phi_2')
zlim([100 300])
caxis([100 300])
view([-125 41])

%% Plot gradient of potential energy
figure('Position',[-836,380,560,420])
[Fx,Fy] = gradient(PE_G);
quiver(Phi1(1:5:end,1:5:end)*180/pi,Phi2(1:5:end,1:5:end)*180/pi,Fx(1:5:end,1:5:end),Fy(1:5:end,1:5:end))

figure('Position',[-836,380,560,420])
[Fx,Fy] = gradient(PE_T);
quiver(Phi1(1:5:end,1:5:end)*180/pi,Phi2(1:5:end,1:5:end)*180/pi,Fx(1:5:end,1:5:end),Fy(1:5:end,1:5:end))

%% RMSE
PEavg_path = mean(pathPE_T);
PE_G_avg_path = mean(pathPE_G);

N=0;
for i = 1:length(pathPE_T)
    N=N+1;
    dev = pathPE_T(i)-PEavg_path;
    dev2(N) = dev^2;
    devG = pathPE_G(i)-PE_G_avg_path;
    devG2(N) = devG^2;
end

sumup = sum(dev2);
RMSE = sqrt(sumup/N);
RMSEgrav = sqrt(sum(devG2)/N);
normRMSE = RMSE/RMSEgrav;