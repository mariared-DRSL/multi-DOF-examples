%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%--------- Two-DOF Watt's Linkage: Continuous Equilibrium -------%%%
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

%% Load continuous equilibrium spring properties (from MAIN_Phi2_90_Path)
aA = 4.3567; aE = 1.6266; aB = 0.2901; aC = 3.2876; aD = 4.1467;
kA = 0.0331; kE = 199.6282; kB = 11.7438; kC = 19.9531; kD = 0.0257;

% Calculate potential energy in continuous equilibrium spring set
PE_A = (1/2)*kA*(ThetaA-aA).^2;
PE_E = (1/2)*kE*(ThetaE-aE).^2;
PE_B = (1/2)*kB*(ThetaB-aB).^2;
PE_C = (1/2)*kC*(ThetaC-aC).^2;
PE_D = (1/2)*kD*(ThetaD-aD).^2;
PE_CE = PE_A+PE_E+PE_B+PE_C+PE_D;

%% Define stable configuration
SS_Phi1 = 210*pi/180;
SS_Phi2 = 90*pi/180;

SS_ThetaA = ThetaA(find(Phi2==SS_Phi2),find(Phi1==SS_Phi1));
SS_ThetaE = ThetaE(find(Phi2==SS_Phi2),find(Phi1==SS_Phi1));
SS_ThetaB = ThetaB(find(Phi2==SS_Phi2),find(Phi1==SS_Phi1));
SS_ThetaC = ThetaC(find(Phi2==SS_Phi2),find(Phi1==SS_Phi1));
SS_ThetaD = ThetaD(find(Phi2==SS_Phi2),find(Phi1==SS_Phi1));

SS_aA = SS_ThetaA; SS_aE = SS_ThetaE; SS_aB = SS_ThetaB; SS_aC = SS_ThetaC; SS_aD = SS_ThetaD;
SS_kA = 20; SS_kE = 20; SS_kB = 20; SS_kC = 20; SS_kD = 20;

% Compute potential energy in stable configuration springs
PE_SS = (1/2)*SS_kA*(ThetaA-SS_aA).^2 + (1/2)*SS_kE*(ThetaE-SS_aE).^2 + (1/2)*SS_kB*(ThetaB-SS_aB).^2 + (1/2)*SS_kC*(ThetaC-SS_aC).^2 + (1/2)*SS_kD*(ThetaD-SS_aD).^2;

%% Total potential energy 
PE_T = PE_CE + PE_SS + PE_G;

%% Potential energy along path 
pathPE_SS = (1/2)*SS_kA*(pathThetaA-SS_aA).^2 + (1/2)*SS_kE*(pathThetaE-SS_aE).^2 + (1/2)*SS_kB*(pathThetaB-SS_aB).^2 + (1/2)*SS_kC*(pathThetaC-SS_aC).^2 + (1/2)*SS_kD*(pathThetaD-SS_aD).^2;

pathPE_A = (1/2)*kA*(pathThetaA-aA).^2;
pathPE_E = (1/2)*kE*(pathThetaE-aE).^2;
pathPE_B = (1/2)*kB*(pathThetaB-aB).^2;
pathPE_C = (1/2)*kC*(pathThetaC-aC).^2;
pathPE_D = (1/2)*kD*(pathThetaD-aD).^2;

pathPE_T = pathPE_A + pathPE_E + pathPE_B + pathPE_C + pathPE_D + pathPE_SS + pathPE_G;

SS_PE_T = PE_T(find(Phi2==SS_Phi2),find(Phi1==SS_Phi1));

%% Plot potential energy surface 
figure()
surf(Phi1*180/pi,Phi2*180/pi,PE_T,'FaceAlpha',0.75)
hold on
contour3(Phi1*180/pi,Phi2*180/pi,PE_T,100:2:1000,'k')
colormap('jet')
colorbar;
surf(Phi1*180/pi,Phi2*180/pi,PE_G,'FaceAlpha',0.75)
contour3(Phi1*180/pi,Phi2*180/pi,PE_G,100:2:1000,'k')
plot3(pathPhi1*180/pi,pathPhi2*180/pi,pathPE_T,'k','LineWidth',2)
scatter3(SS_Phi1*180/pi,SS_Phi2*180/pi,SS_PE_T,'kx')
shading interp
xlabel('\phi_1')
ylabel('\phi_2')
view([-125 41])
zlim([100 300])
caxis([100 300])

%%
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
    p1PE_SS(k) = PE_SS(j,i);
end 

%% Plot potential energy perpendicular to the path 
figure()
area(1:51,[(p1PE_A + p1PE_E + p1PE_B + p1PE_C + p1PE_D + p1PE_SS + p1PE_G)],'FaceColor',[0.5 0.5 0.5])
hold on
area(1:51,[(p1PE_A + p1PE_E + p1PE_B + p1PE_C + p1PE_D + p1PE_SS)],'FaceColor',[0.8 0.8 0.8])
area(1:51,[(p1PE_A + p1PE_E + p1PE_B + p1PE_C + p1PE_D)],'FaceColor','b')
area(1:51,[(p1PE_A + p1PE_E + p1PE_B + p1PE_C)],'FaceColor','g')
area(1:51,[(p1PE_A + p1PE_E + p1PE_B)],'FaceColor','y')
area(1:51,[(p1PE_A + p1PE_E)],'FaceColor','m')
area(1:51,[(p1PE_A)],'FaceColor','r')
xlim([1 51])
ylim([0 300])
xline(find(SS_Phi2==p1Phi2))
xlabel('\phi_2')

%% Plot potential energy along the path
figure()
area(1:51,[(pathPE_A + pathPE_E + pathPE_B + pathPE_C + pathPE_D + pathPE_SS + pathPE_G)],'FaceColor',[0.5 0.5 0.5])
hold on
area(1:51,[(pathPE_A + pathPE_E + pathPE_B + pathPE_C + pathPE_D + pathPE_SS)],'FaceColor',[0.8 0.8 0.8])
area(1:51,[(pathPE_A + pathPE_E + pathPE_B + pathPE_C + pathPE_D)],'FaceColor','b')
area(1:51,[(pathPE_A + pathPE_E + pathPE_B + pathPE_C)],'FaceColor','g')
area(1:51,[(pathPE_A + pathPE_E + pathPE_B)],'FaceColor','y')
area(1:51,[(pathPE_A + pathPE_E)],'FaceColor','m')
area(1:51,[(pathPE_A)],'FaceColor','r')
xlim([1 51])
ylim([0 300])
xlabel('\phi_1')
xline(find(SS_Phi1==Phi1))

