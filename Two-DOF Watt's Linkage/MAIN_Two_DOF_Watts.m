close all
clear 
clc

%% Load two-DOF Watt's linkage kinemtics
load('GetGeo_Two_DOF_Watts','Phi1','Phi2','ThetaA','ThetaB','ThetaC','ThetaD','ThetaE','PE_G')

%% Compute rest angles and stiffnesses that optimize for continuous equilibrium
[aA,aB,aC,aD,aE,kA,kB,kC,kD,kE,fval,flag,output,allmins] = FindRest_Two_DOF_Watts(ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,PE_G);

%% Calculate PE from torsional spring j: 1/2*k_j*(theta_j-alpha_j)^2
PE_A = (1/2)*kA*(ThetaA-aA).^2;
PE_B = (1/2)*kB*(ThetaB-aB).^2;
PE_C = (1/2)*kC*(ThetaC-aC).^2;
PE_D = (1/2)*kD*(ThetaD-aD).^2;
PE_E = (1/2)*kE*(ThetaE-aE).^2;
PE_S = PE_A+PE_B+PE_C+PE_D+PE_E;

% Sum all PE contributions
PE_T = PE_S + PE_G;

%% Plot potential energy surface
figure()
surf(Phi1*180/pi,Phi2*180/pi,PE_T)
hold on
contour3(Phi1*180/pi,Phi2*180/pi,PE_T,0:2:300,'k')
caxis([50 250]);
zlim([50 250]);
shading interp
surf(Phi1*180/pi,Phi2*180/pi,PE_G)
contour3(Phi1*180/pi,Phi2*180/pi,PE_G,0:2:300,'k')
caxis([50 250]);
zlim([50 250]);
colormap jet
view([-145 35])
shading interp
colorbar
xlabel('\phi_1')
ylabel('\phi_2')

%% Plot PE breakdown for Phi_1=180
i = 3;
figure()
area(Phi2*180/pi,PE_A(:,i)+PE_B(:,i)+PE_C(:,i)+PE_D(:,i)+PE_E(:,i)+(PE_G(:,i)),'FaceColor',[0.5 0.5 0.5])
hold on
area(Phi2*180/pi,PE_A(:,i)+PE_B(:,i)+PE_C(:,i)+PE_D(:,i)+PE_E(:,i),'FaceColor','m')
area(Phi2*180/pi,PE_A(:,i)+PE_B(:,i)+PE_C(:,i)+PE_D(:,i),'FaceColor','b')
area(Phi2*180/pi,PE_A(:,i)+PE_B(:,i)+PE_C(:,i),'FaceColor','g')
area(Phi2*180/pi,PE_A(:,i)+PE_B(:,i),'FaceColor','y')
area(Phi2*180/pi,PE_A(:,i),'FaceColor','r')

cross_section_PE_T = PE_A(:,i)+PE_B(:,i)+PE_C(:,i)+PE_D(:,i)+PE_E(:,i)+(PE_G(:,i));
cross_section_PE_G = (PE_G(:,i));
cross_section_PEavg = mean(cross_section_PE_T);
cross_section_PE_Gavg = mean(cross_section_PE_G);

N = numel(cross_section_PE_T);

cross_section_RMSE = sqrt(sum(sum((cross_section_PE_T(i)-cross_section_PEavg).^2))/N);
cross_section_RMSEgrav = sqrt(sum(sum((cross_section_PE_G(i)-cross_section_PE_Gavg).^2))/N);

cs_normRMSE = cross_section_RMSE/cross_section_RMSEgrav;

%% Plot PE breakdown for Phi_2=90
j = 3;
figure()
area(Phi1*180/pi,PE_A(j,:)+PE_B(j,:)+PE_C(j,:)+PE_D(j,:)+PE_E(j,:)+(PE_G(j,:)),'FaceColor',[0.5 0.5 0.5])
hold on
area(Phi1*180/pi,PE_A(j,:)+PE_B(j,:)+PE_C(j,:)+PE_D(j,:)+PE_E(j,:),'FaceColor','m')
area(Phi1*180/pi,PE_A(j,:)+PE_B(j,:)+PE_C(j,:)+PE_D(j,:),'FaceColor','b')
area(Phi1*180/pi,PE_A(j,:)+PE_B(j,:)+PE_C(j,:),'FaceColor','g')
area(Phi1*180/pi,PE_A(j,:)+PE_B(j,:),'FaceColor','y')
area(Phi1*180/pi,PE_A(j,:),'FaceColor','r')

cross_section_PE_T = PE_A(j,:)+PE_B(j,:)+PE_C(j,:)+PE_D(j,:)+PE_E(j,:)+(PE_G(j,:));
cross_section_PE_G = (PE_G(j,:));
cross_section_PEavg = mean(cross_section_PE_T);
cross_section_PE_Gavg = mean(cross_section_PE_G);

N = numel(cross_section_PE_T);

cross_section_RMSE = sqrt(sum(sum((cross_section_PE_T(i)-cross_section_PEavg).^2))/N);
cross_section_RMSEgrav = sqrt(sum(sum((cross_section_PE_G(i)-cross_section_PE_Gavg).^2))/N);

cs_normRMSE = cross_section_RMSE/cross_section_RMSEgrav;

%% Calculate RMSE
PEavg = mean(mean(PE_T));
PE_Gavg = mean(mean(PE_G));

N = numel(PE_T);

RMSE = sqrt(sum(sum((PE_T-PEavg).^2))./N);
RMSEgrav = sqrt(sum(sum((PE_G-PE_Gavg).^2))./N);

normRMSE = RMSE/RMSEgrav;
    
    