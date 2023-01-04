%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%-------- Five-Fold Origami Vertex: Continuous Equilibrium ------%%%
%  "Programming Continuous Equilibrium Motions in Multi-DOF Systems" %
%%% ------------------ Redoutey, M., Filipov, E.T. ----------------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
clc

%% Load five-fold origami vertex kinemtics
load('FiveFoldVertex_Kinematics.mat','PHI1','PHI2','ThetaA','ThetaB','ThetaC','ThetaD','Theta9','PE_G')

%% Compute rest angles and stiffnesses that optimize for continuous equilibrium
[aA,aB,aC,aD,kA,kB,kC,kD,history1,exitflag1,output1,lambda1] = FindRest_FiveFoldVertex(ThetaA,ThetaB,ThetaC,ThetaD,PE_G);

%% Calculate PE from torsional spring j: 1/2*k_j*(theta_j-alpha_j)^2
PE_A = (1/2)*kA*(ThetaA-aA).^2;
PE_B = (1/2)*kB*(ThetaB-aB).^2;
PE_C = (1/2)*kC*(ThetaC-aC).^2;
PE_D = (1/2)*kD*(ThetaD-aD).^2;

% Sum all PE contributions
PE_T = PE_A + PE_B + PE_C + PE_D + PE_G;

%% Plot potential energy surface
inc = 1; inc2=1;
figure()
surf(PHI1([1:inc:end])*180/pi,PHI2([1:inc2:end])*180/pi,PE_G([1:inc:end],[1:inc2:end])','FaceAlpha',0.75)
shading interp
hold on
contour3(PHI1(1:inc:end)*180/pi,PHI2(1:inc2:end)*180/pi,PE_G(1:inc:end,1:inc2:end)',[0:1:16],'k')
xlabel('\phi_1')
ylabel('\phi_2')
zlim([0 16])
caxis([0 16])
surf(PHI1(1:inc:end)*180/pi,PHI2(1:inc2:end)*180/pi,PE_T(1:inc:end,1:inc2:end)','FaceAlpha',0.75)
shading interp
contour3(PHI1(1:inc:end)*180/pi,PHI2(1:inc2:end)*180/pi,PE_T(1:inc:end,1:inc2:end)',[0:1:16],'k')
colorbar
colormap jet
view([-130 30])

%% Plot potential energy breakdown for Phi1 = 150
p2 = 88;

figure()
area(PHI1*180/pi,PE_A(:,p2)+PE_B(:,p2)+PE_C(:,p2)+PE_D(:,p2)+PE_G(:,p2),'FaceColor',[0.5 0.5 0.5])
hold on
area(PHI1*180/pi,PE_A(:,p2)+PE_B(:,p2)+PE_C(:,p2)+PE_D(:,p2),'FaceColor','b')
area(PHI1*180/pi,PE_A(:,p2)+PE_B(:,p2)+PE_C(:,p2),'FaceColor','g')
area(PHI1*180/pi,PE_A(:,p2)+PE_B(:,p2),'FaceColor','y')
area(PHI1*180/pi,PE_A(:,p2),'FaceColor','r')
xlim([min(PHI1*180/pi) max(PHI1*180/pi)])
xlabel('\phi_1')
ylabel('PE')

%% Plot potential energy breakdown for Phi2 = 50
p1 = 141; 

figure()
area(PHI2*180/pi,PE_A(p1,:)+PE_B(p1,:)+PE_C(p1,:)+PE_D(p1,:)+PE_G(p1,:),'FaceColor',[0.5 0.5 0.5])
hold on
area(PHI2*180/pi,PE_A(p1,:)+PE_B(p1,:)+PE_C(p1,:)+PE_D(p1,:),'FaceColor','b')
area(PHI2*180/pi,PE_A(p1,:)+PE_B(p1,:)+PE_C(p1,:),'FaceColor','g')
area(PHI2*180/pi,PE_A(p1,:)+PE_B(p1,:),'FaceColor','y')
area(PHI2*180/pi,PE_A(p1,:),'FaceColor','r')
xlim([min(PHI2*180/pi) max(PHI2*180/pi)])
xlabel('\phi_2')
ylabel('PE')

%% Compute RMSE 
PEavg = mean(mean(PE_T(~isnan(PE_T))));
PE_Gavg = mean(mean(PE_G(~isnan(PE_G))));
    
N=0;
for i = 1:size(PE_T,1)
    for j = 1:size(PE_T,2)
        if ~isnan(PE_T(i,j))
            N=N+1;
            dev = PE_T(i,j)-PEavg;
            dev2(N) = dev^2;
            devG = PE_G(i,j)-PE_Gavg;
            dev2G(N) = devG^2;
        end
    end
end
    
sumup = sum(dev2);
RMSE = sqrt(sumup/N);
RMSEgrav = sqrt(sum(dev2G)/N);
normRMSE = RMSE/RMSEgrav;