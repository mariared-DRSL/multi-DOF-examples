%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%--------------- Four Bar Linkages: Zero Stiffness --------------%%%
%%%---------------------- Maria Redoutey --------------------------%%%
%%%---------------------- 12 April 2021 ---------------------------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; 
clear; clc 

set(0, 'defaultFigureRenderer', 'painters')

%% Path 1
PHI = [50:170]*pi/180; 
PSI = linspace(120,170,length(PHI))*pi/180;

pathPHI = PHI;
pathPSI = PSI;

[ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,PE_G,L1,L2] = GetGeo_Excavator(PSI,PHI);

pathThetaA = diag(ThetaA); pathThetaB = diag(ThetaB); pathThetaC = diag(ThetaC); pathThetaD = diag(ThetaD); pathThetaE = diag(ThetaE); pathL1 = diag(L1); pathL2 = diag(L2); pathPE_G = diag(PE_G);

[rest,history,exitflag,output,lambda] = FindRest_Excavator_Path_1(PHI,PSI,pathThetaA,pathThetaB,pathThetaC,pathThetaD,pathThetaE,pathL1,pathL2,pathPE_G,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,L1,L2,PE_G);
% [rest,output,allmins] = FindRest_Excavator_Path_1_GS(PHI,PSI,pathThetaA,pathThetaB,pathThetaC,pathThetaD,pathThetaE,pathL1,pathPE_G,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,L1,PE_G);

%% Path 2
% PHI = [50:170]*pi/180;
% PSI = linspace(120,170,length(PHI))*pi/180;
% 
% [ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,PE_G,L1] = GetGeo_Excavator(PSI,PHI);
% 
% pathPHI = [170:-1:50]*pi/180;
% pathPSI = PSI;
% 
% for m = 1:length(pathPHI)
%     i = find(pathPHI(m) == PHI);
%     j = find(pathPSI(m) == PSI);
%     pathPE_G(m) = PE_G(j,i);
% end
% 
% pathThetaA = diag(flipud(ThetaA)); pathThetaB = diag(flipud(ThetaB)); pathThetaC = diag(flipud(ThetaC)); pathThetaD = diag(flipud(ThetaD)); pathThetaE = diag(flipud(ThetaE)); pathL1 = diag(flipud(L1)); %pathPE_G = diag(flipud(PE_G));
% 
% [rest,history,exitflag,output,lambda] = FindRest_Excavator_Path_2(PHI,PSI,pathThetaA,pathThetaB,pathThetaC,pathThetaD,pathThetaE,pathL1,pathPE_G,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,L1,PE_G);

%%
% Calculate PE from rotational springs: 1/2*ktheta*(theta-alpha)^2
PE_A = (1/2)*rest(6)*(ThetaA-rest(1)).^2; PE_B = (1/2)*rest(7)*(ThetaB-rest(2)).^2;
PE_C = (1/2)*rest(8)*(ThetaC-rest(3)).^2; PE_D = (1/2)*rest(9)*(ThetaD-rest(4)).^2;
PE_E = (1/2)*rest(10)*(ThetaE-rest(5)).^2; 
PE_1 = (1/2)*rest(13)*(L1-rest(11)).^2; PE_2 = (1/2)*rest(14)*(L2-rest(12)).^2; 
PE_S = PE_A + PE_B + PE_C + PE_D + PE_E + PE_1 + PE_2;

% Sum all PE contributions
PE_T = PE_G + PE_S;

% Calculate PE from rotational springs: 1/2*ktheta*(theta-alpha)^2
pathPE_A = (1/2)*rest(6)*(pathThetaA-rest(1)).^2; pathPE_B = (1/2)*rest(7)*(pathThetaB-rest(2)).^2;
pathPE_C = (1/2)*rest(8)*(pathThetaC-rest(3)).^2; pathPE_D = (1/2)*rest(9)*(pathThetaD-rest(4)).^2;
pathPE_E = (1/2)*rest(10)*(pathThetaE-rest(5)).^2; 
pathPE_1 = (1/2)*rest(13)*(pathL1-rest(11)).^2;  pathPE_2 = (1/2)*rest(14)*(pathL2-rest(12)).^2;
pathPE_S = pathPE_A + pathPE_B + pathPE_C + pathPE_D + pathPE_E + pathPE_1 + pathPE_2;

% Sum all PE contributions
pathPE_T = pathPE_G + pathPE_S;

%%
[RMSE] = GetObjFun_Excavator_Path(rest,pathThetaA,pathThetaB,pathThetaC,pathThetaD,pathThetaE,pathL1,pathL2,pathPE_G);
[RMSEgrav] = GetObjFun_Excavator_Path(zeros(size(rest)),pathThetaA,pathThetaB,pathThetaC,pathThetaD,pathThetaE,pathL1,pathL2,pathPE_G);
normRMSE = RMSE/RMSEgrav;

%%
cmin = 1e04; cmax = 8e05; inc = 1e04;
figure('Position',[-1067,31,1055,801])
subplot(2,2,1)
contourf(PHI*180/pi,PSI*180/pi,PE_G,cmin:inc:cmax)
hold on
plot(pathPHI*180/pi,pathPSI*180/pi,'w--')
xlabel('\phi')
ylabel('\psi')
caxis([cmin cmax])
colormap jet
colorbar

subplot(2,2,2)
% [phi,psi] = meshgrid(PHI,PSI);
contourf(PHI*180/pi,PSI*180/pi,PE_T,cmin:inc:cmax)
hold on
plot(pathPHI*180/pi,pathPSI*180/pi,'w--')
xlabel('\phi')
ylabel('\psi')
caxis([cmin cmax])
colormap jet
colorbar

subplot(2,2,3)
contourf(PHI*180/pi,PSI*180/pi,PE_S)
hold on
plot(pathPHI*180/pi,pathPSI*180/pi,'w--')
xlabel('\phi')
ylabel('\psi')
% caxis([cmin cmax])
colormap jet
colorbar

subplot(2,2,4)
[desX,desY] = GetQuiver_Excavator_1(rest,PHI,PSI,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,L1,PE_G);
% [desX,desY] = GetQuiver_Excavator_2(rest,PHI,PSI,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,L1,PE_G);
quiver(PHI(1:10:end)*180/pi,PSI(1:10:end)*180/pi,desX(1:10:end,1:10:end),desY(1:10:end,1:10:end))
hold on
[Fx,Fy] = gradient(PE_T);
quiver(PHI(1:10:end)*180/pi,PSI(1:10:end)*180/pi,Fx(1:10:end,1:10:end),Fy(1:10:end,1:10:end))

%%
figure()
surf(PHI*180/pi,PSI*180/pi,PE_G)
hold on
contour3(PHI*180/pi,PSI*180/pi,PE_G,cmin:inc:cmax,'k')
surf(PHI*180/pi,PSI*180/pi,PE_T)
contour3(PHI*180/pi,PSI*180/pi,PE_T,cmin:inc:cmax,'k')
plot3(PHI*180/pi,PSI*180/pi,pathPE_T,'w--')
shading interp
xlabel('\phi')
ylabel('\psi')
caxis([cmin cmax])
zlim([cmin cmax])
colormap jet
colorbar
view([-135 35])

%%
figure()
area(pathPE_A + pathPE_B + pathPE_C + pathPE_D + pathPE_E + pathPE_1 + pathPE_2 + pathPE_G,'FaceColor',[0.5 0.5 0.5])
hold on
area(pathPE_A + pathPE_B + pathPE_C + pathPE_D + pathPE_E + pathPE_1 + pathPE_2,'FaceColor',[0.5 0.5 1])
area(pathPE_A + pathPE_B + pathPE_C + pathPE_D + pathPE_E + pathPE_1,'FaceColor','c')
area(pathPE_A + pathPE_B + pathPE_C + pathPE_D + pathPE_E,'FaceColor','m')
area(pathPE_A + pathPE_B + pathPE_C + pathPE_D,'FaceColor','b')
area(pathPE_A + pathPE_B + pathPE_C,'FaceColor','g')
area(pathPE_A + pathPE_B,'FaceColor','y')
area(pathPE_A,'FaceColor','r')
xlim([1 121])
ylim([0 8e05])

%%
perpPE_A = diag(flipud(PE_A)); perpPE_B = diag(flipud(PE_B)); perpPE_C = diag(flipud(PE_C)); perpPE_D = diag(flipud(PE_D)); 
perpPE_E = diag(flipud(PE_E)); perpPE_1 = diag(flipud(PE_1)); perpPE_2 = diag(flipud(PE_2)); perpPE_G = diag(flipud(PE_G));

figure()
area(perpPE_A + perpPE_B + perpPE_C + perpPE_D + perpPE_E + perpPE_1 + perpPE_2 + perpPE_G,'FaceColor',[0.5 0.5 0.5])
hold on
area(perpPE_A + perpPE_B + perpPE_C + perpPE_D + perpPE_E + perpPE_1 + perpPE_2,'FaceColor',[0.5 0.5 1])
area(perpPE_A + perpPE_B + perpPE_C + perpPE_D + perpPE_E + perpPE_1,'FaceColor','c')
area(perpPE_A + perpPE_B + perpPE_C + perpPE_D + perpPE_E,'FaceColor','m')
area(perpPE_A + perpPE_B + perpPE_C + perpPE_D,'FaceColor','b')
area(perpPE_A + perpPE_B + perpPE_C,'FaceColor','g')
area(perpPE_A + perpPE_B,'FaceColor','y')
area(perpPE_A,'FaceColor','r')
xline(60)
xlim([1 121])
ylim([0 8e05])






