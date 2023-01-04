%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%--------- Two-DOF Watt's Linkage: Continuous Equilibrium -------%%%
%%%------------------------- Curved Path --------------------------%%%
%  "Programming Continuous Equilibrium Motions in Multi-DOF Systems" %
%%% ------------------ Redoutey, M., Filipov, E.T. ----------------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

%% Load two-DOF Watt's linkage kinemtics
load('SpringKinematics','Phi1','Phi2','SpringKinematics','PE_G')

Phi1 = Phi1(1:5:end); Phi2 = Phi2(1:5:end);
PE_G = PE_G(1:5:end,1:5:end);

PathPhi1 = Phi1(1:end);
PathPhi2 = (PathPhi1-3.7).^2+1.5;

for m = 1:length(PathPhi1)
    hi = abs(Phi1-PathPhi1(m));
    i = find(hi==min(hi));
    bi = abs(Phi2-PathPhi2(m));
    j = find(bi==min(bi));
    
    pathPE_G(m) = PE_G(j,i);
    myfun = @(SpringKinematics)SpringKinematics(j,i);
    pathExtSpr(:,m) = structfun(myfun,SpringKinematics,'UniformOutput',true); %each of the springs' lengths along the path. each row is a spring
end

SpringIDs = fieldnames(SpringKinematics);

%% Compute rest angles and stiffnesses that optimize for continuous equilibrium
[Stiffness,RestPos,fval,history,exitflag,output,lambda] = FindRest_CurvedPath(pathExtSpr,SpringKinematics,SpringIDs,Phi1,Phi2,PathPhi1,PathPhi2,pathPE_G,PE_G);


%% Calculate PE from torsional spring j: 1/2*k_j*(theta_j-alpha_j)^2
PEsum = 0;
for m = 1:length(SpringIDs)

    k = Stiffness(m); L0a0 = RestPos(m);
    PE_S.(SpringIDs{m}) = (1/2).*k.*(SpringKinematics.(SpringIDs{m})-L0a0).^2;
    PEsum = PEsum + (1/2).*k.*(SpringKinematics.(SpringIDs{m})-L0a0).^2;

end

% Sum all PE contributions
PE_T = PEsum + PE_G;

%% Get potential energy along path 
for m = 1:length(PathPhi1)
    hi = abs(Phi1-PathPhi1(m));
    i = find(hi==min(hi));
    bi = abs(Phi2-PathPhi2(m));
    j = find(bi==min(bi));
    
    pathPE_T(m) = PE_T(j,i);

end

%% Calculate RMSE
PEavg = mean(pathPE_T);
RMSE = sqrt(sum((pathPE_T-PEavg).^2/length(pathPE_T)));
RMSEgrav = sqrt(sum((pathPE_G-mean(pathPE_G)).^2/length(pathPE_G)));
normRMSE = RMSE/RMSEgrav;

%% Plot PE_T
figure('Position',[-799,371,560,420])
contourf(Phi1,Phi2,PE_T,10)
hold on
plot(PathPhi1,PathPhi2,'w-*')
colormap jet
colorbar
title('PE_T')

%% Plot gradient of PE_T
GetQuiver(Phi1,PathPhi1,Phi2,PathPhi2,PE_T)

%% Plot potential energy surface
figure()
surf(Phi1*180/pi,Phi2*180/pi,PE_T,'FaceAlpha',0.75)
hold on
contour3(Phi1*180/pi,Phi2*180/pi,PE_T,0:1:4000,'k')
shading interp
surf(Phi1*180/pi,Phi2*180/pi,PE_G,'FaceAlpha',0.75)
contour3(Phi1*180/pi,Phi2*180/pi,PE_G,0:1:4000,'k')
caxis([100 300]);
zlim([100 300]);
plot3(PathPhi1*180/pi,PathPhi2*180/pi,pathPE_T,'k-*')
colormap jet
view([-115 30])
shading interp
colorbar
xlabel('\phi_1')
ylabel('\phi_2')
