close all
clear 
clc

set(0, 'defaultFigureRenderer', 'painters')

[Phi_Assemble_Vector,Phi_Dome_Vector,Length_LR,Length_TB,PE_G,Length_Span,Length_Rise] = GetGeo_Dome();

Phi_Dome = repmat(Phi_Dome_Vector,36,1)';

% Path: Phi_Assemble = end, Phi_Dome = 0:end
PE_G_Path = PE_G(:,end);
Length_LR_Path = Length_LR(:,end);
Length_TB_Path = Length_TB(:,end);
% Length_Span_Path = Length_Span(:,end);
% Length_Rise_Path = Length_Rise(:,end);
Phi_Dome_Path = Phi_Dome(:,end);

%%
[L0_LR,L0_TB,alpha_Dome,k_LR,k_TB,k_Dome,history,exitflag,output,lambda] = FindRest_Dome(Length_LR_Path,Length_TB_Path,Phi_Dome_Path,PE_G_Path,Length_LR,Length_TB,Phi_Dome,PE_G);

PE_LR = 7*(1/2)*k_LR*(Length_LR-L0_LR).^2;
PE_TB = 7*(1/2)*k_TB*(Length_TB-L0_TB).^2;
% PE_Span = (1/2)*k_Span*(Length_Span-L0_Span).^2;
% PE_Rise = (1/2)*k_Rise*(Length_Rise-L0_Rise).^2;
PE_Dome = 7*(1/2)*k_Dome*(Phi_Dome-alpha_Dome).^2;

PE_T = PE_LR + PE_TB + PE_Dome + PE_G;

PE_LR_Path = 7*(1/2)*k_LR*(Length_LR_Path-L0_LR).^2;
PE_TB_Path = 7*(1/2)*k_TB*(Length_TB_Path-L0_TB).^2;
% PE_Span_Path = (1/2)*k_Span*(Length_Span_Path-L0_Span).^2;
% PE_Rise_Path = (1/2)*k_Rise*(Length_Rise_Path-L0_Rise).^2;
PE_Dome_Path = 7*(1/2)*k_Dome*(Phi_Dome_Path-alpha_Dome).^2;

PE_T_Path = PE_LR_Path + PE_TB_Path + PE_Dome_Path + PE_G_Path;

%%
% PEavg = mean(PE_T_Path);
% PE_Gavg = mean(PE_G_Path);
% 
% N=0;
% for i = 1:length(PE_T_Path)
%     N=N+1;
%     dev = PE_T_Path(i)-PEavg;
%     dev2(N) = dev^2;
%     dev2G(N) = (PE_G_Path(i)-PE_Gavg)^2;
% end
%     
% sumup = sum(dev2);
% RMSE = sqrt(sumup/N);
% RMSEgrav = sqrt(sum(dev2G)/N);
% normRMSE = RMSE/RMSEgrav;

PEavg_Path = mean(PE_T_Path);
PE_Gavg_Path = mean(PE_G_Path);

N=0;
for i = 1:size(PE_T_Path,1)
    for j = 1:size(PE_T_Path,2)
        if ~isnan(PE_T_Path(i,j))
            N=N+1;
            dev = PE_T_Path(i,j)-PEavg_Path;
            dev2(N) = dev^2;
            devG = PE_G_Path(i,j)-PE_Gavg_Path;
            dev2G(N) = devG^2;
        end
    end
end

sumup = sum(dev2);
RMSE = sqrt(sumup/N);
RMSEgrav = sqrt(sum(dev2G)/N);
normRMSE = RMSE/RMSEgrav;

%%
figure('Position',[-800,350,560,420])
contourf(Phi_Assemble_Vector*180/pi,Phi_Dome_Vector*180/pi,PE_G,[0:100:2e5])
colorbar
xlabel('\phi_a')
ylabel('\phi_d')

figure('Position',[-800,850,560,420])
contourf(Phi_Assemble_Vector*180/pi,Phi_Dome_Vector*180/pi,PE_T,[0:100:2e5])
colorbar
xlabel('\phi_a')
ylabel('\phi_d')

%%
figure('Position',[-800,350,560,420])
surf(Phi_Assemble_Vector*180/pi,Phi_Dome_Vector*180/pi,PE_T)
hold on
contour3(Phi_Assemble_Vector*180/pi,Phi_Dome_Vector*180/pi,PE_T,0:100:2e5,'k')
caxis([0 16000]);
zlim([0 16000]);
shading interp
surf(Phi_Assemble_Vector*180/pi,Phi_Dome_Vector*180/pi,PE_G)
contour3(Phi_Assemble_Vector*180/pi,Phi_Dome_Vector*180/pi,PE_G,0:100:2e5,'k')
caxis([0 16000]);
zlim([0 16000]);
colormap jet
view([35 35])
shading interp
colorbar
xlabel('\phi_a')
ylabel('\phi_d')

%%
GetQuiver_SS(Phi_Assemble_Vector,Phi_Dome_Vector,PE_T)
axis equal

%%
seg = 1:35;
figure('Position',[-800,850,560,420])
area([PE_T(1,seg)'; PE_LR_Path+PE_TB_Path+PE_Dome_Path+PE_G_Path],'FaceColor',[0.5 0.5 0.5])
hold on
area([(PE_LR(1,seg)+PE_TB(1,seg)+PE_Dome(1,seg))'; PE_LR_Path+PE_TB_Path+PE_Dome_Path],'FaceColor','y')
area([(PE_LR(1,seg)+PE_TB(1,seg))'; PE_LR_Path+PE_TB_Path],'FaceColor','m')
area([(PE_LR(1,seg))'; PE_LR_Path],'FaceColor','c')
xlim([1 70]);
