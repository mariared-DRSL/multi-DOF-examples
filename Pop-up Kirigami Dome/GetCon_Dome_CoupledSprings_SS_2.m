function [c,ceq] = GetCon_Dome_CoupledSprings_SS_2(r,Length_LR,Length_TB,Phi_Dome,PE_G)

    PE = 7*(1/2)*r(4)*(Length_LR-r(1)).^2 + 7*(1/2)*r(5)*(Length_TB-r(2)).^2 + 7*(1/2)*r(6)*(Phi_Dome-r(3)).^2 + PE_G;
       
    desX = ones(size(PE))*-1;
    desY = zeros(size(PE));
    
    [Fx,Fy] = gradient(PE);
    
    for ii = 1:size(Fx,2)
        for jj = 1:size(Fx,1)
            gradvect = [Fx(jj,ii) Fy(jj,ii)]./(sqrt(Fx(jj,ii).^2+Fy(jj,ii).^2));
            unitvect = [desX(jj,ii) desY(jj,ii)]./(sqrt(desX(jj,ii).^2+desY(jj,ii).^2));
            c(jj,ii) = acos(dot(gradvect,unitvect)) - pi/8;
            ceq = [];
        end
    end
    c(end+1,:) = [-diff(diff(PE(1,:))) max(max(PE))-max(max(PE_G))*3.5 -1];

%     c =  (PE(:,end)' - min(PE')) - 0.1; % dont use
%     ceq = []; % dont use

%     c =  [diff(PE(1,:)) -diff(diff(PE(1,:)))];
%     ceq = min(PE') - PE(:,end)'; 
    
end