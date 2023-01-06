function [ObjFun] = GetObjFun_Dome_CoupledSprings_SS_2(r,Length_LR_Path,Length_TB_Path,Phi_Dome_Path,PE_G_Path)

    PE_Path = 7*(1/2)*r(4)*(Length_LR_Path-r(1)).^2 + 7*(1/2)*r(5)*(Length_TB_Path-r(2)).^2 + (7*(1/2)*r(6)*(Phi_Dome_Path-r(3)).^2) + PE_G_Path;

    PEavg_Path = mean(PE_Path);
    
    N=0;
    for i = 1:size(PE_Path,1)
        for j = 1:size(PE_Path,2)
            if ~isnan(PE_Path(i,j))
                N=N+1;
                dev = PE_Path(i,j)-PEavg_Path;
                dev2(N) = dev^2;
            end
        end
    end
    
    sumup = sum(dev2);
    RMS = sqrt(sumup/N);
    ObjFun = RMS;
    
end