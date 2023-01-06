function [ObjFun] = GetObjFun_RMS_Linear_Phi1_Phi2_Path(r,PHI1,PHI2,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,PE_G)

    PE = (1/2)*r(6)*(ThetaA-r(1)).^2 + (1/2)*r(7)*(ThetaE-r(2)).^2 + (1/2)*r(8)*(ThetaB-r(3)).^2 + (1/2)*r(9)*(ThetaC-r(4)).^2 + (1/2)*r(10)*(ThetaD-r(5)).^2 + PE_G;

    idx = 0; dt = 1/(size(PHI1,1)-1);
    for t = 0:dt:1
        idx = idx+1;
        
        u = (max(max(PHI1))-min(min(PHI1)))*t + min(min(PHI1));
        v = max(max(PHI2)) - (max(max(PHI2))-min(min(PHI2)))*t;

        i = find(round(PHI1(1,:),4)==round(u,4));
        j = find(round(PHI2(:,1),4)==round(v,4));
        
        pathPE(idx) = PE(j,i);         
    end  
    
    PEavg_path = mean(pathPE);
    
    N=0;
    for i = 1:length(pathPE)
        N=N+1;
        dev = pathPE(i)-PEavg_path;
        dev2(N) = dev^2;
    end
    
    sumup = sum(dev2);
    RMS = sqrt(sumup/N);
    ObjFun = RMS;

end