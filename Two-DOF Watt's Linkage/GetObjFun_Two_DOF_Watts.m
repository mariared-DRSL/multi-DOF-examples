function [ObjFun] = GetObjFun_Two_DOF_Watts(r,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,PE_G)

    PE = (1/2)*r(6)*(ThetaA-r(1)).^2 + (1/2)*r(7)*(ThetaB-r(2)).^2 + (1/2)*r(8)*(ThetaC-r(3)).^2 + (1/2)*r(9)*(ThetaD-r(4)).^2 + (1/2)*r(10)*(ThetaE-r(5)).^2 + PE_G;

    PEavg = mean(mean(PE(~isnan(PE))));
    
    N=0;
    for i = 1:size(PE,1)
        for j = 1:size(PE,2)
            if ~isnan(PE(i,j))
                N=N+1;
                dev = PE(i,j)-PEavg;
                dev2(N) = dev^2;
            end
        end
    end
    
    sumup = sum(dev2);
    RMS = sqrt(sumup/N);
    ObjFun = RMS;

end