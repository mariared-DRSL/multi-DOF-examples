function [ObjFun] = GetObjFun_Excavator_Path(r,pathThetaA,pathThetaB,pathThetaC,pathThetaD,pathThetaE,pathL1,pathL2,pathPE_G)

    PEpath = (1/2)*r(6)*(pathThetaA-r(1)).^2 + (1/2)*r(7)*(pathThetaB-(r(2))).^2 + (1/2)*r(8)*(pathThetaC-(r(3))).^2 + (1/2)*r(9)*(pathThetaD-(r(4))).^2 +...
         (1/2)*r(10)*(pathThetaE-(r(5))).^2 + (1/2)*r(13)*(pathL1-r(11)).^2 +  (1/2)*r(14)*(pathL2-r(12)).^2 + pathPE_G;
    
    PEavg = mean(mean(mean(PEpath)));
    dev = PEpath-PEavg;
    dev2 = dev.^2;
    RMS = sqrt(mean(mean(mean(dev2))));
    ObjFun = RMS;
    
end