function [c,ceq] = GetCon_Excavator_1(r,PHI,PSI,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,L1,L2,PE_G)

    
    PE = (1/2)*r(6)*(ThetaA-r(1)).^2 + (1/2)*r(7)*(ThetaB-(r(2))).^2 + (1/2)*r(8)*(ThetaC-(r(3))).^2 + (1/2)*r(9)*(ThetaD-(r(4))).^2 +...
         (1/2)*r(10)*(ThetaE-(r(5))).^2 + (1/2)*r(13)*(L1-r(11)).^2 + (1/2)*r(14)*(L2-r(12)).^2 + PE_G;
     
    unitvect1 = [-1 1]/sqrt(2);
    unitvect2 = [1 -1]/sqrt(2);
    
    [Fx,Fy] = gradient(PE);
    
    for j = 1:size(PE,1)
        for i = 1:size(PE,2)
            
            gradvect = [Fx(j,i) Fy(j,i)]/norm([Fx(j,i) Fy(j,i)]);
            
            if j > i
                desX(j,i) = unitvect1(1);
                desY(j,i) = unitvect1(2);
                c(j,i) = acos(dot(gradvect,unitvect1)) - 3*pi/8;
            elseif j < i
                desX(j,i) = unitvect2(1);
                desY(j,i) = unitvect2(2);
                c(j,i) = acos(dot(gradvect,unitvect2)) - 3*pi/8;
            end
            
        end
    end
    
    ceq = [];

end 