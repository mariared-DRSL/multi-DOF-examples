

function [rest,history,exitflag,output,lambda] = FindRest_Linear_Phi1_Phi2_Path(PHI1,PHI2,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,PE_G)
% Design variables: [alphaA, alphaB, alphaC, alphaD, kthetaA, kthetaB, kthetaC, kthetaD]

% Linear constraints
A = []; b = []; Aeq = []; beq = []; 

% Lower and upper bounds for design variables
lb = [0 0 0 0 0 0 0 0 0 0]; ub = [2*pi 2*pi 2*pi 2*pi 2*pi Inf Inf Inf Inf Inf];

% Initial values for design variables 
r0 = [pi/2 pi/2 pi/2 pi/2 pi/2 10 10 10 10 10];
r0 = [0.0002    1.9195    2.9296    0.3204    1.5125    9.0699   14.2432   10.9270    4.4595   15.7967];
r0 = [0.0000    2.3195    2.8928    2.0512    0.9913    8.0865   19.0849   12.2743    5.9010   14.4386];
r0 = [0.0660 2.2436 3.0195 1.9651 0.1717 9.0420 20.3329 11.6816 13.9004 9.6246];
r0 = [0.0062 2.2214 3.1250 2.0454 0.6142 9.1055 19.2756 10.6611 12.2878 10.7136];
r0 = [0.0295    2.2367    2.5933    2.4435    1.0981    7.5373   18.9598   11.1772   11.2620   11.1475];
r0 = [0.0501 2.2178 2.6317 2.6205 1.3945 7.6992 18.0604 10.2306 9.8422 12.1723];

% Objective function: sum(abs(diff(PE)))
ObjFun = @(r)GetObjFun_RMS_Linear_Phi1_Phi2_Path(r,PHI1,PHI2,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,PE_G);

nonlcon = @(r)mycon(r,PHI1,PHI2,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,PE_G);

% Call Optimization
[rest,~,exitflag,output,lambda,history,~] = runfmincon(ObjFun,r0,A,b,Aeq,beq,lb,ub,nonlcon);
% a1 = rest(1); a2 = rest(2); a3 = rest(3); a4 = rest(4); a5 = rest(5);
% k1 = rest(6); k2 = rest(7); k3 = rest(8); k4 = rest(9); k5 = rest(10);

end 

function [rest,fval,exitflag,output,lambda,history,searchdir] = runfmincon(PEtot,r0,A,b,Aeq,beq,lb,ub,nonlcon)
 
history.x = [];
history.fval = [];
searchdir = [];

% history is a struct with two values: fval and x
% history.x is a n-by-6 matrix where n is the number of optimization steps
% The columns of history.x correspond to variables being optimized: 
% [alphaA, alphaB, alphaC, alphaD, LAC, LBD, kthetaA, kthetaB, kthetaC, kthetaD, kAC, kBD]
 
% Call optimization
options = optimoptions(@fmincon,'OutputFcn',@outfun,'Display','iter','Algorithm','interior-point','FiniteDifferenceStepSize',1e-6,...
    'MaxFunctionEvaluations',6000,'MaxIterations',2000,'StepTolerance',1e-6,'ConstraintTolerance',0.05,...
    'PlotFcn',{'optimplotx','optimplotfirstorderopt','optimplotfvalconstr'}); 
[rest,fval,exitflag,output,lambda] = fmincon(PEtot,r0,A,b,Aeq,beq,lb,ub,nonlcon,options);
 
 function stop = outfun(x,optimValues,state)
     stop = false;
     switch state
         case 'init'
        
         case 'iter'
                history.fval = [history.fval; optimValues.fval];
                history.x = [history.x; x];
         case 'done'
           
         otherwise
     end
 end

end

function [c,ceq] = mycon(r,PHI1,PHI2,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,PE_G)

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
    
    [Fx,Fy] = gradient(PE);
    
    unitvect1 = [-1/sqrt(2), -1/sqrt(2)];
    unitvect2 = [1/sqrt(2), 1/sqrt(2)];
    
    dotpro = ones(50,50);
    for i = 1:51 %20:40
        for j = 1:51 %20:40 %51-(i-10):-1:51-(i+10) %20:40 %
            gradvect = [Fx(j,i) Fy(j,i)]./sqrt(Fx(j,i)^2+Fy(j,i)^2);
            if j < 51 - i + 1
                dotpro(j,i) = dot(gradvect,unitvect1); 
                c(j,i) = acos(dot(gradvect,unitvect1)) - pi/2;
                fx(j,i) = unitvect1(1); fy(j,i) = unitvect1(2);
            elseif j > 51 - i + 1
                dotpro(j,i) = dot(gradvect,unitvect2); 
                c(j,i) = acos(dot(gradvect,unitvect2)) - pi/2;
                fx(j,i) = unitvect2(1); fy(j,i) = unitvect2(2);
            end
        end
    end
           
    dotpro = dotpro(2:end,2:end); 
    %c = 1-dotpro; 
    ceq=[];
    
end
