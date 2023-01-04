

function [aA,aB,aC,aD,aE,aF,kA,kB,kC,kD,kE,kF,history,exitflag,output,lambda] = FindRest_Three_DOF_Watts(ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,ThetaF,PE_G)

% Design variables: [alphaA, alphaB, alphaC, alphaD, kthetaA, kthetaB, kthetaC, kthetaD]

% Linear constraints
A = []; b = []; Aeq = []; beq = []; 

% Lower and upper bounds for design variables
lb = [0,0,0,0,0,0,0,0,0,0,0,0]; ub = [2*pi,2*pi,2*pi,2*pi,2*pi,2*pi,Inf,Inf,Inf,Inf,Inf,Inf];

% Initial values for design variables 
r0 = [(ub(1:6)-lb(1:6))/2+lb(1:6),10,10,10,10,10,10]; 

% Objective function
ObjFun = @(r)GetObjFun_RMS_Three_DOF_Watts(r,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,ThetaF,PE_G);

% Call Optimization
[rest,~,exitflag,output,lambda,history,~] = runfmincon(ObjFun,r0,A,b,Aeq,beq,lb,ub);
aA = rest(1); aB = rest(2); aC = rest(3); aD = rest(4); aE = rest(5); aF = rest(6);
kA = rest(7); kB = rest(8); kC = rest(9); kD = rest(10); kE = rest(11); kF = rest(12);

end 

function [rest,fval,exitflag,output,lambda,history,searchdir] = runfmincon(PEtot,r0,A,b,Aeq,beq,lb,ub)
 
history.x = [];
history.fval = [];
searchdir = [];

% history is a struct with two values: fval and x
% history.x is a n-by-6 matrix where n is the number of optimization steps
% The columns of history.x correspond to variables being optimized: 
% [alphaA, alphaB, alphaC, alphaD, LAC, LBD, kthetaA, kthetaB, kthetaC, kthetaD, kAC, kBD]
 
% Call optimization
options = optimoptions(@fmincon,'OutputFcn',@outfun,'Display','iter','Algorithm','interior-point','FiniteDifferenceStepSize',1e-6,'MaxFunctionEvaluations',6000,'MaxIterations',2000,'StepTolerance',1e-6); 
[rest,fval,exitflag,output,lambda] = fmincon(PEtot,r0,A,b,Aeq,beq,lb,ub,[],options);
 
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

function [ObjFun] = GetObjFun_RMS_Three_DOF_Watts(r,ThetaA,ThetaB,ThetaC,ThetaG,ThetaJ,ThetaK,PE_G)

    PE = (1/2)*r(7)*(ThetaA-r(1)).^2 + (1/2)*r(8)*(ThetaB-r(2)).^2 + (1/2)*r(9)*(ThetaC-r(3)).^2 + (1/2)*r(10)*(ThetaG-r(4)).^2 + (1/2)*r(11)*(ThetaJ-r(5)).^2 + (1/2)*r(12)*(ThetaK-r(6)).^2 + PE_G;

    for j = size(PE,1):-1:1
        for k = size(PE,2):-1:1
            for i = size(PE,3):-1:1
                if isnan(PE(j,k,i))
                    PE(j,k,i) = PE(j,k,i+1);
                end
            end
        end
    end
    

    PEavg = mean(mean(mean(PE)));
    dev = PE-PEavg;
    dev2 = dev.^2;
    RMS = sqrt(mean(mean(mean(dev2))));
    ObjFun = RMS;

end