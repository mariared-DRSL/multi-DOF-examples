%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%--------------- Miura-ori Cell: Zero Stiffness -----------------%%%
%%%---------------------- Maria Redoutey --------------------------%%%
%%%---------------------- 10 March 2021 ---------------------------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [as,ks,history,exitflag,lambda] = FindRest_One_DOF_Watts(ThetaA,ThetaB,ThetaC,ThetaD,PE_G) 

% Design variables: [alphaA, alphaB, alphaC, alphaD, kthetaA, kthetaB, kthetaC, kthetaD]

% Linear constraints
A = []; b = []; Aeq = []; beq = []; 

% Lower and upper bounds for design variables
lb = [0 0 0 0 0 0 0 0]; ub = [2*pi 2*pi 2*pi 2*pi 5 5 5 5];

% Initial values for design variables
r0 = [pi pi pi pi 1 1 1 1]; 

% Objective function: sum(abs(diff(PE)))
[ObjFun] = @(r)GetObjFun_RMS_One_DOF_Watts(r,ThetaA,ThetaB,ThetaC,ThetaD,PE_G);

% Call Optimization
[rest,~,exitflag,lambda,history,~] = runfmincon(ObjFun,r0,A,b,Aeq,beq,lb,ub);
aA = rest(1); aB = rest(2); aC = rest(3); aD = rest(4);
kA = rest(5); kB = rest(6); kC = rest(7); kD = rest(8);
as = [aA aB aC aD]; ks = [kA kB kC kD];

end 

function [rest,fval,exitflag,lambda,history,searchdir] = runfmincon(PEtot,r0,A,b,Aeq,beq,lb,ub)
 
history.x = [];
history.fval = [];
searchdir = [];

% history is a struct with two values: fval and x
% history.x is a n-by-6 matrix where n is the number of optimization steps
% The columns of history.x correspond to variables being optimized: 
% [alphaA, alphaB, alphaC, alphaD, LAC, LBD, kthetaA, kthetaB, kthetaC, kthetaD, kAC, kBD]
 
% Call optimization
options = optimoptions(@fmincon,'OutputFcn',@outfun,'Display','iter','Algorithm','interior-point','FiniteDifferenceStepSize',1e-6,'MaxFunctionEvaluations',12000,'MaxIterations',2000,'StepTolerance',1e-6); 
[rest,fval,exitflag,~,lambda] = fmincon(PEtot,r0,A,b,Aeq,beq,lb,ub,[],options);
 
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

function [ObjFun] = GetObjFun_RMS_One_DOF_Watts(r,ThetaA,ThetaB,ThetaC,ThetaD,PE_G)

    PE = (1/2)*r(5)*(ThetaA-r(1)).^2 + (1/2)*r(6)*(ThetaB-r(2)).^2 + (1/2)*r(7)*(ThetaC-r(3)).^2 + (1/2)*r(8)*(ThetaD-r(4)).^2 + PE_G;

    PEavg = mean(PE(~isnan(PE)));
    
    N = length(PE);

    RMS = sqrt(sum((PE-PEavg).^2)/N);
    
    ObjFun = RMS;

end