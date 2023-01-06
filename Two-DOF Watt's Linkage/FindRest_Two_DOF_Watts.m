

function [aA,aB,aC,aD,aE,kA,kB,kC,kD,kE,fval,flag,output,allmins] = FindRest_Two_DOF_Watts(ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,PE_G) 

% Design variables: [alphaA, alphaB, alphaC, alphaD, kthetaA, kthetaB, kthetaC, kthetaD]

% Linear constraints
A = []; b = []; Aeq = []; beq = []; 

% Lower and upper bounds for design variables
lb = [0,0,0,0,0,0,0,0,0,0]; ub = [2*pi,2*pi,2*pi,2*pi,2*pi,Inf,Inf,Inf,Inf,Inf];

% Initial values for design variables 
r0 = [pi,pi,pi,pi,pi,10,10,10,10,10]; 

% Objective function: sum(abs(diff(PE)))
[ObjFun] = @(r)GetObjFun_Two_DOF_Watts(r,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,PE_G);

rng default % For reproducibility
opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter','MaxFunctionEvaluations',6000); % 'PlotFcn','optimplotfval'   'PlotFcn','optimplotx'
problem = createOptimProblem('fmincon','objective',ObjFun,'x0',r0,'lb',lb,'ub',ub,'options',opts);
% ms = MultiStart;
% [rest,fval] = run(ms,problem,20);
gs = GlobalSearch('NumStageOnePoints',500);
[rest,fval,flag,output,allmins] = run(gs,problem);

aA = rest(1); aB = rest(2); aC = rest(3); aD = rest(4); aE = rest(5);
kA = rest(6); kB = rest(7); kC = rest(8); kD = rest(9); kE = rest(10);

end 
