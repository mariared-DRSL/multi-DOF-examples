

function [L0_LR,L0_TB,alpha_Dome,k_LR,k_TB,k_Dome,history,exitflag,output,lambda] = FindRest_Dome(Length_LR_Path,Length_TB_Path,Phi_Dome_Path,PE_G_Path,Length_LR,Length_TB,Phi_Dome,PE_G)

% Design variables: [L0_LR,L0_TB,k_LR,k_TB]

% Linear constraints
A = []; b = []; Aeq = []; beq = []; 

% Lower and upper bounds for design variables: L0_LR, L0_TB, alpha_dome, k_LR, k_TB, k_dome
lb = [0,0,0,0,0,0]; ub = [max(max(Length_LR)),max(max(Length_TB)),2*pi,200,200,200];
% lb = [0,0,0,0,0,0,0,0,0,0]; ub = [5,5,Inf,5,2*pi,500,500,500,500,500];

% Initial values for design variables 
r0 = [1,1,pi,10,10,10]; 
r0 = [0.6547 2.3023 1.6866 105.9255 197.4521 198.9299];

% Objective function: sum(abs(diff(PE)))
[ObjFun] = @(r)GetObjFun_Dome_CoupledSprings_SS_2(r,Length_LR_Path,Length_TB_Path,Phi_Dome_Path,PE_G_Path);

[nonlcon] = @(r)GetCon_Dome_CoupledSprings_SS_2(r,Length_LR,Length_TB,Phi_Dome,PE_G);

% Call Optimization
[rest,~,exitflag,output,lambda,history,~] = runfmincon(ObjFun,r0,A,b,Aeq,beq,lb,ub,nonlcon);
L0_LR = rest(1); L0_TB = rest(2); alpha_Dome = rest(3); 
k_LR = rest(4); k_TB = rest(5); k_Dome = rest(6);

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
options = optimoptions(@fmincon,'OutputFcn',@outfun,'Display','iter','Algorithm','interior-point','ConstraintTolerance',0.05,'FiniteDifferenceStepSize',1e-6,'MaxFunctionEvaluations',6000,'MaxIterations',2000,'StepTolerance',1e-6,'PlotFcn',{'optimplotx','optimplotfirstorderopt','optimplotfvalconstr'}); % 'PlotFcn','optimplotfval'   'PlotFcn','optimplotx'; 
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

