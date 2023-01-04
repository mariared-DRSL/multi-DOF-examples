

function [rest,history,exitflag,output,lambda] = FindRest_Phi2_90_Path(PHI1,PHI2,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,PE_G) 

% Design variables: [alphaA, alphaB, alphaC, alphaD, kthetaA, kthetaB, kthetaC, kthetaD]

% Linear constraints
A = []; b = []; Aeq = []; beq = []; 

% Lower and upper bounds for design variables
lb = [0,0,0,0,0,0,0,0,0,0]; ub = [2*pi,2*pi,2*pi,2*pi,2*pi,200,200,200,200,200];

% Initial values for design variables 
r0 = [pi,pi,pi,pi,pi,10,10,10,10,10]; 

% Objective function: sum(abs(diff(PE)))
ObjFun = @(r)GetObjFun_Phi2_90_Path(r,PHI1,PHI2,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,PE_G);

nonlcon = @(r)mycon(r,PHI1,PHI2,ThetaA,ThetaE,ThetaB,ThetaC,ThetaD,PE_G);

% Call Optimization
[rest,~,exitflag,output,lambda,history,~] = runfmincon(ObjFun,r0,A,b,Aeq,beq,lb,ub,nonlcon);
% aA = rest(1); aE = rest(2); aB = rest(3); aC = rest(4); aD = rest(5);
% kA = rest(6); kE = rest(7); kB = rest(8); kC = rest(9); kD = rest(10);

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
options = optimoptions(@fmincon,'OutputFcn',@outfun,'Display','iter','Algorithm','interior-point','FiniteDifferenceStepSize',1e-3,'MaxFunctionEvaluations',12000,'MaxIterations',2000,'StepTolerance',1e-12,'ConstraintTolerance',0.05); 
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

function [c,ceq] = mycon(r,PHI1,PHI2,ThetaA,ThetaE,ThetaB,ThetaC,ThetaD,PE_G)

    PE = (1/2)*r(6)*(ThetaA-r(1)).^2 + (1/2)*r(7)*(ThetaE-r(2)).^2 + (1/2)*r(8)*(ThetaB-r(3)).^2 + (1/2)*r(9)*(ThetaC-r(4)).^2 + (1/2)*r(10)*(ThetaD-r(5)).^2 + PE_G;

    idx = 0; dt = 1/(size(PHI1,1)-1);
    for t = 0:dt:1
        idx = idx+1;
        
        u = (max(max(PHI1))-min(min(PHI1)))*t + min(min(PHI1));
        v = 1*pi/2;

        i = find(round(PHI1(1,:),4)==round(u,4));
        j = find(round(PHI2(:,1),4)==round(v,4));
        
        pathPE(idx) = PE(j,i);
        pathPE_1(idx) = (1/2)*r(6)*(PHI1(j,i)-r(1)).^2;
        
        unitvect1 = [0, -1];
        unitvect2 = [0, 1];
        
        [Fx,Fy] = gradient(PE);
        gradvect1 = [Fx(j-1,i) Fy(j-1,i)]./sqrt(Fx(j-1,i)^2+Fy(j-1,i)^2);
        gradvect2 = [Fx(j+1,i) Fy(j+1,i)]./sqrt(Fx(j+1,i)^2+Fy(j+1,i)^2);
        
        dotprod1(idx) = dot(gradvect1,unitvect1);
        dotprod2(idx) = dot(gradvect2,unitvect2);
         
    end
    
    p1Phi1 = repmat(pi,1,41); p1Phi2 = [130:-1:90]*pi/180;
    p2Phi1 = [180:220]*pi/180; p2Phi2 = repmat(pi/2,1,41);

    for k = 1:length(p1Phi1)
        i = find(round(PHI1(1,:),4)==round(p1Phi1(k),4));
        j = find(round(PHI2(:,1),4)==round(p1Phi2(k),4));
        p1PE_1(k) = (1/2)*r(6)*(PHI1(j,i)-r(1)).^2;
    end 

    
    c = [];
    ceq = [1-dotprod1 1-dotprod2 p1PE_1]; 

end


function [ObjFun] = GetObjFun_Phi2_90_Path(r,PHI1,PHI2,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,PE_G)

    PE = (1/2)*r(6)*(ThetaA-r(1)).^2 + (1/2)*r(7)*(ThetaE-r(2)).^2 + (1/2)*r(8)*(ThetaB-r(3)).^2 + (1/2)*r(9)*(ThetaC-r(4)).^2 + (1/2)*r(10)*(ThetaD-r(5)).^2 + PE_G;

    idx = 0; dt = 1/(size(PHI1,1)-1);
    for t = 0:dt:1
        idx = idx+1;
        
        u = (max(max(PHI1))-min(min(PHI1)))*t + min(min(PHI1));
        v = pi/2;
        
        i = find(round(PHI2(:,1),4)==round(v,4));
        j = find(round(PHI1(1,:),4)==round(u,4));

        pathPE(idx) = PE(i,j);%%%%%%%%%%%%%%%%%%%%%%%%%check i,j
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