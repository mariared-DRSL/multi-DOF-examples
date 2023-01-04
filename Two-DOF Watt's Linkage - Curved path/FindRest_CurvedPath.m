

function [Stiffness,RestPos,fval,history,exitflag,output,lambda] = FindRest_CurvedPath(pathExtSpr,SpringKinematics,SpringIDs,Phi1,Phi2,PathPhi1,PathPhi2,pathPE_G,PE_G) 

% Design variables: [alphaA, alphaB, alphaC, alphaD, kthetaA, kthetaB, kthetaC, kthetaD]

% Linear constraints
A = []; b = []; Aeq = []; beq = []; 

% Lower and upper bounds for design variables
lb = [repmat(0,1,length(SpringIDs)) repmat(0,1,length(SpringIDs))]; 
ub = [repmat(Inf,1,length(SpringIDs)) repmat(Inf,1,length(SpringIDs))];

% Initial values for design variables 
r0 = [repmat(1,1,length(SpringIDs)) repmat(10,1,length(SpringIDs))];
r0 = [0.5835    1.3601    2.8534    1.3454   42.8253   44.6974   29.7481  208.2773];


% Objective function: sum(abs(diff(PE)))
[ObjFun] = @(r)GetObjFun(r,pathExtSpr,SpringIDs,pathPE_G);

% Constraint
[nonlcon] = @(r)GetCon(r,SpringKinematics,SpringIDs,Phi1,Phi2,PathPhi1,PathPhi2,PE_G);

% Call Optimization
[rest,fval,exitflag,output,lambda,history,~] = runfmincon(ObjFun,r0,A,b,Aeq,beq,lb,ub,nonlcon);

RestPos = rest(1:length(SpringIDs));
Stiffness = rest(length(SpringIDs)+1:end);

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
    'MaxFunctionEvaluations',3000,'MaxIterations',2000,'StepTolerance',1e-6,'ConstraintTolerance',0.05,'PlotFcn',{'optimplotx','optimplotfirstorderopt','optimplotfvalconstr'}); 
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

function [ObjFun] = GetObjFun(r,pathExtSpr,SpringIDs,pathPE_G)


PE_S_Path = (1/2).*r(length(SpringIDs)+1:end)'.*(pathExtSpr-r(1:length(SpringIDs))').^2;

PE_T_Path = sum(PE_S_Path) + pathPE_G;

PEavg = mean(PE_T_Path);

RMS = sqrt(sum((PE_T_Path-PEavg).^2/length(PE_T_Path)));

ObjFun = RMS;

end 

function [c,ceq] = GetCon(r,SpringKinematics,SpringIDs,Phi1,Phi2,PathPhi1,PathPhi2,PE_G)

    ceq = [];

    PEsum = 0;

    for m = 1:length(SpringIDs)
    
        k = r(length(SpringIDs)+m); L0a0 = r(m);
    
%         PE_S.(SpringIDs{m}) = (1/2).*k.*(SpringKinematics.(SpringIDs{m})-L0a0).^2;

        PEsum = PEsum + (1/2).*k.*(SpringKinematics.(SpringIDs{m})-L0a0).^2;
    
    end

    PE =  PEsum + PE_G;

    for m = 1:length(PathPhi1)
        hi = abs(Phi1-PathPhi1(m));
        is(m) = find(hi==min(hi));
        bi = abs(Phi2-PathPhi2(m));
        js(m) = find(bi==min(bi));

    end
    
    % calculate y component tangent to path at each point
    dPy = gradient(PathPhi2);
    
    % x component of tangent always = 1
    dPx = gradient(PathPhi1);
    
    % tangent vector
    dP = [dPx' dPy'];
    
    % compute normal to path at each point (rot by 90)
    nP = [dP(:,2) -dP(:,1)];
    
    % assign whole space to follow normal of path point closest to it
    desX = -ones(size(PE))/1000; desY = -ones(size(PE))/1000;
    for nn = 1:size(PE,2)
        for mm = 1:size(PE,1)
            if mm < js(nn)
                desX(mm,nn) = nP(nn,1);
                desY(mm,nn) = nP(nn,2);
            elseif mm > js(nn) % (rot by 180)
                desX(mm,nn) = -nP(nn,1);
                desY(mm,nn) = -nP(nn,2);
            end
        end
    end
    
    [Fx,Fy] = gradient(PE);
    
    for ii = 1:size(Fx,2)
        for jj = 1:size(Fx,1)
            gradvect = [Fx(jj,ii) Fy(jj,ii)]./(sqrt(Fx(jj,ii).^2+Fy(jj,ii).^2));
            unitvect = [desX(jj,ii) desY(jj,ii)]./(sqrt(desX(jj,ii).^2+desY(jj,ii).^2));
            c(jj,ii) = acos(dot(gradvect,unitvect)) - 4*pi/8;
        end
    end
     
%     c=real(c);
    

end