

function [a1,a2,a3,a4,k1,k2,k3,k4,history,exitflag,output,lambda] = FindRest_FiveFoldVertex(Theta3,Theta4,Theta7,Theta9,PE_G) 

% Design variables: [alphaA, alphaB, alphaC, alphaD, kthetaA, kthetaB, kthetaC, kthetaD]

% Linear constraints
A = []; b = []; Aeq = []; beq = []; 

% Lower and upper bounds for design variables
lb = [min(min(Theta3)),min(min(Theta4)),min(min(Theta7)),min(min(Theta9)),0,0,0,0]; ub = [max(max(Theta3)),max(max(Theta4)),max(max(Theta7)),max(max(Theta9)),Inf,Inf,Inf,Inf];

% Initial values for design variables 
r0 = [(ub(1)-lb(1))/2+lb(1),(ub(2)-lb(2))/2+lb(2),(ub(3)-lb(3))/2+lb(3),(ub(4)-lb(4))/2+lb(4),10,10,10,10]; 

% Objective function: sum(abs(diff(PE)))
% [ObjFun] = @(r)GetObjFun_SAD(r,Theta3,Theta4,Theta7,Theta9,PE_G);
[ObjFun] = @(r)GetObjFun_FiveFoldVertex(r,Theta3,Theta4,Theta7,Theta9,PE_G);

[nonlcon] = @(r)GetCon(r,Theta3,Theta4,Theta7,Theta9,PE_G);

% Call Optimization
[rest,~,exitflag,output,lambda,history,~] = runfmincon(ObjFun,r0,A,b,Aeq,beq,lb,ub,[]);
a1 = rest(1); a2 = rest(2); a3 = rest(3); a4 = rest(4); k1 = rest(5); k2 = rest(6); k3 = rest(7); k4 = rest(8);

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
options = optimoptions(@fmincon,'OutputFcn',@outfun,'Display','iter','Algorithm','interior-point','FiniteDifferenceStepSize',1e-6,'MaxFunctionEvaluations',12000,'MaxIterations',2000,'StepTolerance',1e-6); 
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

% function [ObjFun] = GetObjFun(r,Theta3,Theta4,Theta7,Theta9,PE_G)
% 
%     PE = (1/2)*r(5)*(Theta3-r(1)).^2 + (1/2)*r(6)*(Theta4-r(2)).^2 + (1/2)*r(7)*(Theta7-r(3)).^2 + (1/2)*r(8)*(Theta9-r(4)).^2 + PE_G;
%     PE=PE(2:end,:);
% 
% %     for i = 1:size(PE,1)
% %         for j = 1:size(PE,2)
% %             if isnan(PE(i,j))
% %                 PE(i,j) = PE(i,j-1);
% %             end
% %         end
% %     end
%     
%     for i = 2:size(PE,1)
%         dPEd1 = diff(PE(~isnan(PE(i,:))));
%         SAdPEd1(i,:) = sum(abs(dPEd1));
%     end
%     dSAdPEd1d2 = diff(SAdPEd1);
%     SAdSAdPEd1d2 = sum(abs(dSAdPEd1d2));
%     
%     for j = 1:size(PE,2)-1
%         dPEd2 = diff(PE(~isnan(PE(:,j))));
%         SAdPEd2(:,j) = sum(abs(dPEd2));
%     end
%     dSAdPEd2d1 = diff(SAdPEd2);
%     SAdSAdPEd2d1 = sum(abs(dSAdPEd2d1));
%     
%     ObjFun = sum([SAdSAdPEd1d2,SAdSAdPEd2d1]);
%     
% %     dPEd1 = diff(PE);
% %     dPEd2 = diff(PE');
% 
% %     SAdPEd1 = sum(abs(dPEd1));
% %     SAdPEd2 = sum(abs(dPEd2));
% 
% %     dSAdPEd1d2 = diff(SAdPEd1);
% %     dSAdPEd2d1 = diff(SAdPEd2);
% 
% %     SAdSAdPEd1d2 = sum(abs(dSAdPEd1d2));
% %     SAdSAdPEd2d1 = sum(abs(dSAdPEd2d1));
% 
%     % I think this isn't accurate because of diff and how I'm currently
%     % dealing with nans. go thru row by row (or col by col) and take diff,
%     % sum abs, etc. end up with a big for loop but essentially the same
%     % thing, just separated out bc diff row/col sizes as we go thru.
%     
%     for i = 1:size(PE,1)
%         for j = 1:size(PE,2)
%             if isnan(PE(i,j))
%                 PE(i,j) = PE(i,j-1);
% %                 PE(i,j) = 0;
%             end
%         end
%     end
% 
%     PEavg = mean(mean(PE));
%     dev = PE-PEavg;
%     dev2 = dev.^2;
%     RMS = sqrt(mean(mean(dev2)));
%     ObjFun = RMS;
% end

function [c,ceq] = GetCon(r,Theta3,Theta4,Theta7,Theta9,PE_G)

    PE = (1/2)*r(5)*(Theta3-r(1)).^2 + (1/2)*r(6)*(Theta4-r(2)).^2 + (1/2)*r(7)*(Theta7-r(3)).^2 + (1/2)*r(8)*(Theta9-r(4)).^2 + PE_G;
    
    gPE = gradient(PE); gPE_G = gradient(PE_G);
%     c = sum(abs(gPE(~isnan(gPE)))) - sum(abs(gPE_G(~isnan(gPE_G))));

%     ceq = sum(abs(gPE(~isnan(gPE)))) ;

ceq=[]; 
c = [];

    [fx,fy] = gradient(PE);
  
    fx(isnan(fx))=0; fy(isnan(fy))=0;
    
%     ceq = mean(sum(abs(fx))) + mean(sum(abs(fy))) - 0.8;

end

function [ObjFun] = GetObjFun_FiveFoldVertex(r,Theta3,Theta4,Theta7,Theta9,PE_G)

    PE = (1/2)*r(5)*(Theta3-r(1)).^2 + (1/2)*r(6)*(Theta4-r(2)).^2 + (1/2)*r(7)*(Theta7-r(3)).^2 + (1/2)*r(8)*(Theta9-r(4)).^2 + PE_G;

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
