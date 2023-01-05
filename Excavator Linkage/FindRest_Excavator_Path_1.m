%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%--------------- Miura-ori Cell: Zero Stiffness -----------------%%%
%%%---------------------- Maria Redoutey --------------------------%%%
%%%---------------------- 10 March 2021 ---------------------------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rest,history,exitflag,output,lambda] = FindRest_Excavator_Path_1(PHI,PSI,pathThetaA,pathThetaB,pathThetaC,pathThetaD,pathThetaE,pathL1,pathL2,pathPE_G,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,L1,L2,PE_G) 

% Design variables: [alphaA, alphaB, alphaC, alphaD, kthetaA, kthetaB, kthetaC, kthetaD, L0, k]

% Linear constraints
A = []; b = []; Aeq = []; beq = []; 

% Lower and upper bounds for design variables
lb = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]; ub = [2*pi,2*pi,2*pi,2*pi,2*pi,Inf,Inf,Inf,Inf,Inf,15,15,Inf,Inf]; 

% Initial values for design variables 
r0 = [pi/2,pi/2,pi/2,pi/2,pi/2,10,10,10,10,10,1,1,10,10]; 

% Excavator . mat
r0 = [1.78444017481114e-07,3.00344857071251,6.28,6.23459558048198,1.16439515974214,926.254983958261,12626.0235222638,1.25021925385406e-07,14627.3692689629,120383.764744539,0.798689672169768,1,13358.8955692643,10];
% r0 = [1.3338e-06 2.9770 8.2814e-07 6.2726 1.1465 895.6119 1.2400e+04 0.1007 1.4376e+04 1.1947e+05 0.7877 10.0000 1.2937e+04 9.7225e-05];
% r0 = [6.28 3.4445 6.28 6.28 1.1404 908.8929 1.2396e+04 46.4516 1.4255e+04 1.1953e+05 0.7276 9.9293 1.2796e+04 215.7496];
% r0 = [1.98692815333144,1.25771749998059,6.28318193195158,3.90633015164068,0.584016707296362,4345.93126438677,4119.78398723900,1.29754965931997e-05,19985.7762507079,132085.182701253,2.92511824688227,9.99999999705444,4640.12105382187,50067.8637717928];
% r0 = [0.0465617304149103,6.28133027662302,6.27442685648626,2.64021535551548,0.647032864710834,4158.70194584023,3955.56549514774,0.327461299119964,19780.9595504521,132117.974674019,3.15983001800099,9.99999922717435,4334.90692456585,50279.8547448374];
% r0 = [0.1726 6.0044 5.9368 2.8866 0.6827 4.1586e+03 3.9555e+03 0.4979 1.9781e+04 1.3212e+05 2.6861 9.9 4.3348e+03 5.0280e+04];
% r0 = [0.4535 5.7840 5.4496 3.0778 0.7286 4.1560e+03 3.9529e+03 0.4202 1.9779e+04 1.3212e+05 2.3390 9.8958 4.3291e+03 5.0281e+04];
% r0 = [0.6613 5.5589 5.4067 3.2984 0.7742 4.1568e+03 3.9537e+03 0.3245 1.9780e+04 1.3212e+05 1.7393 9.8380 4.3290e+03 5.0281e+04];
% r0 = [0.9163 5.3848 5.2517 3.3903 0.7761 4.1562e+03 3.9531e+03 0.3467 1.9779e+04 1.3212e+05 1.6526 9.8314 4.3114e+03 5.0284e+04];
% r0 = [0.000220621533866321,4.50491221997101,6.28312902542478,3.65252598673370,0.682998485800713,8532.22038732978,8344.16120688234,4294.23540802333,24116.3573669565,142412.344591786,3.64252164952834e-07,9.99999998140965,4826.77150731224,63805.4577809973];

r0 = [4.80016846697659e-09,2.54063875579716,6.28,6.28,0.234824060292108,4.44501140893598e-07,24506.3950190332,5.94331221192344e-08,15733.1279475424,88463.7140765357,1.95878856846134,14.9999999962999,20758.6117933905,4385.48956843238];


%Excavator 2nd results . mat
% r0 = [[0.0740 3.8357 6.1311 5.5341 1.1587 914.4794 1.2614e+04 0.0172 1.4616e+04 1.2039e+05 0.7879 1.3343e+04]];

% Objective function - sum(abs(diff(PE)))
ObjFun = @(r)GetObjFun_Excavator_Path(r,pathThetaA,pathThetaB,pathThetaC,pathThetaD,pathThetaE,pathL1,pathL2,pathPE_G);

nonlcon = @(r)GetCon_Excavator_1(r,PHI,PSI,ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,L1,L2,PE_G);

% Call Optimization
[rest,~,exitflag,output,lambda,history,~] = runfmincon(ObjFun,r0,A,b,Aeq,beq,lb,ub,nonlcon);

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
    'MaxFunctionEvaluations',6000,'MaxIterations',2000,'StepTolerance',1e-6,'ConstraintTolerance',0.71,'PlotFcn',{'optimplotx','optimplotfirstorderopt','optimplotfvalconstr'}); 
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