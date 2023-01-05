%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%------------------------ Knee Model ----------------------------%%%
%%%---------------------- Maria Redoutey --------------------------%%%
%%%--------------------- 22 February 2021 -------------------------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ThetaA,ThetaB,ThetaC,ThetaD,ThetaE,PEgrav,L1,L2] = GetGeo_Excavator(PSI,PHI)


    % Bar Lengths [m]
    L = 1 ; 
    
    % Boom Length [m]
    Lb = 5.7; %6
    
    %  Arm Length [m]
    La = 2.9; % 4
       
    % Pin Support Coordinates
    A = [Lb-L 0]; D = [Lb 0];
    
    for j = 1:length(PSI)
        for i = 1:length(PHI)

            C = D + [L*cos(pi-PHI(i)) L*sin(pi-PHI(i))];

            Bucket = D + [La*cos(pi-PHI(i)) La*sin(pi-PHI(i))];
            
            % Cab Connection Coordinates
            Cab = [0 0];

            B = C - [L 0];

            coords = [A;B;C;D];

            R = [cos(PSI(j)) -sin(PSI(j));
                    sin(PSI(j)) cos(PSI(j))];

            pts = [coords(:,1) coords(:,2)];
            pts2 = R*transpose(pts);
            coords = transpose([pts2(1,:); pts2(2,:)]);

            Bucket2 = R*transpose(Bucket);
            Bucket = transpose([Bucket2(1,:); Bucket2(2,:)]);
            
%             Cab2 = R*transpose(Cab);
%             Cab = transpose([Cab2(1,:); Cab2(2,:)]);

            L1(j,i) = sqrt(coords(2,1)^2+coords(2,2)^2);
            Back = Cab + [3.65 0];
            L2(j,i) = sqrt((coords(4,1)-Back(1))^2+(coords(4,2)-Back(2))^2);
%             L2(j,i) = sqrt((coords(2,1)-Back(1))^2+(coords(2,2)-Back(2))^2);

            
            coords(:,2) = coords(:,2)+2;
            Bucket(2) = Bucket(2)+2;
            Cab(2) = Cab(2)+2;
            Back(2) = Back(2)+2;


%             if i==j
%                 figure()
%                 plot([Cab(1); coords(4,1)],[Cab(2); coords(4,2)],'k',[Bucket(1); coords(4,1)],[Bucket(2); coords(4,2)],'k','LineWidth',2)
%                 hold on
%                 %plot([coords(:,1); coords(1,1)],[coords(:,2); coords(1,2)],'bo-',[coords(1,1);coords(4,1)],[coords(1,2);coords(4,2)],'b',[coords(3,1);coords(4,1)],[coords(3,2);coords(4,2)],'b','LineWidth',1)
%                 %plot([Cab(1); coords(2,1)],[Cab(2), coords(2,2)],'r',[Back(1); coords(2,1)],[Back(2), coords(2,2)],'r')
% %                 scatter([Back(1) Cab(1) coords(4,1) Bucket(1)],[Back(2) Cab(2) coords(4,2) Bucket(2)])
%                 axis equal
%                 axis off
%             end
                %drawnow
%             elseif i==1 && j==1
% %                 figure()
%                 plot([Cab(1); coords(4,1)],[Cab(2); coords(4,2)],'k',[Bucket(1); coords(4,1)],[Bucket(2); coords(4,2)],'k','LineWidth',2)
%                 hold on
%                 plot([coords(:,1); coords(1,1)],[coords(:,2); coords(1,2)],'bo-',[coords(1,1);coords(4,1)],[coords(1,2);coords(4,2)],'b',[coords(3,1);coords(4,1)],[coords(3,2);coords(4,2)],'b','LineWidth',1)
%                 plot([Cab(1); coords(2,1)],[Cab(2), coords(2,2)],[Back(1); coords(4,1)],[Back(2), coords(4,2)],'r')                
%                 axis equal
%                 axis off
%             elseif i==length(PHI) && j==length(PSI)
% %                 figure()
%                 plot([Cab(1); coords(4,1)],[Cab(2); coords(4,2)],'k',[Bucket(1); coords(4,1)],[Bucket(2); coords(4,2)],'k','LineWidth',2)
%                 hold on
%                 plot([coords(:,1); coords(1,1)],[coords(:,2); coords(1,2)],'bo-',[coords(1,1);coords(4,1)],[coords(1,2);coords(4,2)],'b',[coords(3,1);coords(4,1)],[coords(3,2);coords(4,2)],'b','LineWidth',1)
%                 plot([Cab(1); coords(2,1)],[Cab(2), coords(2,2)],[Back(1); coords(4,1)],[Back(2), coords(4,2)],'r')
%                 axis equal
%                 axis off
%             elseif i==length(PHI) && j==1
% %                 figure()
%                 plot([Cab(1); coords(4,1)],[Cab(2); coords(4,2)],'k',[Bucket(1); coords(4,1)],[Bucket(2); coords(4,2)],'k','LineWidth',2)
%                 hold on
%                 plot([coords(:,1); coords(1,1)],[coords(:,2); coords(1,2)],'bo-',[coords(1,1);coords(4,1)],[coords(1,2);coords(4,2)],'b',[coords(3,1);coords(4,1)],[coords(3,2);coords(4,2)],'b','LineWidth',1)
%                 plot([Cab(1); coords(2,1)],[Cab(2), coords(2,2)],[Back(1); coords(4,1)],[Back(2), coords(4,2)],'r')
%                 axis equal
%                 axis off
%             end

            boom_mass = 1710; %[kg] 1000
            arm_mass = 1300; %[kg] boom_mass*La/Lb
            bucket_mass = 787; %[kg] 1000

            PEgrav(j,i) = 9.81*(boom_mass*((Cab(2)+coords(4,2))/2) + arm_mass*((Bucket(2)+coords(4,2))/2) + bucket_mass*Bucket(2));
            
            ThetaD(j,i) = PHI(i);
            ThetaB(j,i) = PHI(i);
            ThetaA(j,i) = pi-PHI(i);
            ThetaC(j,i) = pi-PHI(i);
            ThetaE(j,i) = PSI(j);

        end 
        
    end
end 