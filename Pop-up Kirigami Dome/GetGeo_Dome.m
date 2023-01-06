function [Phi_Assemble_Vector,Phi_Dome_Vector,Length_LR,Length_TB,PE_G,Length_Span,Length_Rise] = GetGeo_Dome()

% Assume panels have uniform mass density of 1 kg/m^2. So area = mass

L = 1.5; g = 20*pi/180;

tc = 2*acos(tan(pi/3)*tan(g));

Phi_Assemble_Vector = 0:0.05:tc;
Phi_Dome_Vector = 0:0.01:0.3424; % obtained from bar & hinge simulations

PE_G_Assemble = zeros(length(Phi_Dome_Vector),length(Phi_Assemble_Vector));
PE_G_Dome = zeros(size(PE_G_Assemble));
Length_LR = zeros(size(PE_G_Assemble));
Length_TB = zeros(size(PE_G_Assemble));

for i = 1:length(Phi_Assemble_Vector)
    
    Phi_Assemble = Phi_Assemble_Vector(i);

    s=1; w=sqrt(3)*s; h=2*s; hor = w; vert = h*3/4;

    base_coords = [0 s 0;
                   hor/2 vert/3 0;
                   hor/2 -vert/3 0;
                   0 -s 0;
                   -hor/2 -vert/3 0;
                   -hor/2 vert/3 0];

    side_coords =  [base_coords(2,1)+L*cos(Phi_Assemble/2) base_coords(2,2)+L*tan(g) L*sin(Phi_Assemble/2);
                    base_coords(3,1)+L*cos(Phi_Assemble/2) base_coords(3,2)-L*tan(g) L*sin(Phi_Assemble/2);
                    base_coords(5,1)-L*cos(Phi_Assemble/2) base_coords(5,2)-L*tan(g) L*sin(Phi_Assemble/2);
                    base_coords(6,1)-L*cos(Phi_Assemble/2) base_coords(6,2)+L*tan(g) L*sin(Phi_Assemble/2)];

    outer_coords = [side_coords(1:2,:); 
                    side_coords(1:2,:)*[cosd(60) sind(60) 0; -sind(60) cosd(60) 0; 0 0 1];
                    side_coords(1:2,:)*[cosd(120) sind(120) 0; -sind(120) cosd(120) 0; 0 0 1];
                    side_coords(3:4,:); 
                    side_coords(3:4,:)*[cosd(60) sind(60) 0; -sind(60) cosd(60) 0; 0 0 1];
                    side_coords(3:4,:)*[cosd(120) sind(120) 0; -sind(120) cosd(120) 0; 0 0 1]];

    top_coords = base_coords; top_coords(:,3) = 2*L*sin(Phi_Assemble/2);
    
    coords = [base_coords; side_coords; outer_coords; top_coords];
    
    panels = [1 2 14 13; 
                2 3 12 11; 
                3 4 22 21; 
                4 5 20 19; 
                5 6 18 17;
                6 1 16 15;
                23 24 14 13; 
                24 25 12 11; 
                25 26 22 21; 
                26 27 20 19; 
                27 28 18 17;
                28 23 16 15];

    Length_LR(:,i) = w + 2*L*cos(Phi_Assemble/2);
    Length_TB(:,i) = 2*L*sin(Phi_Assemble/2);

    % Panel Areas
    HexArea = 3*sqrt(3)/2*s^2;
    TrapArea = 0.5*L*(s+s+2*L*tan(g));
    
    % Height of panel center of mass
    HexZ = 2*L*sin(Phi_Assemble/2);
    TrapZBot = (L-(L*(2*s+s+2*L*tan(g)))/(3*(s+s+2*L*tan(g))))*sin(Phi_Assemble/2);
    TrapZTop = L*sin(Phi_Assemble/2) + (L*(2*s+s+2*L*tan(g)))/(3*(s+s+2*L*tan(g)))*sin(Phi_Assemble/2);

    % PE_G = m*g*h
    PE_G_Assemble(:,i) = 7*9.81*(HexArea*HexZ + 6*TrapArea*TrapZBot + 6*TrapArea*TrapZTop);
        
    for j = 1:length(Phi_Dome_Vector)
        
        HexZCenterTop = 2*L*sin(Phi_Assemble/2) + 2*(w/2 + L*sin(Phi_Assemble/2))*sin(Phi_Dome_Vector(j));
        HexZCenterBot = 2*(w/2 + L*sin(Phi_Assemble/2))*sin(Phi_Dome_Vector(j));

        TrapZBotCenter = (L-(L*(2*s+s+2*L*tan(g)))/(3*(s+s+2*L*tan(g))))*sin(Phi_Assemble/2) + 2*(w/2 + L*sin(Phi_Assemble/2))*sin(Phi_Dome_Vector(j));
        TrapZTopCenter = L*sin(Phi_Assemble/2) + (L*(2*s+s+2*L*tan(g)))/(3*(s+s+2*L*tan(g)))*sin(Phi_Assemble/2) + 2*(w/2 + L*sin(Phi_Assemble/2))*sin(Phi_Dome_Vector(j));

        HexZSideTop = 2*L*sin(Phi_Assemble/2) + (L*sin(Phi_Assemble/2)+w/2)*sin(Phi_Dome_Vector(j));
        HexZSideBot = (L*sin(Phi_Assemble/2)+w/2)*sin(Phi_Dome_Vector(j));

        TrapZBotInside = (L-(L*(2*s+s+2*L*tan(g)))/(3*(s+s+2*L*tan(g))))*sin(Phi_Assemble/2) + (L*sin(Phi_Assemble/2)+w+L*sin(Phi_Assemble/2)-(L*(2*s+s+2*L*tan(g)))/(3*(s+s+2*L*tan(g))))*sin(Phi_Dome_Vector(j));
        TrapZTopInside = L*sin(Phi_Assemble/2) + (L*(2*s+s+2*L*tan(g)))/(3*(s+s+2*L*tan(g)))*sin(Phi_Assemble/2) + (L*sin(Phi_Assemble/2)+w+L*sin(Phi_Assemble/2)-(L*(2*s+s+2*L*tan(g)))/(3*(s+s+2*L*tan(g))))*sin(Phi_Dome_Vector(j));

        TrapZBotOutside = (L-(L*(2*s+s+2*L*tan(g)))/(3*(s+s+2*L*tan(g))))*sin(Phi_Assemble/2)+  (L*(2*s+s+2*L*tan(g)))/(3*(s+s+2*L*tan(g)))*sin(Phi_Dome_Vector(j));
        TrapZTopOutside = L*sin(Phi_Assemble/2) + (L*(2*s+s+2*L*tan(g)))/(3*(s+s+2*L*tan(g)))*sin(Phi_Assemble/2)+  (L*(2*s+s+2*L*tan(g)))/(3*(s+s+2*L*tan(g)))*sin(Phi_Dome_Vector(j));

        PE_G_Center = 9.81*(HexArea*HexZCenterTop + HexArea*HexZCenterBot + 6*TrapArea*TrapZBotCenter + 6*TrapArea*TrapZTopCenter);

        PE_G_Side = 9.81*(6*HexArea*HexZSideTop + 6*HexArea*HexZSideBot + 3*6*TrapArea*TrapZBotInside + 3*6*TrapArea*TrapZTopInside + 3*6*TrapArea*TrapZBotOutside + 3*6*TrapArea*TrapZTopOutside);

        PE_G_Dome(j,i) = PE_G_Center + PE_G_Side;
                
        Length_Span(j,i) = 2*((0.5*Length_LR(j,i))+Length_LR(j,i)*cos(Phi_Dome_Vector(j)));
        Length_Rise(j,i) = HexZCenterBot;
        
%         figure()
%         patch('Faces',panels,'Vertices',coords)
%         axis equal
    
    end
end

PE_G = PE_G_Dome;

end