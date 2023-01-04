function [] = GetQuiver(Phi1,PathPhi1,Phi2,PathPhi2,PE_T)

    for m = 1:length(PathPhi1)
        hi = abs(Phi1-PathPhi1(m));
        is(m) = find(hi==min(hi));
        bi = abs(Phi2-PathPhi2(m));
        js(m) = find(bi==min(bi));

        pathPE(m) = PE_T(js(m),is(m));
    end
    
    % calculate y component tangent to path at each point
    dPy = gradient(PathPhi2);
    
    % x component of tangent 
    dPx = gradient(PathPhi1);
    
    % tangent vector
    dP = [dPx' dPy'];
    
    % compute normal to path at each point (rot by 90)
    nP = [dP(:,2) -dP(:,1)];
    
    % assign whole space to follow normal of path point closest to it
    desX = -ones(size(PE_T))/1000; desY = -ones(size(PE_T))/1000;
    for nn = 1:size(PE_T,2)
        for mm = 1:size(PE_T,1)
            if mm < js(nn)
                desX(mm,nn) = nP(nn,1);
                desY(mm,nn) = nP(nn,2);
            elseif mm > js(nn) % (rot by 180)
                desX(mm,nn) = -nP(nn,1);
                desY(mm,nn) = -nP(nn,2);
            end
        end
    end
    
    [Fx,Fy] = gradient(PE_T);
    
%     for ii = 1:size(Fx,2)
%         for jj = 1:size(Fx,1)
%             gradvect = [Fx(jj,ii) Fy(jj,ii)];
%             unitvect = [desX(jj,ii) desY(jj,ii)];
%             c(jj,ii) = acos(dot(gradvect,unitvect)) - pi/2;
%         end
%     end
%     
    figure('Position',[-798,-163,560,420])
    quiver(Phi1*180/pi,Phi2*180/pi,desX,desY,'b')
    hold on
    quiver(Phi1*180/pi,Phi2*180/pi,Fx,Fy,'r')
    


end