function [] = GetQuiver_SS(Phi_A_1,Phi_D_2,PE_T)
    
    desX = ones(size(PE_T))*-1;
    desY = zeros(size(PE_T));
    
    [Fx,Fy] = gradient(PE_T);
     
    figure('Position',[-798,-163,560,420])
    quiver(Phi_A_1(1:5:end)*180/pi,Phi_D_2(1:5:end)*180/pi,desX(1:5:end,1:5:end),desY(1:5:end,1:5:end),'Color',[0.5 0.5 0.5])
    hold on
    quiver(Phi_A_1(1:5:end)*180/pi,Phi_D_2(1:5:end)*180/pi,Fx(1:5:end,1:5:end),Fy(1:5:end,1:5:end),'r')
%     xlim([min(Phi_A_1*180/pi) max(Phi_A_1*180/pi)]);
%     ylim([min(Phi_D_2*180/pi) max(Phi_D_2*180/pi)]);

end