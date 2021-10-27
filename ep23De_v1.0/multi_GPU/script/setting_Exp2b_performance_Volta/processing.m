%% MADMAX
% Copyright (C) 2020: Emmanuel Wyser      , emmanuel.wyser[at]unil.ch
%                     Yury Alkhimenkov    , yury.alkhimenkov[at]unil.ch
%                     Michel Jaboyedoff   , michel.jaboyedoff[at]unil.ch
%                     Yury Y. Podladchikov, yury.podladchikov[at]unil.ch
% -------------------------------------------------------------------------%
% version    : v1.0
% date       : february, 2021
% description: explicit mpm (GIMPM) solver based on an updated Lagrangian
% frame for elasto-plastic problem discretized on a 4-noded quadrilateral 
% mesh
% -------------------------------------------------------------------------%
clear,close all                                                ;%
typeD = 'single'
numel = [320];
nGPU  = 4;

 for i=1:nGPU   

    
         fid  = fopen(['xp_rank_',num2str(i-1),'.dat']);xp = fread(fid,typeD);fclose(fid);
         fid  = fopen(['up_rank_',num2str(i-1),'.dat']);up = fread(fid,typeD);fclose(fid);
         fid  = fopen(['epII_rank_',num2str(i-1),'.dat']);epII = fread(fid,typeD);fclose(fid);
        fid  = fopen(['sig_rank_',num2str(i-1),'.dat']);sig = fread(fid,typeD);fclose(fid);
        s    = reshape(sig,length(sig)/6,6)                  ;% stresses
        xp = reshape(xp,length(xp)/3,3);
        up = reshape(up,length(up)/3,3);
        p   = (s(:,1)+s(:,2)+s(:,3))/3;
        du   = sqrt(up(:,2).^2);
        
        figure(642346)
        hold on
        scatter3(xp(:,1),xp(:,2),xp(:,3),5,epII,'filled')
        axis equal
        hold off
        colormap(jet(1000))
        colorbar
        axis tight
         view(47,13);
         view(0,0)
        drawnow
 end
  




% writematrix([x p0],'myData.txt','Delimiter',';'); 
% type myData.txt;
