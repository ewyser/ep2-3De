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
delete('*.dat','*.out','*.exe','*.lib','*.exp');
clear,close all                                                ;%
typeD = 'single'
numel = [80];

mpiprocess = [1 2 4 8 16 32 64 128]
mpiprocess = [128]
L(1,:) = [0.25 0.5 1 2 4 8 16 19]; %L(1,:) = [0.25 0.5 1 2 4 8 16]
% L(2,:) = mpiprocess(2).*L(1,:)
% L(3,:) = mpiprocess(3).*L(1,:)
% L(4,:) = mpiprocess(4).*L(1,:)
% L(5,:) = mpiprocess(5).*L(1,:)
% L(6,:) = mpiprocess(6).*L(1,:)
% L(7,:) = mpiprocess(7).*L(1,:)
% L(8,:) = mpiprocess(8).*L(1,:)

for k = 1;%1:length(mpiprocess)
    
    nGPU = mpiprocess(k);
    l = L(k,:);
for sim=1:length(l)
    disp('------------------------')                                      ;%
    disp(['Run ',num2str(sim),': nel = ',num2str(numel),''])         ;%
    disp('------------------------')                                      ;%
    %% NON-DIMENSIONAL CONSTANT
    nu      = 0.3                                                         ;% poisson ratio
    ni      = 2                                                           ;% number of mp in h(1) direction
    %---------------------------------------------------------------------%
    %% INDEPENDANT PHYSICAL CONSTANT
    g       = 9.81                                                        ;% gravitationnal acceleration [m/s^2]
    E       = 1.0e6                                                       ;% young's modulus             [Pa]
    Gc      = E/(2*(1+nu))                                                ;% shear modulus               [Pa]
    Kc      = E/(3*(1-2*nu))                                              ;% bulk modulus                [Pa]
    rho0    = 2700                                                        ;% density                     [kg/m^3]
    n0      = 0.0                                                         ;% initial porosity            [-]
    yd      = sqrt((Kc+4/3*Gc)/rho0)                                      ;% elastic wave speed          [m/s]
    coh0    = 20.0e3                                                      ;% cohesion                    [Pa]
    phi0    = 20.0*pi/180                                                 ;% friction angle              [Rad]
    psi0    =  0.0*pi/180                                                 ;%
    H       = -60.0e3                                                      ;% softening modulus           [Pa]
    cohr    =  4.0e3                                                      ;% residual cohesion           [Pa]
    phir    = 7.0*pi/180                                                  ;% residual friction angle     [Rad]
    t       = 1.1                                                        ;% simulation time             [s]
    te      = 1.1                                                         ;% elastic loading             [s]
    tg      = te/1.5                                                      ;% gravity increase            [s]
    %---------------------------------------------------------------------%
    %% MESH & MP INITIALIZATION
    lx      = 64.1584                                                     ;% mesh size along x direction
    ly      = lx                                                          ;% mesh size along y direction
    ly      = l(sim)*lx                                                        ;% mesh size along y direction
    lz      = 12.80                                                       ;% mesh size along z direction
    [h]     = setupMPI(nGPU,numel,lx,ly,lz,typeD,ni,coh0,cohr,phi0,phir,n0,rho0)
    % ISOTROPIC ELASTIC MATRIX
    Del     = [ Kc+4/3*Gc,Kc-2/3*Gc,Kc-2/3*Gc,0.0,0.0,0.0;...
                Kc-2/3*Gc,Kc+4/3*Gc,Kc-2/3*Gc,0.0,0.0,0.0;...
                Kc-2/3*Gc,Kc-2/3*Gc,Kc+4/3*Gc,0.0,0.0,0.0;...
                0.0      ,0.0      ,0.0      ,Gc ,0.0,0.0;...                             
                0.0      ,0.0      ,0.0      ,0.0,Gc ,0.0;...
                0.0      ,0.0      ,0.0      ,0.0,0.0,Gc ;]               ;%     
    Hp      = H*h(1)                                                      ;%
    %---------------------------------------------------------------------%
    % physics
    fid = fopen( 'Del.dat','wb')                                          ;% export
    fwrite(fid, Del(:),typeD)                                             ;% export
    fclose(fid)                                                           ;% export
    iDel = inv(Del)                                                       ;%
    fid = fopen( 'iDel.dat','wb')                                         ;% export
    fwrite(fid, iDel(:),typeD)                                            ;% export
    fclose(fid)                                                           ;% export
    p   = [g;rho0;psi0;nu;E;Kc;Gc;cohr;Hp;t;te;tg]                        ;%
    fid = fopen( 'phys.dat','wb')                                         ;% export
    fwrite(fid,p ,typeD)                                                  ;% export
    fclose(fid)                                                           ;% export
    %% MPM DM ALGORITHM EXPLICIT SOLVER FOR INFINITESIMAL STRAIN
BCK1  = 16;
BCK2  = 256;
% compile gpu code
arg  = [' ' '',num2str(sim),' float 1 0.1 ',num2str(BCK1),' ',num2str(BCK2),' ',num2str(nGPU),' 9'];   
system(['chmod +x ' pwd '/shellscript_prof.sh']);
pathToScript = fullfile([pwd '/shellscript_prof.sh']);
shell_run = [pathToScript arg];
 system(shell_run);
    %---------------------------------------------------------------------%    
    %% DISPLAY
np = 0;
for i=1:nGPU
    fid  = fopen(['GPUinfo_' num2str(i-1) '.dat']); GPUtime = fread(fid,typeD); fclose(fid);
    fid = fopen( ['param_rank_',num2str(i-1),'.dat']); p = fread(fid,typeD); fclose(fid);
    time(:,i) = [GPUtime];
    np = np+p(1);
end
  



Perf(:,sim) = [np;mean(time(1:2,:),2);sum(time(3:4,:),2)];
% np                          
% runtime                     |average
% it/s                              |average
% GBs-1                         |sum
% memory footprint  |average
end



if(nGPU==1)
    filename = 'Performance_oneGPU'
    save(filename,'Perf')
else
    filename = ['Performance_',num2str(nGPU),'GPUs']
    save(filename,'Perf')
end
end


