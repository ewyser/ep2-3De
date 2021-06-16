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
clear,close,clf                                                           ;%
delete('*.dat','*.out','*.exe','*.avi','*.lib','*.exp');
run('preprocessing');%
run('compileExp1a');%

numel = [40] 
for sim=1:length(numel)
    disp('------------------------')                                      ;%
    disp(['Run ',num2str(sim),': nel = ',num2str(numel(sim)),''])         ;%
    disp('------------------------')                                      ;%
    %% NON-DIMENSIONAL CONSTANT
    nu      = 0.0                                                         ;% poisson ratio
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
    coh0    = 0.0e3                                                       ;% cohesion                    [Pa]
    phi0    = 19.8*pi/180                                                 ;% friction angle              [Rad]
    psi0    =  0.0*pi/180                                                 ;%
    H       = -0.0e3                                                      ;% softening modulus           [Pa]
    cohr    =  0.0e3                                                      ;% residual cohesion           [Pa]
    phir    = 7.0*pi/180                                                  ;% residual friction angle     [Rad]
    t       = 1.0                                                         ;% simulation time             [s]
    te      = 0.001                                                       ;% elastic loading             [s]
    tg      = te/1.5                                                      ;% gravity increase            [s]
    %---------------------------------------------------------------------%
    %% MESH & MP INITIALIZATION
    length  = 200e-3
    height  = 0.5*length
    width   = length/4
    lx      = length                                                      ;% mesh size along x direction
    ly      = width                                                       ;% mesh size along y direction
    lz      = 1.1*height                                                  ;% mesh size along z direction
    [meD,bc]= meSetup(numel(sim),length,height,width,typeD)                          ;% - see function                                                      ;%
    [mpD]   = mpSetup(meD,ni,width,height,length,coh0,cohr,phi0,phir,n0,rho0,typeD)    ;% - see function
    
%     figure
%     plot3(mpD.x(:,1),mpD.x(:,2),mpD.x(:,3),'o',meD.x([bc.x1;bc.x2;bc.y1;bc.y2]),meD.y([bc.x1;bc.x2;bc.y1;bc.y2]),meD.z([bc.x1;bc.x2;bc.y1;bc.y2]),'s'),axis equal
%     figure
%     plot(mpD.x(:,1),mpD.x(:,3),'o',meD.x([bc.x1;bc.x2;bc.y1;bc.y2]),meD.z([bc.x1;bc.x2;bc.y1;bc.y2]),'s'),axis equal

    
    
    % ISOTROPIC ELASTIC MATRIX
    Del     = [ Kc+4/3*Gc,Kc-2/3*Gc,Kc-2/3*Gc,0.0,0.0,0.0;...
                Kc-2/3*Gc,Kc+4/3*Gc,Kc-2/3*Gc,0.0,0.0,0.0;...
                Kc-2/3*Gc,Kc-2/3*Gc,Kc+4/3*Gc,0.0,0.0,0.0;...
                0.0      ,0.0      ,0.0      ,Gc ,0.0,0.0;...                             
                0.0      ,0.0      ,0.0      ,0.0,Gc ,0.0;...
                0.0      ,0.0      ,0.0      ,0.0,0.0,Gc ;]               ;%     
    Hp      = H*meD.h(1)                                                  ;%
    %---------------------------------------------------------------------%
    %% MPM DM ALGORITHM EXPLICIT SOLVER FOR INFINITESIMAL STRAIN
    run('export');
    % CUDA CODE
    if(ismac || isunix)
        system('./gpu.out');
    elseif(ispc)
        system('gpu.exe');
    else
        disp('Outer space OS');
    end
    %---------------------------------------------------------------------%    
    %% DISPLAY
    fid  = fopen('sig.dat')    ; sig     = fread(fid,typeD); fclose(fid);
    fid  = fopen('epII.dat')   ; epII    = fread(fid,typeD); fclose(fid);   
    fid  = fopen('xp.dat')     ; xp      = fread(fid,typeD); fclose(fid);
    fid  = fopen('xp0.dat')    ; xp0     = fread(fid,typeD); fclose(fid);
    fid  = fopen('lp.dat')     ; lp      = fread(fid,typeD); fclose(fid);
    fid  = fopen('up.dat')     ; up      = fread(fid,typeD); fclose(fid);
    fid  = fopen('GPUinfo.dat'); GPUtime = fread(fid,typeD); fclose(fid);
    s    = reshape(sig,mpD.n,6)                  ;% stresses
    epII = epII                                  ;% accumulated plastic strain
    x    = reshape(xp ,mpD.n,3)                  ;% coordinates
    x0   = reshape(xp0,mpD.n,3)                  ;% coordinates
    up   = reshape(up ,mpD.n,3)                  ;% displacement
    l    = reshape(lp ,mpD.n,3)                  ;% domain lengths
    p0   = -((s(:,1)+s(:,2)+s(:,3))./3)./1e3     ;% pressure
    du   = sqrt(up(:,1).^2+up(:,2).^2+up(:,3).^2);% L2 norm of displacement
    
    name=['Exp1a_D_' num2str(mpD.n) 'np.mat'];
    save(name,'lz','mpD','meD','s','epII','x','du','GPUtime','-v7.3');
    run('postprocessing')
end

delete('*.dat','*.out','*.exe','*.avi','*.lib','*.exp');