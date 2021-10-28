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
run('preprocessing')

numEl = [80] 
lambda0 = ones(size(numEl));
for sim=1:length(numEl)
    %% COMPILE GPU CODE
    run('compileExp1b')
    %% MATLAB STUFFS
    clf
    disp('------------------------')                                      ;%
    disp(['Run ',num2str(sim),': nel = ',num2str(numEl),''])              ;%
    disp('------------------------')                                      ;%
    %% NON-DIMENSIONAL CONSTANT
    ar      = 0.8                                                         ;% aspect ratio thickness/width
    nu      = 0.3                                                         ;% poisson ratio
    ni      = 2                                                           ;% number of mp in h(1) direction
    nstr    = 3                                                           ;% number of stress components
    %---------------------------------------------------------------------%
    %% INDEPENDANT PHYSICAL CONSTANT
    g       = 9.81                                                        ;% gravitationnal acceleration [m/s^2]
    E       = 0.7e6                                                       ;% young's modulus             [Pa]
    Gc      = E/(2*(1+nu))                                                ;% shear modulus               [Pa]
    Kc      = E/(3*(1-2*nu))                                              ;% bulk modulus                [Pa]
    rho0    = 2650                                                        ;% density                     [kg/m^3]
    n0      = 0.0                                                         ;% initial porosity            [-]
    yd      = sqrt((Kc+4/3*Gc)/rho0)                                      ;% elastic wave speed          [m/s]
    coh0    = 0.0e3                                                       ;% cohesion                    [Pa]
    phi0    = 19.8*pi/180                                                 ;% friction angle              [Rad]
    psi0    =  0.0*pi/180                                                 ;%
    H       = -0.0e3                                                      ;% softening modulus           [Pa]
    cohr    =  0.0e3                                                      ;% residual cohesion           [Pa]
    phir    = 7.0*pi/180                                                  ;% residual friction angle     [Rad]
    t       = 1.00                                                        ;% simulation time             [s]
    te      = 0.001                                                       ;% elastic loading             [s]
    tg      = te/1.5                                                      ;% gravity increase            [s]
    %---------------------------------------------------------------------%
    %% MESH & MP INITIALIZATION
    colWidth  = 200e-3;
    colHeight = colWidth*lambda0(sim)*0.5;
    lz        = 1.1*colHeight;
    colElem   = numEl(sim);
    [meD,bc]  = meSetupExp1b(typeD,colHeight,colWidth,colElem);
    [mpD]     = mpSetupExp1b(meD,ni,colHeight,colWidth,coh0,cohr,phi0,phir,n0,rho0,nstr,typeD);% - see function
    mpD.coh   = ones(mpD.n,1)*coh0;
    % ISOTROPIC ELASTIC MATRIX
    Del       = [ Kc+4/3*Gc,Kc-2/3*Gc,Kc-2/3*Gc,0.0;...
                  Kc-2/3*Gc,Kc+4/3*Gc,Kc-2/3*Gc,0.0;...
                  Kc-2/3*Gc,Kc-2/3*Gc,Kc+4/3*Gc,0.0;...
                  0.0      ,0.0      ,0.0      ,Gc]                       ;%                
    Hp        = H*meD.h(1)                                                ;%
    %---------------------------------------------------------------------%
    %% MPM DM ALGORITHM EXPLICIT SOLVER FOR INFINITESIMAL STRAIN
    run('exportExp1b');
    % GPU CODE
    if(ismac || isunix)
        system('./gpu.out');
    elseif(ispc)
        system('gpu.exe');
    else
        disp('Outer space OS');
    end
    %---------------------------------------------------------------------%    
    %% DISPLAY
    fid  = fopen('GPUinfo.dat'); GPUtime = fread(fid,typeD); fclose(fid);
    fid  = fopen(['sig.dat']);   sig     = fread(fid,typeD); fclose(fid);
    fid  = fopen(['epII.dat']);  epII    = fread(fid,typeD); fclose(fid);
    fid  = fopen(['xp.dat']);    xp      = fread(fid,typeD); fclose(fid);
    fid  = fopen(['up.dat']);    up      = fread(fid,typeD); fclose(fid);
    s    = reshape(sig,mpD.n,4);
    epII = epII;
    x    = reshape(xp,mpD.n,2);
    u    = reshape(up,mpD.n,2);
    unorm= sqrt(u(:,1).^2+u(:,2).^2);
    p0   = -((s(:,1)+s(:,2)+s(:,3))./3)./1e3;
    
    run('postprocessing');
    
    name=['Exp1a_D_' num2str(mpD.n) 'np.mat'];
    save(name,'lz','mpD','meD','s','epII','x','u','-v7.3');
end
delete('*.dat','*.out','*.exe','*.avi','*.lib','*.exp');






