function [mpD] = mpSetupGRFS(meD,ni,lz,coh0,cohr,phi0,phir,n0,rho0,typeD)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
set(0,'defaulttextinterpreter','latex')                                   ;%
%% MPM DISCRETIZATION
% layer initialization
xL      = meD.xB(1)+(0.5*meD.h(1)/ni):meD.h(1)/ni:meD.xB(2)   ;
zL      = meD.xB(3)+(0.5*meD.h(2)/ni):meD.h(2)/ni:  lz-(0.5*meD.h(2)/ni)       ;
yL      = meD.xB(5)+(0.5*meD.h(1)/ni):meD.h(1)/ni:meD.xB(6)   ;
[xl,zl,yl] = meshgrid(xL,zL,yL)                                                 ;
wl      = 0.15*lz                                                         ; % 0.15

run('GRFS_gauss');
%run('GRFS_exp');
xl  = xl(:);
yl  = yl(:);
zl  = zl(:);
coh = coh(:);

x=linspace(min(xl),max(xl),200);
a= -1.25;
z= a.*x;
x= x+meD.L(1)/2;

xlt = [];
zlt = [];
ylt = [];
cot = [];
for mp=1:length(xl)
    for p = 1:length(z)
        DX = xl(mp)-x(p);
        DZ = zl(mp)-z(p);
        nx = a;
        ny = -1;
        s = DX*nx+DZ*ny;
        if(s>0)
            pos = 1;
        else
            pos = 0;
        end
        if(zl(mp)<wl)
            pos = 1;
        end
    end
    if(pos==1)
        xlt = [xlt xl(mp)];
        zlt = [zlt zl(mp)];
        ylt = [ylt yl(mp)];
        cot = [cot coh(mp)];
    end
end
xl = xlt;
zl = zlt;
yl = ylt;
%% MATERIAL POINT:
% SCALAR & VECTOR
mpD.n   = length(xl(:))                                                   ;% number of material point
mpD.n0  = ones(mpD.n,1,typeD).*n0                                         ;% porosity
mpD.l   = ones(mpD.n,3,typeD).*(meD.h(1)/ni)./2                           ;% reference domain dimension
mpD.V   = ones(mpD.n,1,typeD).*(2.*mpD.l(:,1).*2.*mpD.l(:,2).*2.*mpD.l(:,3)) ;% reference volume
mpD.m   = rho0.*mpD.V                                                     ;% mass
mpD.x   = [xl(:) yl(:) zl(:)]                                             ;% coordinate
mpD.coh = cot(:)                                                          ;% cohesion
mpD.cohr= cohr.*ones(mpD.n,1,typeD)                                       ;% cohesion
mpD.phi = ones(1,mpD.n,typeD).*phi0                                       ;% friction
mpD.phi(zl<2*wl)= phir                                                    ;% -
%% DISPLAY
try
    fig1 = figure(1);
    a = ones(mpD.n,1); a(mpD.coh<coh0+sf & mpD.coh>coh0-sf) = 0.0;
    set(fig1,'Units','pixels','Position',[681 761 577 218]);
    colormap(parula(1000));
    caxis=([round(min(mpD.coh(:)/1e3)) fix(max(mpD.coh(:)/1e3))]);
    s = scatter3(mpD.x(:,1),mpD.x(:,2),mpD.x(:,3),10.0,mpD.coh(:)/1e3,'filled');
    s.MarkerFaceAlpha = 'flat';
    s.AlphaData       = a;
    s.AlphaDataMapping= 'none';
    cb=colorbar('FontSize',12,'TickLabelInterpreter','latex','Fontsize',12,'Location','eastoutside');
    cb.Position         =[0.85 0.4 0.0370 0.36];
    cb.Label.String     ='$c_{\mathrm{peak}}$ [kPa]';
    cb.Label.Units      ='normalized';
    cb.Label.Position   =[-1 0.5];
    cb.Label.FontSize   =12;
    cb.Label.Interpreter='Latex';
    cb.Ticks            =[round(min(mpD.coh(:)/1e3)) fix(max(mpD.coh(:)/1e3))];
    axis equal;
    xlabel('$x$ [m]');
    ylabel('$y$ [m]');
    zlabel('$z$ [m]');
    yticks([0 5 10]);
    zticks([0 5 10]);
    zlim([0 lz]);
    set(gca,'FontSize',12,'TickLabelInterpreter','latex');
    box on;
    ax = gca;
    ax.BoxStyle = 'full';
    ax.Position = [[0.05 0.20 0.7750 0.7234]];
    view(47,13);
    name = strcat('GRFS_cohesion_pmsigma');
    print(fig1,name,'-dpng','-r600');
    
    fig2 = figure(2);
    a = ones(mpD.n,1);
    set(fig2,'Units','pixels','Position',[681 761 577 218]);
    colormap(parula(1000));
    caxis=([round(min(mpD.coh(:)/1e3)) fix(max(mpD.coh(:)/1e3))]);
    s = scatter3(mpD.x(:,1),mpD.x(:,2),mpD.x(:,3),10.0,mpD.coh(:)/1e3,'filled');
    s.MarkerFaceAlpha = 'flat';
    s.AlphaData       = a;
    s.AlphaDataMapping= 'none';
    cb=colorbar('FontSize',12,'TickLabelInterpreter','latex','Fontsize',12,'Location','eastoutside');
    cb.Position         =[0.85 0.4 0.0370 0.36];
    cb.Label.String     ='$c_{\mathrm{peak}}$ [kPa]';
    cb.Label.Units      ='normalized';
    cb.Label.Position   =[-1 0.5];
    cb.Label.FontSize   =12;
    cb.Label.Interpreter='Latex';
    cb.Ticks            =[round(min(mpD.coh(:)/1e3)) fix(max(mpD.coh(:)/1e3))];
    axis equal;
    xlabel('$x$ [m]');
    ylabel('$y$ [m]');
    zlabel('$z$ [m]');
    yticks([0 5 10]);
    zticks([0 5 10]);
    zlim([0 lz]);
    set(gca,'FontSize',12,'TickLabelInterpreter','latex');
    box on;
    ax = gca;
    ax.BoxStyle = 'full';
    ax.Position = [[0.05 0.20 0.7750 0.7234]];
    view(47,13);
    name = strcat('GRFS_cohesion_all');
    print(fig2,name,'-dpng','-r600');
catch
end
end

