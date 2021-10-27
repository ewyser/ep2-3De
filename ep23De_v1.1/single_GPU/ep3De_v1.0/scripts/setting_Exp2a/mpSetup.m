function [mpD] = mpSetup(meD,ni,lz,coh0,cohr,phi0,phir,n0,rho0,typeD)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% MPM DISCRETIZATION
% layer initialization
xL      = meD.xB(1)+(0.5*meD.h(1)/ni):meD.h(1)/ni:meD.xB(2)   ;
zL      = meD.xB(3)+(0.5*meD.h(2)/ni):meD.h(2)/ni:  lz-(0.5*meD.h(2)/ni)       ;
yL      = meD.xB(5)+(0.5*meD.h(1)/ni):meD.h(1)/ni:meD.xB(6)   ;
[xl,zl,yl] = meshgrid(xL,zL,yL)                                                 ;
wl      = 0.15*lz                                                         ; % 0.15

xl  = xl(:);
yl  = yl(:);
zl  = zl(:);

x=linspace(min(xl),max(xl),200);
a= -1.25;
z= a.*x;
x= x+meD.L(1)/2;

xlt = [];
zlt = [];
ylt = [];
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
mpD.coh = coh0.*ones(mpD.n,1,typeD)                                       ;% cohesion
mpD.cohr= cohr.*ones(mpD.n,1,typeD)                                       ;% cohesion
mpD.phi = ones(1,mpD.n,typeD).*phi0                                       ;% friction 
mpD.phi(zl<2*wl)= phir                                                    ;% -


end

