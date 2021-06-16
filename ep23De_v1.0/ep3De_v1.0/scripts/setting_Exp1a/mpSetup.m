function [mpD,p2] = mpSetup(meD,ni,ly,height,width,coh0,cohr,phi0,phir,n0,rho0,typeD)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% MPM DISCRETIZATION
% layer initialization
xL      = (0.0)+(0.5*meD.h(1)/ni):meD.h(1)/ni:(0.0+width);
zL      = meD.xB(3)+(0.5*meD.h(2)/ni):meD.h(2)/ni:height-(0.5*meD.h(2)/ni)       ;
yL      = (0.0)+(0.5*meD.h(2)/ni):meD.h(2)/ni:(ly-(0.5*meD.h(2)/ni));



[xl,zl,yl] = meshgrid(xL,zL,yL)                                               ;
xl = xl(:);
yl = yl(:);

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
mpD.phi = ones(1,mpD.n,typeD).*phi0                                       ;% friction angle 
end

