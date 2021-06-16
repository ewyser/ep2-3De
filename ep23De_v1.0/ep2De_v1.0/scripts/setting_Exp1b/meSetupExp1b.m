function [meD,bc] = meSetupBui(typeD,columnHeight,columnWidth,columnElem)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% MESH INITIALIZATION
dx      = (columnWidth)/(columnElem);
[xn,yn] = meshgrid(0-2*dx:dx:4*columnWidth+2*dx,0.0-2*dx:dx:columnHeight+3*dx)               ;%
meD.h   = [dx,dx]        ;                 
meD.L   = [max(xn)-min(xn),max(yn)-min(yn)]     ; 
               
xn      = flip(xn)                                                        ;%
yn      = flip(yn)                                                        ;%
               
               
meD.nnx = size(xn,2)                                                      ;% number of nodes along x
meD.nny = size(yn,1)                                                      ;% number of nodes along y
meD.no  = meD.nnx*meD.nny                                                 ;% number of nodes
meD.nex = meD.nnx-1;
meD.ney = meD.nny-1;
meD.nel = meD.nex*meD.ney                                                 ;%
meD.nn  = 16                                                              ;% number of node per element
meD.DoF = 2                                                               ;% degree of freedom
meD.nDoF = meD.DoF.*[meD.nn meD.no]                                       ;% local and global number of degree of freedom
meD.x   = xn(:)                                                           ;% x coordinate
meD.y   = yn(:)                                                           ;% y coordinate
% NODAL VECTOR INITIALIZATION                                              
meD.m   = zeros(meD.no    ,1,typeD)                                       ;% mass
meD.mr  = zeros(meD.nDoF(2),1,typeD)                                      ;% repmated mass
meD.f   = zeros(meD.nDoF(2),1,typeD)                                      ;% force balance                                                       
meD.fi  = zeros(meD.nDoF(2),1,typeD)                                      ;% internal force
meD.d   = zeros(meD.nDoF(2),1,typeD)                                      ;% damping force
meD.a   = zeros(meD.nDoF(2),1,typeD)                                      ;% acceleration
meD.p   = zeros(meD.nDoF(2),1,typeD)                                      ;% momentum
meD.v   = zeros(meD.nDoF(2),1,typeD)                                      ;% velocity
meD.u   = zeros(meD.nDoF(2),1,typeD)                                      ;% displacement
%--------------------------------------------------------------------------%
%% ELEMENT-NODE CONNECTIVITY
[meD.e2N] = e2N(meD.nny,meD.nnx,meD.nex,meD.ney,meD.nn)                   ;% element to node topology
%--------------------------------------------------------------------------%
% BOUNDARY CONDITIONS
meD.xB = [0 columnWidth 0.0 Inf]                                          ;%
[row]  = find(meD.y<=meD.xB(3))                                           ;%
bc.y   = row                                                              ;%
[row]  = find(meD.x<=meD.xB(1))                                           ;%
bc.x1  = row                                                              ;%
[row]  = find(meD.x>=meD.xB(2))                                           ;%
bc.x2  = row                                                              ;%
                                                
end
function [g_num] = e2N(nny,nnx,nelx,nely,nnpe)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
gnumbers= flip(reshape(1:(nnx*nny),nny ,nnx ))                            ;%
g_num   = zeros(nelx*nely,nnpe,'int32')                                   ;%
iel     = 1                                                               ;%
disp('------------------------')                                          ;%
disp('16 nodes per quadrilateral element')                                ;%
disp('------------------------')                                          ;%
for i = 1:nelx
    for j = 1:nely
        if(i>1 && i<nelx && j>1 && j<nely)
            g_num(iel,1 ) = gnumbers(j-1,i-1)                             ;%
            g_num(iel,2 ) = gnumbers(j-0,i-1)                             ;%
            g_num(iel,3 ) = gnumbers(j+1,i-1)                             ;%
            g_num(iel,4 ) = gnumbers(j+2,i-1)                             ;%
            
            g_num(iel,5 ) = gnumbers(j-1,i  )                             ;%
            g_num(iel,6 ) = gnumbers(j-0,i  )                             ;%
            g_num(iel,7 ) = gnumbers(j+1,i  )                             ;%
            g_num(iel,8 ) = gnumbers(j+2,i  )                             ;%
            
            g_num(iel,9 ) = gnumbers(j-1,i+1)                             ;%
            g_num(iel,10) = gnumbers(j-0,i+1)                             ;%
            g_num(iel,11) = gnumbers(j+1,i+1)                             ;%
            g_num(iel,12) = gnumbers(j+2,i+1)                             ;%
                
            g_num(iel,13) = gnumbers(j-1,i+2)                             ;%
            g_num(iel,14) = gnumbers(j-0,i+2)                             ;%
            g_num(iel,15) = gnumbers(j+1,i+2)                             ;%
            g_num(iel,16) = gnumbers(j+2,i+2)                             ;%
        end
        iel = iel+1;
    end
end
end
