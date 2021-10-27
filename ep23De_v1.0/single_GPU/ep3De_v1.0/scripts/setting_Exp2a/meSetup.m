function [meD,bc] = meSetup(nEx,lx,ly,lz,typeD)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% MESH INITIALIZATION
meD.vd  = 0.1                                                             ;% viscous damping coefficient
Lx      = lx                                                              ;% mesh dimension (x,y)
Ly      = ly                                                              ;% mesh dimension (x,y)
Lz      = ceil(lz)                                                        ;% number of element in x direction
nEz     = nEx                                                             ;% number of element in z direction
nEy     = nEx                                                             ;% number of element in y direction
meD.L   = [Lx Lz Ly]                                                      ;% mesh length in x            [m]
meD.h   = [meD.L(1)/nEx meD.L(1)/nEx meD.L(1)/nEx]                        ;% [dx dy dz]
[x,z,y] = meshgrid(0.0-2*meD.h(1):meD.h(1):meD.L(1)+2*meD.h(1),...
                   0.0-2*meD.h(3):meD.h(3):meD.L(2)+2*meD.h(3),...
                   0.0-2*meD.h(2):meD.h(2):meD.L(3)+2*meD.h(2))           ;%
x      = flip(x)                                                          ;%
z      = flip(z)                                                          ;%
y      = flip(y)                                                          ;%               
               
meD.nnx = size(x,2)                                                       ;% number of nodes along x
meD.nnz = size(z,1)                                                       ;% number of nodes along z
meD.nny = size(y,3)                                                       ;% number of nodes along z
meD.no  = meD.nnx*meD.nnz*meD.nny                                         ;% number of nodes
meD.nex = meD.nnx-1;
meD.nez = meD.nnz-1;
meD.ney = meD.nny-1;
meD.nel = meD.nex*meD.nez*meD.ney                                         ;%
meD.nn  = 64                                                              ;% number of node per element
meD.DoF = 3                                                               ;% degree of freedom
meD.nDoF = meD.DoF.*[meD.nn meD.no]                                       ;% local and global number of degree of freedom
meD.x   = x(:)                                                            ;% x coordinate
meD.z   = z(:)                                                            ;% y coordinate
meD.y   = y(:)                                                            ;% z coordinate
% NODAL VECTOR INITIALIZATION                                              
meD.m   = zeros(meD.no     ,1,typeD)                                      ;% mass
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
[meD.e2N] = e2N(meD.nnz,meD.nnx,meD.nny,meD.nex,meD.nez,meD.ney,meD.nn)   ;% element to node topology
%--------------------------------------------------------------------------%
% BOUNDARY CONDITIONS
meD.xB = [min(meD.x)+2*meD.h(1) max(meD.x)-2*meD.h(1) 0.0 Inf min(meD.y)+2*meD.h(1) max(meD.y)-2*meD.h(1)]            ;%

[row]  = find(meD.x<=meD.xB(1))                                           ;%
bc.x   = row                                                              ;%
[row]  = find(meD.x>=meD.xB(2))                                           ;%
bc.x   = [bc.x;row]                                                       ;%

[row]  = find(meD.z<=meD.xB(3))                                           ;%
bc.z   = row                                                              ;%

[row]  = find(meD.y<=meD.xB(5))                                           ;%
bc.y   = row                                                              ;%
[row]  = find(meD.y>=meD.xB(6))                                           ;%
bc.y   = [bc.y;row]                                                       ;%
end
function [g_num] = e2N(nnz,nnx,nny,nelx,nelz,nely,nnpe)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
gnumbers= flip(reshape(1:(nnx*nny*nnz),nnz,nnx,nny))                            ;%
g_num   = zeros(nelx*nelz*nely,nnpe,'int32')                                   ;%
iel     = 1                                                               ;%
disp('------------------------')                                          ;%
disp('16 nodes per quadrilateral element')                                ;%
disp('------------------------')                                          ;%
    for k = 1:nely
        for i = 1:nelx
            for j = 1:nelz
                if(i>1 && i<nelx && j>1 && j<nelz && k>1 && k<nely )
                    g_num(iel,1 )= gnumbers(j-1,i-1,k-1);
                    g_num(iel,2 )= gnumbers(j-0,i-1,k-1);
                    g_num(iel,3 )= gnumbers(j+1,i-1,k-1);
                    g_num(iel,4 )= gnumbers(j+2,i-1,k-1);
                    g_num(iel,5 )= gnumbers(j-1,i  ,k-1);
                    g_num(iel,6 )= gnumbers(j-0,i  ,k-1);
                    g_num(iel,7 )= gnumbers(j+1,i  ,k-1);
                    g_num(iel,8 )= gnumbers(j+2,i  ,k-1);
                    g_num(iel,9 )= gnumbers(j-1,i+1,k-1);
                    g_num(iel,10)= gnumbers(j-0,i+1,k-1);
                    g_num(iel,11)= gnumbers(j+1,i+1,k-1);
                    g_num(iel,12)= gnumbers(j+2,i+1,k-1);
                    g_num(iel,13)= gnumbers(j-1,i+2,k-1);
                    g_num(iel,14)= gnumbers(j-0,i+2,k-1);
                    g_num(iel,15)= gnumbers(j+1,i+2,k-1);
                    g_num(iel,16)= gnumbers(j+2,i+2,k-1);
                    
                    g_num(iel,17)= gnumbers(j-1,i-1,k  );
                    g_num(iel,18)= gnumbers(j-0,i-1,k  );
                    g_num(iel,19)= gnumbers(j+1,i-1,k  );
                    g_num(iel,20)= gnumbers(j+2,i-1,k  );
                    g_num(iel,21)= gnumbers(j-1,i  ,k  );
                    g_num(iel,22)= gnumbers(j-0,i  ,k  );
                    g_num(iel,23)= gnumbers(j+1,i  ,k  );
                    g_num(iel,24)= gnumbers(j+2,i  ,k  );
                    g_num(iel,25)= gnumbers(j-1,i+1,k  );
                    g_num(iel,26)= gnumbers(j-0,i+1,k  );
                    g_num(iel,27)= gnumbers(j+1,i+1,k  );
                    g_num(iel,28)= gnumbers(j+2,i+1,k  );
                    g_num(iel,29)= gnumbers(j-1,i+2,k  );
                    g_num(iel,30)= gnumbers(j-0,i+2,k  );
                    g_num(iel,31)= gnumbers(j+1,i+2,k  );
                    g_num(iel,32)= gnumbers(j+2,i+2,k  );
                    
                    g_num(iel,33)= gnumbers(j-1,i-1,k+1);
                    g_num(iel,34)= gnumbers(j-0,i-1,k+1);
                    g_num(iel,35)= gnumbers(j+1,i-1,k+1);
                    g_num(iel,36)= gnumbers(j+2,i-1,k+1);
                    g_num(iel,37)= gnumbers(j-1,i  ,k+1);
                    g_num(iel,38)= gnumbers(j-0,i  ,k+1);
                    g_num(iel,39)= gnumbers(j+1,i  ,k+1);
                    g_num(iel,40)= gnumbers(j+2,i  ,k+1);
                    g_num(iel,41)= gnumbers(j-1,i+1,k+1);
                    g_num(iel,42)= gnumbers(j-0,i+1,k+1);
                    g_num(iel,43)= gnumbers(j+1,i+1,k+1);
                    g_num(iel,44)= gnumbers(j+2,i+1,k+1);
                    g_num(iel,45)= gnumbers(j-1,i+2,k+1);
                    g_num(iel,46)= gnumbers(j-0,i+2,k+1);
                    g_num(iel,47)= gnumbers(j+1,i+2,k+1);
                    g_num(iel,48)= gnumbers(j+2,i+2,k+1);
                    
                    g_num(iel,49)= gnumbers(j-1,i-1,k+2);
                    g_num(iel,50)= gnumbers(j-0,i-1,k+2);
                    g_num(iel,51)= gnumbers(j+1,i-1,k+2);
                    g_num(iel,52)= gnumbers(j+2,i-1,k+2);
                    g_num(iel,53)= gnumbers(j-1,i  ,k+2);
                    g_num(iel,54)= gnumbers(j-0,i  ,k+2);
                    g_num(iel,55)= gnumbers(j+1,i  ,k+2);
                    g_num(iel,56)= gnumbers(j+2,i  ,k+2);
                    g_num(iel,57)= gnumbers(j-1,i+1,k+2);
                    g_num(iel,58)= gnumbers(j-0,i+1,k+2);
                    g_num(iel,59)= gnumbers(j+1,i+1,k+2);
                    g_num(iel,60)= gnumbers(j+2,i+1,k+2);
                    g_num(iel,61)= gnumbers(j-1,i+2,k+2);
                    g_num(iel,62)= gnumbers(j-0,i+2,k+2);
                    g_num(iel,63)= gnumbers(j+1,i+2,k+2);
                    g_num(iel,64)= gnumbers(j+2,i+2,k+2);
                end
                iel = iel+1;
            end
        end
        %disp(['element = ',num2str(iel/(nelx*nelz*nely)),'']);
    end


end
