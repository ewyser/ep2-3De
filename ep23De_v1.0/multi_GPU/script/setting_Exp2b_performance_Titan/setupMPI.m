function [h] = setupMPI(nGPU,nEx,lx,ly,lz,typeD,ni,coh0,cohr,phi0,phir,n0,rho0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
L       = linspace(0,ly,nGPU+1)                                           ;% GPU parametric space (linear 1D)
l       = L(1:nGPU)                                                       ;% lower limit
u       = l+ly/nGPU                                                       ;% upper limit
nEz     = nEx                                                             ;% number of element in z direction
nEy     = nEx                                                             ;% number of element in y direction
meD.L   = [lx ceil(lz) ly]                                                      ;% mesh dimension [m]
meD.h   = [meD.L(1)/nEx meD.L(1)/nEx meD.L(1)/nEx]                        ;% [dx dy dz]
for k=1:nGPU
    disp(['GPU ',num2str(k),' / ',num2str(nGPU),': generating.'])         ;%
    %% mesh initialization
    meD.lim = [l(k),u(k)]                                                 ;%
    [X,Z,Y] = meshgrid(0.0-4*meD.h(1):meD.h(1):meD.L(1)+4*meD.h(1),...
                       0.0-4*meD.h(3):meD.h(3):meD.L(2)+4*meD.h(3),...
                      l(k)-4*meD.h(2):meD.h(2):u(k)    +4*meD.h(2))       ;%
    x       = flip(X)                                                     ;%
    z       = flip(Z)                                                     ;%
    y       = flip(Y)                                                     ;%
    meD.nnx = size(x,2)                                                   ;% number of nodes along x
    meD.nnz = size(z,1)                                                   ;% number of nodes along z
    meD.nny = size(y,3)                                                   ;% number of nodes along z
    meD.no  = meD.nnx*meD.nnz*meD.nny                                     ;% number of nodes
    meD.nex = meD.nnx-1                                                   ;%
    meD.nez = meD.nnz-1                                                   ;%
    meD.ney = meD.nny-1                                                   ;%
    meD.nel = meD.nex*meD.nez*meD.ney                                     ;%
    meD.nn  = 64                                                          ;% number of node per element
    meD.x   = x(:)                                                        ;% x coordinate
    meD.z   = z(:)                                                        ;% y coordinate
    meD.y   = y(:)                                                        ;% z coordinate
    %----------------------------------------------------------------------%
    % element to nodes topology
    [meD.e2N] = e2N(meD.nnz,meD.nnx,meD.nny,meD.nex,meD.nez,meD.ney,meD.nn);% element to node topology
    %----------------------------------------------------------------------%
    % boundary condition
    meD.xB = [0.0,meD.L(1),0.0,Inf,0.0,meD.L(3)]                          ;%
    [row]  = find(meD.x<=meD.xB(1))                                       ;%
    bc.x   = row                                                          ;%
    [row]  = find(meD.x>=meD.xB(2))                                       ;%
    bc.x   = [bc.x;row]                                                   ;%
    [row]  = find(meD.z<=meD.xB(3))                                       ;%
    bc.z   = row                                                          ;%
    [row]  = find(meD.y<=meD.xB(5))                                       ;%
    bc.y   = row                                                          ;%
    [row]  = find(meD.y>=meD.xB(6))                                       ;%
    bc.y   = [bc.y;row]                                                   ;%
    %----------------------------------------------------------------------%
    %% mpm initialization
    xL      = meD.xB(1) +(0.5*meD.h(1)/ni):meD.h(1)/ni:meD.xB(2)   ;
    zL      = meD.xB(3) +(0.5*meD.h(2)/ni):meD.h(2)/ni:lz-(0.5*meD.h(2)/ni)       ;
    yL      = meD.lim(1)+(0.5*meD.h(1)/ni):meD.h(1)/ni:meD.lim(2)   ;
    [xl,zl,yl] = meshgrid(xL,zL,yL)                                                 ;
    wl      = 0.15*lz                                                         ; % 0.15
    
    
    
    xl  = xl(:);yl  = yl(:);zl  = zl(:);
    coh = coh0.*ones(size(xl),typeD);
    % % %run('GRFS_gauss');
% % run('GRFS_exp');
% % xl  = xl(:);
% % yl  = yl(:);
% % zl  = zl(:);
% % coh = coh(:);
    
    
    
    
    x=linspace(min(xl),max(xl),200);
    a= -1.25;
    z= a.*x;
    x= x+meD.L(1)/2;
    
    xlt = [];zlt = [];ylt = [];cot = [];
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
    xl = xlt;zl = zlt;yl = ylt;
    %----------------------------------------------------------------------%
    % mp's state variables
    mpD.n   = length(xl(:))                                               ;% number of material point
    mpD.n0  = ones(mpD.n,1,typeD).*n0                                     ;% porosity
    mpD.l   = ones(mpD.n,3,typeD).*(meD.h(1)/ni)./2                       ;% reference domain dimension
    mpD.V   = ones(mpD.n,1,typeD).*(2.*mpD.l(:,1).*2.*mpD.l(:,2).*2.*mpD.l(:,3)) ;% reference volume
    mpD.m   = rho0.*mpD.V                                                 ;% mass
    mpD.x   = [xl(:) yl(:) zl(:)]                                         ;% coordinate
    mpD.coh = cot(:)                                                      ;% cohesion
    mpD.cohr= cohr.*ones(mpD.n,1,typeD)                                   ;% cohesion
    mpD.phi = ones(1,mpD.n,typeD).*phi0                                   ;% friction
    mpD.phi(zl<2*wl)= phir                                                ;% -
    %----------------------------------------------------------------------%
    %% export
    disp(['GPU ',num2str(k),' / ',num2str(nGPU),': exporting'])         ;%
    % mesh, topology & boundary conditions
    x   = meD.x                                            ;% x-th component
    y   = meD.y                                            ;% y-th component
    z   = meD.z                                            ;% z-th component
    fid = fopen(['xn_rank_',num2str(k-1),'.dat'] ,'wb')           ;% export
    fwrite(fid, [x(:);y(:);z(:)],typeD)                           ;% export
    fclose(fid)                                                   ;% export
    e2n = meD.e2N                                          ;% element to nodes topology list
    fid = fopen(['e2n_rank_',num2str(k-1),'.dat'],'wb')           ;% export
    fwrite(fid, e2n(:)-1,'int32')                                   ;% export
    fclose(fid)                                                   ;% export
    Bx  = int32(bc.x)                                     ;% boundary number in x
    By  = int32(bc.y)                                     ;% boundary number in y
    Bz  = int32(bc.z)                                     ;% boundary number in z
    Bxx = int32(Bx)+0*meD.no                               ;% -
    Byy = int32(By)+1*meD.no                               ;% -
    Bzx = int32(Bz)+0*meD.no                               ;% -
    Bzy = int32(Bz)+1*meD.no                               ;% -
    Bzz = int32(Bz)+2*meD.no                               ;% -
    BC  = ones(meD.no*3,1,'int32')                         ;% -
    BC([Bxx;Byy;Bzz]) = 0                                             ;% logical vector of active bc's (=0)
    fid = fopen(['bcs_rank_',num2str(k-1),'.dat'],'wb')           ;% export
    fwrite(fid, BC(:),'int32')                                    ;% export
    fclose(fid)                                                   ;% export
    %--------------------------------------------------------------------------%
    % material points
    nmp  = mpD.n                                                       ;
    xp  = mpD.x                                                 ;% coordinates
    fid = fopen(['xp_rank_',num2str(k-1),'.dat'] ,'wb')           ;% export
    fwrite(fid, xp(:),typeD)                                      ;% export
    fclose(fid)                                                   ;% export
    mp  = mpD.m                                                   ;% mass
    fid = fopen(['mp_rank_',num2str(k-1),'.dat'] ,'wb')           ;% export
    fwrite(fid, mp(:),typeD)                                      ;% export
    fclose(fid)                                                   ;% export
    lp  = mpD.l                                                 ;% domain size
    fid = fopen(['lp_rank_',num2str(k-1),'.dat'] ,'wb')           ;% export
    fwrite(fid, lp(:),typeD)                                      ;% export
    fclose(fid)                                                   ;% export
    vol = mpD.V                                                 ;% initial volume
    fid = fopen(['vol_rank_',num2str(k-1),'.dat'],'wb')           ;% export
    fwrite(fid, vol(:),typeD)                                  ;% export
    fclose(fid)                                                   ;% export
    coh = mpD.coh                                               ;% cohesion
    fid = fopen(['cohp_rank_',num2str(k-1),'.dat'],'wb')          ;% export
    fwrite(fid,coh(:),typeD)                              ;% export
    fclose(fid)                                                   ;% export
    phi = mpD.phi                                               ;% friction angle
    fid = fopen(['phip_rank_',num2str(k-1),'.dat'],'wb')          ;% export
    fwrite(fid, mpD.phi(:),typeD)                              ;% export
    fclose(fid)                                                   ;% export
    %--------------------------------------------------------------------------%
    % 
% %     figure(4352)
% %     plot3(x,y,z,'s',xp(:,1),xp(:,2),xp(:,3),'x')
% %     axis equal
    
    % numerical parameters
    p = [nmp;...
        meD.no;...
        meD.nex*meD.ney*meD.nez;...
        meD.nex;...
        meD.ney;...
        meD.nez;...
        meD.nn;...
        meD.h(1);...
        meD.h(2);...
        meD.h(3);...
        min(x(:));...
        min(y(:));...
        min(z(:))];
    fid = fopen( ['param_rank_',num2str(k-1),'.dat'],'wb')        ;%
    fwrite(fid,p(:),typeD)                                        ;%
    fclose(fid)                                                   ;%
    disp(['GPU ',num2str(k),' / ',num2str(nGPU),': done'])         ;%
end
h = meD.h;
end
function [g_num] = e2N(nnz,nnx,nny,nelx,nelz,nely,nnpe)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
gnumbers= flip(reshape(1:(nnx*nny*nnz),nnz,nnx,nny))                            ;%
g_num   = zeros(nelx*nelz*nely,nnpe,'int32')                                   ;%
iel     = 1                                                               ;%
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
