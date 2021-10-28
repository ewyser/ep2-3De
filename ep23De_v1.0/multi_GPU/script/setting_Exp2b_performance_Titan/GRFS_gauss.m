% =====================================================================
% GRFS: Gaussian Random Field Simulator - Gaussian covariance
% 
% Copyright (C) 2019  Ludovic Raess, Dmitriy Kolyukhin and Alexander Minakov.
% 
% GRFS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% GRFS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with GRFS. If not, see <http://www.gnu.org/licenses/>.
% =====================================================================
reset(RandStream.getGlobalStream);
% physics
sf  = 5e3;                        % standard deviation
If  = 2.5;                        % correlation lengths in [x,y,z]
% numerics
Nh  = 5000;                       % inner parameter, number of harmonics
k_m = 100;                        % maximum value of the wave number
nx  = size(xl,2);                 % numerical grid resolution in x
ny  = size(xl,3);                 % numerical grid resolution in y
nz  = size(xl,1);                 % numerical grid resolution in z
dx  = 0.5*meD.h(1)/ni;            % numerical grid step size in x
dy  = 0.5*meD.h(1)/ni;            % numerical grid step size in y
dz  = 0.5*meD.h(1)/ni;            % numerical grid step size in z
% preprocessing
C   = sf/sqrt(Nh); 
coh = zeros(nz,nx,ny);
tmp = zeros(nz,nx,ny);
% action
for ih = 1:Nh
    fi = 2*pi*rand;
    lf = 2*If/sqrt(pi);
    %   Gaussian spectrum
    flag = true;
    while flag
        k = k_m*rand;
        d = k*k*exp(-0.5*k*k);
        if (rand*2*exp(-1)<d)
            flag = false;
        end
    end   
    k     = sqrt(2)*k/lf;    
    theta = acos(1-2*rand);
    V1 = k*sin(fi)*sin(theta);
    V2 = k*cos(fi)*sin(theta);
    V3 = k*cos(theta); 
    a  = randn;
    b  = randn;
    for iz=1:nz
    for iy=1:ny
    for ix=1:nx
        tmp(iz,ix,iy) = dx*(ix-0.5)*V1 + dy*(iy-0.5)*V2 + dz*(iz-0.5)*V3;
    end
    end
    end 
    coh = coh + a*sin(tmp) + b*cos(tmp);
end
coh = coh0+C*coh;coh(coh<cohr)=cohr;