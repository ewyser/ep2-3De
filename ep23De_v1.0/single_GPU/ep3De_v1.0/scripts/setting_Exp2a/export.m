%% EXPORT
%NODES
xn  = [meD.x,meD.y,meD.z];
fid = fopen( 'xn.dat' ,'wb'); fwrite(fid, xn(:),typeD); fclose(fid);
fid = fopen( 'e2n.dat','wb'); fwrite(fid, meD.e2N(:)-1,'int32');  fclose(fid);
fid = fopen( 'mp.dat' ,'wb'); fwrite(fid, mpD.m,typeD); fclose(fid);
%PHYSICS
fid = fopen( 'Del.dat','wb'); fwrite(fid, Del(:),typeD); fclose(fid);
iDel = inv(Del);
fid = fopen( 'iDel.dat','wb'); fwrite(fid, iDel(:),typeD); fclose(fid);
fid = fopen( 'cohp.dat','wb'); fwrite(fid, mpD.coh(:),typeD); fclose(fid);
fid = fopen( 'phip.dat','wb'); fwrite(fid, mpD.phi(:),typeD); fclose(fid);
%PARAMETERS
xp  = mpD.x;
vol = mpD.V;
lp  = mpD.l;

p   = [mpD.n;meD.nn;meD.no;meD.h(:);min(meD.x);min(meD.y);min(meD.z);meD.nnx;meD.nny;meD.nnz];
fid = fopen( 'param.dat','wb'); fwrite(fid,p ,typeD)     ; fclose(fid);
p   = [g;rho0;psi0;nu;E;Kc;Gc;cohr;Hp;t;te;tg];
fid = fopen( 'phys.dat','wb') ; fwrite(fid,p ,typeD)     ; fclose(fid);
fid = fopen( 'mp.dat' ,'wb')  ; fwrite(fid, mpD.m ,typeD); fclose(fid);
fid = fopen( 'xp.dat','wb')   ; fwrite(fid, xp(:),typeD) ; fclose(fid);
fid = fopen( 'vol.dat','wb')  ; fwrite(fid, vol(:),typeD); fclose(fid);
fid = fopen( 'lp.dat','wb')   ; fwrite(fid, lp(:),typeD); fclose(fid);

bcx = int32(bc.x);
bcy = int32(bc.y);
bcz = int32(bc.z);
bcx1 = int32(bcx)+0*meD.no;
% bcx2 = int32(bcx)+1*meD.no;
% bcx3 = int32(bcx)+2*meD.no;
% bcy1 = int32(bcy)+0*meD.no;
bcy2 = int32(bcy)+1*meD.no;
% bcy3 = int32(bcy)+2*meD.no;
% bcz1 = int32(bcz)+0*meD.no;
% bcz2 = int32(bcz)+1*meD.no;
bcz3 = int32(bcz)+2*meD.no;
BC  = ones(meD.no*3,1,'int32');
BC([bcx1;bcy2;bcz3]) = 0;
fid = fopen( 'bcs.dat','wb')  ; fwrite(fid, BC(:),'int32'); fclose(fid);



bcx = int32(bc.x);
bcy = int32(bc.y);
bcz = int32(bc.z);
bcx1 = int32(bcx)+0*meD.no;
bcx2 = int32(bcx)+1*meD.no;
bcx3 = int32(bcx)+2*meD.no;
bcy1 = int32(bcy)+0*meD.no;
bcy2 = int32(bcy)+1*meD.no;
bcy3 = int32(bcy)+2*meD.no;
bcz3 = int32(bcz)+2*meD.no;
BC  = ones(meD.no*3,1,'int32');
BC([bcx1;bcx2;bcx3;bcy1;bcy2;bcy3;bcz3]) = 0;
fid = fopen( 'bcN.dat','wb')  ; fwrite(fid, BC(:),'int32'); fclose(fid);

