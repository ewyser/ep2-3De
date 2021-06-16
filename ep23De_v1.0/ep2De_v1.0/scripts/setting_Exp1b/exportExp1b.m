%% EXPORT
%NODES
xn  = [meD.x,meD.y];
fid = fopen( 'xn.dat' ,'wb'); fwrite(fid, xn(:),typeD); fclose(fid);
fid = fopen( 'e2n.dat','wb'); fwrite(fid, meD.e2N(:)-1,'int32');  fclose(fid);
fid = fopen( 'mp.dat' ,'wb'); fwrite(fid, mpD.m,typeD); fclose(fid);
%PHYSICS
fid = fopen( 'Del.dat','wb'); fwrite(fid, Del(:),typeD); fclose(fid);
iDel = inv(Del);
fid = fopen( 'iDel.dat','wb'); fwrite(fid, iDel(:),typeD); fclose(fid);
fid = fopen( 'cohp.dat','wb'); fwrite(fid, mpD.coh(:),typeD); fclose(fid);
fid = fopen( 'phip.dat','wb'); fwrite(fid, mpD.phi(:),typeD); fclose(fid);
fid = fopen( 'epII.dat','wb'); fwrite(fid, mpD.epII(:),typeD); fclose(fid);
%PARAMETERS
vpx = mpD.v(:,1);
vpy = mpD.v(:,2);
xp  = mpD.x;
sig = mpD.s;
vol = mpD.V;
lp  = mpD.l;

p   = [mpD.n;meD.nn;meD.no;meD.h(:);min(meD.x);min(meD.y);meD.nnx;meD.nny];
fid = fopen( 'param.dat','wb'); fwrite(fid,p ,typeD)     ; fclose(fid);
p   = [g;rho0;psi0;nu;E;Kc;Gc;cohr;Hp;t;te;tg];
fid = fopen( 'phys.dat','wb') ; fwrite(fid,p ,typeD)     ; fclose(fid);
fid = fopen( 'mp.dat' ,'wb')  ; fwrite(fid, mpD.m ,typeD); fclose(fid);
fid = fopen( 'xp.dat','wb')   ; fwrite(fid, xp(:),typeD) ; fclose(fid);
fid = fopen( 'sig.dat','wb')  ; fwrite(fid, sig(:),typeD); fclose(fid);
fid = fopen( 'vol.dat','wb')  ; fwrite(fid, vol(:),typeD); fclose(fid);
fid = fopen( 'lp.dat','wb')   ; fwrite(fid, lp(:),typeD); fclose(fid);


bcx1 = int32(bc.x1)+0*meD.no;
bcx2 = int32(bc.x1)+0*meD.no;
bcy1 = int32(bc.y)+0*meD.no;
bcy2 = int32(bc.y)+1*meD.no;
BC  = ones(meD.no*2,1,'int32');
BC([bcx1;bcx2;bcy1;bcy2]) = 0;
fid = fopen( 'bcs.dat','wb')  ; fwrite(fid, BC(:),'int32'); fclose(fid);
BC  = ones(meD.no*2,1,'int32');
BC([bcx1;bcy1;bcy2]) = 0;
fid = fopen( 'bcN.dat','wb')  ; fwrite(fid, BC(:),'int32'); fclose(fid);