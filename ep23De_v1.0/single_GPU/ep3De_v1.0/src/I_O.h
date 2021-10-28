// header I_O.h to import data and set initial conditions
// GEOMETRY
    load(param,12,1,"param.dat",DAT);
    int nmp   = (int) param_h[0 ];
    int nn    = (int) param_h[1 ];
    int no    = (int) param_h[2 ];
    DAT dx    = (DAT) param_h[3 ];
    DAT dy    = (DAT) param_h[4 ];
    DAT dz    = (DAT) param_h[5 ];
    DAT xnmin = (DAT) param_h[6 ];
    DAT ynmin = (DAT) param_h[7 ];
    DAT znmin = (DAT) param_h[8 ];
    int nnx   = (int) param_h[9 ];
    int nny   = (int) param_h[10];
    int nnz   = (int) param_h[11];
    int nex   = (int) nnx-1;
    int ney   = (int) nny-1;
    int nez   = (int) nnz-1;
    int nel   = nex*ney*nez;
// PHYSICS
    load(phys,12,1,"phys.dat",DAT);
    DAT g     = (DAT) phys_h[0 ];
    DAT rho0  = (DAT) phys_h[1 ];
    DAT psi0  = (DAT) phys_h[2 ];
    DAT nu    = (DAT) phys_h[3 ];
    DAT E     = (DAT) phys_h[4 ];
    DAT Kc    = (DAT) phys_h[5 ];
    DAT Gc    = (DAT) phys_h[6 ];
    DAT cohr  = (DAT) phys_h[7 ];
    DAT Hp    = (DAT) phys_h[8 ];
    DAT t     = (DAT) phys_h[9 ];
    DAT te    = (DAT) phys_h[10];
    DAT tg    = (DAT) phys_h[11];
    DAT dt    = (DAT) 0.0;
    DAT yd    = (DAT) sqrt((Kc+1.333*Gc)*(DAT)1.0/(DAT)rho0);
    DAT tw    = (DAT) 0.0;
    int it    = (int) 1;
// IMPORT DATA
    load(mp  ,nmp       ,1,"mp.dat"  ,DAT);
    load(e2n ,nel*nn    ,1,"e2n.dat" ,int);
    load(xn  ,no*3      ,1,"xn.dat"  ,DAT);
    load(xp  ,3*nmp     ,1,"xp.dat"  ,DAT);
    load(vol ,nmp       ,1,"vol.dat" ,DAT);
    load(lp  ,nmp*3     ,1,"lp.dat"  ,DAT);
    load(Del ,36        ,1,"Del.dat" ,DAT);
    load(iDel,36        ,1,"iDel.dat",DAT);
    load(cohp,nmp       ,1,"cohp.dat",DAT);
    load(phip,nmp       ,1,"phip.dat",DAT);
// INITIALIZE
// elements
    zeros(pel ,2*nel    ,1           ,DAT);
// nodes
    zeros(mn  ,no       ,1           ,DAT);
    zeros(pn  ,3*no     ,1           ,DAT);
    zeros(fen ,3*no     ,1           ,DAT);
    zeros(fin ,3*no     ,1           ,DAT);
    zeros(fn  ,3*no     ,1           ,DAT);
    zeros(an  ,3*no     ,1           ,DAT);
    zeros(vn  ,3*no     ,1           ,DAT);
    zeros(un  ,3*no     ,1           ,DAT);
// material points
    zeros(p2e ,nmp      ,1           ,int);
    zeros(p2n ,nmp*nn   ,1           ,int);
    zeros(sig ,nmp*6    ,1           ,DAT);
    zeros(epII,nmp      ,1           ,DAT);
    zeros(N   ,nmp*nn   ,1           ,DAT);
    zeros(dNx ,nmp*nn   ,1           ,DAT);
    zeros(dNy ,nmp*nn   ,1           ,DAT);
    zeros(dNz ,nmp*nn   ,1           ,DAT);
    zeros(vp  ,nmp*3    ,1           ,DAT);
    zeros(up  ,nmp*3    ,1           ,DAT);
    zeros(dF  ,nmp*9    ,1           ,DAT);
    zeros(eps ,nmp*6    ,1           ,DAT);
    zeros(ome ,nmp*3    ,1           ,DAT);
    zeros(dev ,nmp*6    ,1           ,DAT);    