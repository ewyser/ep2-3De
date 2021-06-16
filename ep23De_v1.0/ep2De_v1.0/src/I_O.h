// header I_O.h to import data and set initial conditions
// GEOMETRY
    load(param,9,1,"param.dat",DAT);
    int nmp   = (int) param_h[0];
    int nn    = (int) param_h[1];
    int no    = (int) param_h[2];
    DAT dx    = (DAT) param_h[3];
    DAT dy    = (DAT) param_h[4];
    DAT xnmin = (DAT) param_h[5];
    DAT ynmin = (DAT) param_h[6];
    int nnx   = (int) param_h[7];
    int nny   = (int) param_h[8];
    int nex   = (int) nnx-1;
    int ney   = (int) nny-1;
    int nel   = nex*ney;
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
    load(xn  ,no*2      ,1,"xn.dat"  ,DAT);
    load(xp  ,2*nmp     ,1,"xp.dat"  ,DAT);
    load(vol ,nmp       ,1,"vol.dat" ,DAT);
    load(lp  ,nmp*2     ,1,"lp.dat"  ,DAT);
    load(Del ,16        ,1,"Del.dat" ,DAT);
    load(iDel,16        ,1,"iDel.dat",DAT);
    load(cohp,nmp       ,1,"cohp.dat",DAT);
    load(phip,nmp       ,1,"phip.dat",DAT);
// INITIALIZE
// elements
    zeros(pel ,2*nel    ,1           ,DAT);
// nodes
    zeros(mn  ,no       ,1           ,DAT);
    zeros(pn  ,2*no     ,1           ,DAT);
    zeros(fen ,2*no     ,1           ,DAT);
    zeros(fin ,2*no     ,1           ,DAT);
    zeros(fn  ,2*no     ,1           ,DAT);
    zeros(an  ,2*no     ,1           ,DAT);
    zeros(vn  ,2*no     ,1           ,DAT);
    zeros(un  ,2*no     ,1           ,DAT);
// material points
    zeros(p2e ,nmp      ,1           ,int);
    zeros(p2n ,nmp*nn   ,1           ,int);
    zeros(sig ,nmp*4    ,1           ,DAT);
    zeros(epII,nmp      ,1           ,DAT);
    zeros(N   ,nmp*nn   ,1           ,DAT);
    zeros(dNx ,nmp*nn   ,1           ,DAT);
    zeros(dNy ,nmp*nn   ,1           ,DAT);
    zeros(vp  ,nmp*2    ,1           ,DAT);
    zeros(up  ,nmp*2    ,1           ,DAT);
    zeros(dF  ,nmp*4    ,1           ,DAT);
    zeros(eps ,nmp*3    ,1           ,DAT);
    zeros(ome ,nmp      ,1           ,DAT);
    zeros(dev ,nmp*4    ,1           ,DAT);    