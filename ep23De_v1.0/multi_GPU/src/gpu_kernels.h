// header gpu_kernels.h : set of all functions call by the device
//-----------------------------------------------------------------------//
//SOLVER-----------------------------------------------------------------//
//-----------------------------------------------------------------------//
__global__ void initD(DAT* mn, DAT* pn, DAT* fen, DAT* fin, DAT* pel, int no, int nel){
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    if(ix<no){
        mn [ix      ] = (DAT)0.0; 
        pn [ix+0*no ] = (DAT)0.0;
        pn [ix+1*no ] = (DAT)0.0;
        pn [ix+2*no ] = (DAT)0.0;
        fen[ix+0*no ] = (DAT)0.0;
        fen[ix+1*no ] = (DAT)0.0;
        fen[ix+2*no ] = (DAT)0.0;
        fin[ix+0*no ] = (DAT)0.0;
        fin[ix+1*no ] = (DAT)0.0;
        fin[ix+2*no ] = (DAT)0.0;
    }
    if(ix<nel){        
        pel[ix+0*nel] = (DAT)0.0;
        pel[ix+1*nel] = (DAT)0.0;
    }
}// WRITE: 8*no+2*nel || READ: - ||  TOTAL IO: 8*no+2*nel
//-----------------------------------------------------------------------//
__global__ void basisD(DAT* N, DAT* dNx, DAT* dNy, DAT* dNz, int* p2n, int* p2e, int* e2N, DAT* xp, DAT* lp, DAT* xn, DAT dx, DAT dy, DAT dz, int nmp, int nn, int nex, int ney, int nez, int nel, int no, DAT xnmin, DAT ynmin, DAT znmin){
#define iD  ix+iy*nmp
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    int iy = blockIdx.y*blockDim.y + threadIdx.y; // thread ID, dimension y
    DAT Nxt,Nyt,Nzt,dNxt,dNyt,dNzt,xi,et,ze;
    if ((ix<nmp)&&(iy<nn)){
        p2e[ix] = (int)((floor((xp[ix+2*nmp]-znmin)*(1.0/dz)))+(nez)*floor((xp[ix+0*nmp]-xnmin)*(1.0/dx))+nez*nex*(floor((xp[ix+1*nmp]-ynmin)*(1.0/dy))));
        // Find connectivity
        p2n[iD] = e2N[p2e[ix]+iy*nel];
        // Calculate xi,eta and zeta
        xi = xp[ix+0*nmp] - xn[p2n[iD]+0*no];
        et = xp[ix+1*nmp] - xn[p2n[iD]+1*no];
        ze = xp[ix+2*nmp] - xn[p2n[iD]+2*no];
        // Initialize basis functions and derivative
        Nxt = Nyt = Nzt = dNxt = dNyt = dNzt = (DAT)0.0;
        // Compute basis functions and derivatives
        if( fabs(xi)< (   lp[ix+0*nmp])                            ){
            Nxt =(DAT)1.0-((DAT)4.0*(xi*xi)+(((DAT)2.0*lp[ix+0*nmp])*((DAT)2.0*lp[ix+0*nmp])))*((DAT)1.0/((DAT)8.0*dx*lp[ix+0*nmp]));
            dNxt=(-(DAT)8.0*xi)*((DAT)1.0/((DAT)8.0*dx*lp[ix+0*nmp]));
        }
        if((fabs(xi)>=(   lp[ix+0*nmp])) && (fabs(xi)<(dx-lp[ix+0*nmp]))){
            Nxt=1.0-(fabs(xi)*((DAT)1.0/dx));
            if(xi<(DAT)0.0){dNxt=( (DAT)1.0/dx);}
            if(xi>(DAT)0.0){dNxt=(-(DAT)1.0/dx);}
        }
        if((fabs(xi)>=(dx-lp[ix+0*nmp])) && (fabs(xi)<(dx+lp[ix+0*nmp]))){
            Nxt=((dx+lp[ix+0*nmp]-fabs(xi))*(dx+lp[ix+0*nmp]-fabs(xi)))*((DAT)1.0/((DAT)4.0*dx*lp[ix+0*nmp]));
            if(xi<(DAT)0.0){dNxt=(  dx+lp[ix+0*nmp]-fabs(xi))*((DAT)1.0/((DAT)2.0*dx*lp[ix+0*nmp]));}
            if(xi>(DAT)0.0){dNxt=(-(dx+lp[ix+0*nmp]-fabs(xi))*((DAT)1.0/((DAT)2.0*dx*lp[ix+0*nmp])));}
        }
        
        
        
        
        if( fabs(et)< (   lp[ix+1*nmp])                            ){
            Nyt =(DAT)1.0-((DAT)4.0*(et*et)+(((DAT)2.0*lp[ix+1*nmp])*((DAT)2.0*lp[ix+1*nmp])))*((DAT)1.0/((DAT)8.0*dy*lp[ix+1*nmp]));
            dNyt=(-(DAT)8.0*et)*((DAT)1.0/((DAT)8.0*dy*lp[ix+1*nmp]));
        }
        if((fabs(et)>=(   lp[ix+1*nmp])) && (fabs(et)<(dy-lp[ix+1*nmp]))){
            Nyt=(DAT)1.0-(fabs(et)*((DAT)1.0/dy));
            if(et<(DAT)0.0){dNyt=( (DAT)1.0/dy);}
            if(et>(DAT)0.0){dNyt=(-(DAT)1.0/dy);}
        }
        if((fabs(et)>=(dy-lp[ix+1*nmp])) && (fabs(et)<(dy+lp[ix+1*nmp]))){
            Nyt=((dy+lp[ix+1*nmp]-fabs(et))*(dy+lp[ix+1*nmp]-fabs(et)))*((DAT)1.0/((DAT)4.0*dy*lp[ix+1*nmp]));
            if(et<(DAT)0.0){dNyt=(  dy+lp[ix+1*nmp]-fabs(et))*((DAT)1.0/((DAT)2.0*dy*lp[ix+1*nmp]));}
            if(et>(DAT)0.0){dNyt=(-(dy+lp[ix+1*nmp]-fabs(et))*((DAT)1.0/((DAT)2.0*dy*lp[ix+1*nmp])));}
        }
        
        
        
        
        if( fabs(ze)< (   lp[ix+2*nmp])                            ){
            Nzt =(DAT)1.0-((DAT)4.0*(ze*ze)+(((DAT)2.0*lp[ix+2*nmp])*((DAT)2.0*lp[ix+2*nmp])))*((DAT)1.0/((DAT)8.0*dz*lp[ix+2*nmp]));
            dNzt=(-(DAT)8.0*ze)*((DAT)1.0/((DAT)8.0*dz*lp[ix+2*nmp]));
        }
        if((fabs(ze)>=(   lp[ix+2*nmp])) && (fabs(ze)<(dz-lp[ix+2*nmp]))){
            Nzt=(DAT)1.0-(fabs(ze)*((DAT)1.0/dz));
            if(ze<(DAT)0.0){dNzt=( (DAT)1.0/dz);}
            if(ze>(DAT)0.0){dNzt=(-(DAT)1.0/dz);}
        }
        if((fabs(ze)>=(dz-lp[ix+2*nmp])) && (fabs(ze)<(dz+lp[ix+2*nmp]))){
            Nzt=((dz+lp[ix+2*nmp]-fabs(ze))*(dz+lp[ix+2*nmp]-fabs(ze)))*((DAT)1.0/((DAT)4.0*dz*lp[ix+2*nmp]));
            if(ze<(DAT)0.0){dNzt=(  dz+lp[ix+2*nmp]-fabs(ze))*((DAT)1.0/((DAT)2.0*dz*lp[ix+2*nmp]));}
            if(ze>(DAT)0.0){dNzt=(-(dz+lp[ix+2*nmp]-fabs(ze))*((DAT)1.0/((DAT)2.0*dz*lp[ix+2*nmp])));}
        }
        // Convolution of basis functions and derivatives
        N  [iD] =  Nxt* Nyt* Nzt;
        dNx[iD] = dNxt* Nyt* Nzt;
        dNy[iD] =  Nxt*dNyt* Nzt;
        dNz[iD] =  Nxt* Nyt*dNzt;
    }
#undef iD
}// WRITE: 4*nn*nmp || READ: 6*nmp+3*nn*nmp  ||  TOTAL IO: nmp*(7*nn+6)
//-----------------------------------------------------------------------//
__global__ void accumD(DAT* mn, DAT* pn, DAT* fen, DAT* fin, DAT* N, DAT* dNx, DAT* dNy, DAT* dNz, DAT* sig, DAT* vp, DAT* mp, DAT* vol, int* p2n, DAT g, int nmp, int nn, int no){
#define iD  p2n[ix+iy*nmp]
#define iDx iD+0*no
#define iDy iD+1*no
#define iDz iD+2*no
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    int iy = blockIdx.y*blockDim.y + threadIdx.y; // thread ID, dimension y
    int id;
    DAT cache;
    if ((ix<nmp)&&(iy<nn)){
        id = ix+iy*nmp                                                                                  ;
        cache = N[id]*mp[ix];
        atomicAdd(&mn [iD ],cache                                                               ); // map material point mass on nodes
        atomicAdd(&pn [iDx],cache*vp[ix+0*nmp]                                                  );
        atomicAdd(&pn [iDy],cache*vp[ix+1*nmp]                                                  );
        atomicAdd(&pn [iDz],cache*vp[ix+2*nmp]                                                  );
        atomicAdd(&fen[iDz],cache*(-g)                                                          );
        atomicAdd(&fin[iDx],vol[ix]*(dNx[id]*sig[ix+0*nmp]+dNy[id]*sig[ix+3*nmp]+dNz[id]*sig[ix+5*nmp]));// xx xy xz 0 3 5
        atomicAdd(&fin[iDy],vol[ix]*(dNx[id]*sig[ix+3*nmp]+dNy[id]*sig[ix+1*nmp]+dNz[id]*sig[ix+4*nmp]));// yx yy yz 3 1 4
        atomicAdd(&fin[iDz],vol[ix]*(dNx[id]*sig[ix+5*nmp]+dNy[id]*sig[ix+4*nmp]+dNz[id]*sig[ix+2*nmp]));// zx zy zz 5 4 2
    }
#undef iD
#undef iDx
#undef iDy
#undef iDz
}// WRITE: 8*nmp*nn || READ: 5*nmp+4*nmp*nn+6*nmp || TOTAL IO: nmp*(11+12*nn)
//-----------------------------------------------------------------------//
__global__ void solveD(DAT* fen, DAT* fin, DAT* mn, DAT* an, DAT* pn, DAT* vn, DAT* un, DAT dt, int no, int* bc){       
#define iDx ix+0*no  
#define iDy ix+1*no 
#define iDz ix+2*no 
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    DAT dmp,m,fx,fy,fz,px,py,pz,vx,vy,vz;
    // solve momentum equation on the background mesh
    if(ix<no){
        an[iDx]=(DAT)0.0;
        vn[iDx]=(DAT)0.0;
        un[iDx]=(DAT)0.0;
        an[iDy]=(DAT)0.0;
        vn[iDy]=(DAT)0.0;
        un[iDy]=(DAT)0.0;
        an[iDz]=(DAT)0.0;
        vn[iDz]=(DAT)0.0;
        un[iDz]=(DAT)0.0;
        if(mn[ix]>(DAT)0.0){
            m   = (DAT)1.0/(DAT)mn[ix]            ;
            fx  = fen[iDx]-fin[iDx]               ;
            fy  = fen[iDy]-fin[iDy]               ;
            fz  = fen[iDz]-fin[iDz]               ;
            px  = pn[iDx]                         ;
            py  = pn[iDy]                         ;
            pz  = pn[iDz]                         ;
            vx  = px*m                            ;
            vy  = py*m                            ;
            vz  = pz*m                            ;
            dmp = sqrt(fx*fx+fy*fy+fz*fz)         ; 
            fx  = fx-D*dmp*((DAT)vx/(DAT)fabs(vx));
            fy  = fy-D*dmp*((DAT)vy/(DAT)fabs(vy));
            fz  = fz-D*dmp*((DAT)vz/(DAT)fabs(vz));
            if(fabs(vx)<(DAT)1E-6){
                fx = fen[iDx]-fin[iDx];
            }
            if(fabs(vy)<(DAT)1E-6){
                fy = fen[iDy]-fin[iDy];
            }
            if(fabs(vz)<(DAT)1E-6){
                fz = fen[iDz]-fin[iDz];
            }
            an[iDx] = (fx        )*m*(DAT)bc[iDx]    ;
            vn[iDx] = (px + dt*fx)*m*(DAT)bc[iDx]    ;            
            an[iDy] = (fy        )*m*(DAT)bc[iDy]    ;
            vn[iDy] = (py + dt*fy)*m*(DAT)bc[iDy]    ;
            an[iDz] = (fz        )*m*(DAT)bc[iDz]    ;
            vn[iDz] = (pz + dt*fz)*m*(DAT)bc[iDz]    ;
        }
    }
#undef iDx
#undef iDy
#undef iDz
}// WRITE: 9*no || READ: 9*no || TOTAL IO: 18*no
//-----------------------------------------------------------------------//
__global__ void projectD(DAT* an, DAT* vn, DAT* N, DAT* vp, DAT* xp, DAT* up, int* p2n, DAT dt, int nmp, int nn, int no){
#define iD  ix+k*nmp
#define iDx p2n[iD]+0*no
#define iDy p2n[iD]+1*no
#define iDz p2n[iD]+2*no
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    DAT dax,day,daz,dvx,dvy,dvz,vxp,vyp,vzp;
    // initialize cache
    dax=day=daz=dvx=dvy=dvz=(DAT)0.0;
    if ((ix<nmp)){
        // parallel reduction
        for(int k=0;k<nn;k++){
            dax+=N[iD]*an[iDx];
            day+=N[iD]*an[iDy];
            daz+=N[iD]*an[iDz];
            dvx+=N[iD]*vn[iDx];
            dvy+=N[iD]*vn[iDy];
            dvz+=N[iD]*vn[iDz];
        }
        // update velocities
        vp[ix+0*nmp] += dt*dax;
        vp[ix+1*nmp] += dt*day;
        vp[ix+2*nmp] += dt*daz;
        // update coordinates
        xp[ix+0*nmp] += dt*dvx;
        xp[ix+1*nmp] += dt*dvy;
        xp[ix+2*nmp] += dt*dvz;
        // update incremental displacements
        up[ix+0*nmp] += dt*dvx;
        up[ix+1*nmp] += dt*dvy;
        up[ix+2*nmp] += dt*dvz;
    }
#undef iD
#undef iDx
#undef iDy
#undef iDz
}// WRITE: 9*nmp || READ: 3*nn*nmp || TOTAL IO: 3*nmp*(3+nn)
//-----------------------------------------------------------------------//
__global__ void DMD(DAT* mn, DAT* un, DAT* N, DAT* vp, DAT* mp, int* p2n, DAT dt, int nmp, int nn, int* bc, int no){
#define iD  p2n[ix+iy*nmp]
#define iDx iD+0*no
#define iDy iD+1*no
#define iDz iD+2*no
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    int iy = blockIdx.y*blockDim.y + threadIdx.y; // thread ID, dimension y
    DAT prod;
    if ((ix<nmp)&&(iy<nn)&&(mn[iD]>(DAT)0.0)){
        prod = dt*N[ix+iy*nmp]*mp[ix]*((DAT)1.0/(DAT)mn[iD]);
        atomicAdd(&un[iDx],prod*vp[ix+0*nmp]*(DAT)bc[iDx]);
        atomicAdd(&un[iDy],prod*vp[ix+1*nmp]*(DAT)bc[iDy]);
        atomicAdd(&un[iDz],prod*vp[ix+2*nmp]*(DAT)bc[iDz]);
    }
#undef iD    
#undef iDx
#undef iDy
#undef iDz
}// WRITE: 3*nmp*nn || READ: 1*nmp*nn+3*nmp || TOTAL IO: nmp*(4*nn+3)
                            __global__ void DMD1(DAT* un, DAT* N, DAT* vp, DAT* mp, int* p2n, int nmp, int nn, int no){
                            #define iD  p2n[ix+iy*nmp]
                            #define iDx iD+0*no
                            #define iDy iD+1*no
                            #define iDz iD+2*no
                                int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
                                int iy = blockIdx.y*blockDim.y + threadIdx.y; // thread ID, dimension y
                                DAT prod;
                                if ((ix<nmp)&&(iy<nn)){
                                    prod = N[ix+iy*nmp]*mp[ix];
                                    atomicAdd(&un[iDx],prod*vp[ix+0*nmp]);
                                    atomicAdd(&un[iDy],prod*vp[ix+1*nmp]);
                                    atomicAdd(&un[iDz],prod*vp[ix+2*nmp]);
                                }
                            #undef iD    
                            #undef iDx
                            #undef iDy
                            #undef iDz
                            }// WRITE: 3*nmp*nn || READ: 1*nmp*nn+3*nmp || TOTAL IO: nmp*(4*nn+3)
                            __global__ void DMD2(DAT* mn, DAT* un, DAT dt, int* bc, int no){
                            #define iDx ix+0*no
                            #define iDy ix+1*no
                            #define iDz ix+2*no
                                int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
                                DAT m=(DAT)0.0;
                                if(ix<no){
                                    if(mn[ix]>(DAT)0.0){
                                        m   = (DAT)1.0/(DAT)mn[ix]            ;
                                        un[iDx] = (dt*un[iDx])*m*(DAT)bc[iDx]    ;            
                                        un[iDy] = (dt*un[iDy])*m*(DAT)bc[iDy]    ;
                                        un[iDz] = (dt*un[iDz])*m*(DAT)bc[iDz]    ;
                                    }
                                    else{
                                        un[iDx]=(DAT)0.0;
                                        un[iDy]=(DAT)0.0;
                                        un[iDz]=(DAT)0.0;
                                    }
                                }
                            #undef iDx
                            #undef iDy
                            #undef iDz
                            }// WRITE: 3*nmp*nn || READ: 1*nmp*nn+3*nmp || TOTAL IO: nmp*(4*nn+3)
//-----------------------------------------------------------------------//
__global__ void getdFD(DAT* un, DAT* dF, DAT* dNx, DAT* dNy, DAT* dNz, int* p2n, int nmp, int nn, int no){
#define iD  ip+k*nmp
#define iDx p2n[iD]+0*no
#define iDy p2n[iD]+1*no
#define iDz p2n[iD]+2*no
    int ip = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    DAT dFxx,dFxy,dFxz,dFyx,dFyy,dFyz,dFzx,dFzy,dFzz;
    // initialize cache
    dFxx=dFxy=dFxz=dFyx=dFyy=dFyz=dFzx=dFzy=dFzz=(DAT)0.0;
    if ((ip<nmp)){
        // parallel reduction
        for(int k=0;k<nn;k++){
            dFxx+=dNx[iD]*un[iDx];
            dFxy+=dNy[iD]*un[iDx];
            dFxz+=dNz[iD]*un[iDx];
            dFyx+=dNx[iD]*un[iDy];
            dFyy+=dNy[iD]*un[iDy];
            dFyz+=dNz[iD]*un[iDy];
            dFzx+=dNx[iD]*un[iDz];
            dFzy+=dNy[iD]*un[iDz];
            dFzz+=dNz[iD]*un[iDz];
        }
        // incremental deformation gradient
        dF[ip+0*nmp]=(DAT)1.0+dFxx;
        dF[ip+1*nmp]=(DAT)0.0+dFxy;
        dF[ip+2*nmp]=(DAT)0.0+dFxz;
        dF[ip+3*nmp]=(DAT)0.0+dFyx;
        dF[ip+4*nmp]=(DAT)1.0+dFyy;
        dF[ip+5*nmp]=(DAT)0.0+dFyz;
        dF[ip+6*nmp]=(DAT)0.0+dFzx;
        dF[ip+7*nmp]=(DAT)0.0+dFzy;
        dF[ip+8*nmp]=(DAT)1.0+dFzz;
    }
#undef iD
#undef iDx
#undef iDy
#undef iDz
}// WRITE: nmp*9 || READ: 6*nmp*nn || TOTAL IO: 3*nmp*(3+6*nn)
//-----------------------------------------------------------------------//
__global__ void elastD(DAT* dF, DAT* eps, DAT* sig, DAT* ome, DAT* lp, DAT* vol, int nmp, DAT* Del, DAT dt){
    int p = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    DAT DT,dexx,deyy,dezz,dexy,dexz,deyz,doxy,doxz,doyz,sxx0,syy0,szz0,sxy0,sxz0,syz0,J;
    if (p<nmp){
        // store previous stresses
        sxx0 = sig[p+0*nmp];
        syy0 = sig[p+1*nmp];
        szz0 = sig[p+2*nmp];
        sxy0 = sig[p+3*nmp];
        syz0 = sig[p+4*nmp];
        sxz0 = sig[p+5*nmp];
        // compute strains 
        DT   = (DAT)1.0/(DAT)dt;
        dexx = (dF[p+0*nmp]-(DAT)1.0   )*DT;
        deyy = (dF[p+4*nmp]-(DAT)1.0   )*DT;
        dezz = (dF[p+8*nmp]-(DAT)1.0   )*DT;
        dexy = (dF[p+1*nmp]+dF[p+3*nmp])*DT;
        deyz = (dF[p+7*nmp]+dF[p+5*nmp])*DT;
        dexz = (dF[p+6*nmp]+dF[p+2*nmp])*DT;
        // compute spin rate
        doxy = (DAT)0.5*(dF[p+3*nmp]-dF[p+1*nmp])*DT;  
        doyz = (DAT)0.5*(dF[p+5*nmp]-dF[p+7*nmp])*DT; 
        doxz = (DAT)0.5*(dF[p+2*nmp]-dF[p+6*nmp])*DT; 
        // update objective stress        
        sig[p+0*nmp] +=  (DAT)2*dt*(sxy0*doxy+sxz0*doxz);
        sig[p+1*nmp] += -(DAT)2*dt*(sxy0*doxy-syz0*doyz);
        sig[p+2*nmp] += -(DAT)2*dt*(sxz0*doxz+syz0*doyz);
        sig[p+3*nmp] +=         dt*(doxy*(syy0-sxx0)+syz0*doxz+sxz0*doyz);
        sig[p+4*nmp] +=         dt*(doyz*(szz0-syy0)-sxy0*doxz-sxz0*doxy);
        sig[p+5*nmp] +=         dt*(doxz*(szz0-sxx0)+syz0*doxy-sxy0*doyz);
        // incremental strain
        sig[p+0*nmp] += dt*(Del[0]*dexx+Del[6 ]*deyy+Del[12]*dezz+Del[18]*dexy+Del[24]*deyz+Del[30]*dexz);
        sig[p+1*nmp] += dt*(Del[1]*dexx+Del[7 ]*deyy+Del[13]*dezz+Del[19]*dexy+Del[25]*deyz+Del[31]*dexz);
        sig[p+2*nmp] += dt*(Del[2]*dexx+Del[8 ]*deyy+Del[14]*dezz+Del[20]*dexy+Del[26]*deyz+Del[32]*dexz);
        sig[p+3*nmp] += dt*(Del[3]*dexx+Del[9 ]*deyy+Del[15]*dezz+Del[21]*dexy+Del[27]*deyz+Del[33]*dexz);
        sig[p+4*nmp] += dt*(Del[4]*dexx+Del[10]*deyy+Del[16]*dezz+Del[22]*dexy+Del[28]*deyz+Del[34]*dexz);
        sig[p+5*nmp] += dt*(Del[5]*dexx+Del[11]*deyy+Del[17]*dezz+Del[23]*dexy+Del[29]*deyz+Del[35]*dexz);
        // update material point volume and domain lengths
        J           = (DAT)1.0+dt*(dexx+deyy+dezz);
        vol[p]      = J*vol[p];
        if(sizeof(DAT)==(int)8){
            J           = pow(J,(DAT)0.3333);
        }
        else if(sizeof(DAT)==(int)4){
            J           = powf(J,(DAT)0.3333);
        }
        lp[p+0*nmp] = J*lp[p+0*nmp];
        lp[p+1*nmp] = J*lp[p+1*nmp];
        lp[p+2*nmp] = J*lp[p+2*nmp];
    } 
}// WRITE: 9*nmp || READ: 12*nmp || TOTAL IO: 21*nmp
//-----------------------------------------------------------------------//
__global__ void volLockD1(DAT* pel,DAT* sig, DAT* dev, DAT* vol, int* p2e, int nmp, int nel){
#define iD  p2e[ix]
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    DAT pr,sxx,syy,szz,sxy,syz,sxz;
    // accumulate material point pressure on element
    if (ix<nmp){
        sxx          =sig[ix+0*nmp]                     ;
        syy          =sig[ix+1*nmp]                     ;
        szz          =sig[ix+2*nmp]                     ;
        sxy          =sig[ix+3*nmp]                     ;
        syz          =sig[ix+4*nmp]                     ;
        sxz          =sig[ix+5*nmp]                     ;
        pr           =-(sxx+syy+szz)*((DAT)1.0/(DAT)3.0);
        dev[ix+0*nmp]=sxx+pr                            ;
        dev[ix+1*nmp]=syy+pr                            ;
        dev[ix+2*nmp]=szz+pr                            ;
        dev[ix+3*nmp]=sxy                               ;
        dev[ix+4*nmp]=syz                               ;
        dev[ix+5*nmp]=sxz                               ;
        atomicAdd(&pel[iD+0*nel],pr*vol[ix])            ; 
        atomicAdd(&pel[iD+1*nel],   vol[ix])            ;
    }
}// WRITE: 8*nmp || READ: 6*nmp || TOTAL IO: 14*nmp
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
__global__ void volLockD2(DAT* pel,DAT* sig, DAT* dev, int* p2e, int nmp, int nel){
#define iD  p2e[ix]
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    DAT pr;
    // assign element pressure to material point
    if (ix<nmp){
        pr           =pel[iD+0*nel]*((DAT)1.0/(DAT)pel[iD+1*nel]);
        sig[ix+0*nmp]=dev[ix+0*nmp]-pr                           ;
        sig[ix+1*nmp]=dev[ix+1*nmp]-pr                           ;
        sig[ix+2*nmp]=dev[ix+2*nmp]-pr                           ;
        sig[ix+3*nmp]=dev[ix+3*nmp]                              ;
        sig[ix+4*nmp]=dev[ix+4*nmp]                              ;
        sig[ix+5*nmp]=dev[ix+5*nmp]                              ;
    }
#undef iD
}// WRITE: 6*nmp || READ: 7*nmp || TOTAL IO: 13*nmp
//-----------------------------------------------------------------------//
//UTILITIES--------------------------------------------------------------//
DAT checkGpuMem(int me){
    double free_m,total_m,used_m;
    size_t free_t,total_t;
    cudaMemGetInfo(&free_t,&total_t);
    free_m = free_t          /(double)(1024.0*1024.0*1024.0);
    total_m= total_t         /(double)(1024.0*1024.0*1024.0);
    used_m = (total_t-free_t)/(double)(1024.0*1024.0*1024.0);
    if(me==0){
    printf("\n-----------------------------------");
    if(sizeof(DAT)==8){
        printf("\n-----------------------------------");
        printf("\nGPU memory: double arithmetic");
    }
    if(sizeof(DAT)==4){
        printf("\nGPU memory: float arithmetic");
    }
    printf("\n-----------------------------------");
    printf("\n  total: %.2f GB --> 100 %%" ,total_m                    );
    printf("\n  free : %.2f GB --> %.2f %%",free_m,free_m/total_m*100.0);
    printf("\n  used : %.2f GB --> %.2f %%",used_m,used_m/total_m*100.0);
    }
    return((DAT)((total_m-free_m)/total_m));   
}
//-----------------------------------------------------------------------//
void checkCudaErrors(){
    cudaError_t err = cudaGetLastError();
    if ( err != cudaSuccess )
    {
        printf("\n /!\\ GPU error: %s\n", cudaGetErrorString(err));
        exit(-1);
    }
}
//MISC-------------------------------------------------------------------//
//-----------------------------------------------------------------------//
__device__ DAT getN(DAT xi,DAT lp,DAT h){
    DAT N;
    if     (xi>(-h-lp) && xi<=(-h+lp)){N = (h+lp+xi)*(h+lp+xi)*((DAT)1.0/((DAT)4.0*h*lp));}
    else if(xi>(-h+lp) && xi<=(  -lp)){N = 1.0+xi             *((DAT)1.0/((DAT)h       ));}
    else if(xi>(  -lp) && xi<=(   lp)){N = 1.0-((xi*xi)+lp*lp)*((DAT)1.0/((DAT)2.0*h*lp));}
    else if(xi>(   lp) && xi<=( h-lp)){N = 1.0-xi             *((DAT)1.0/((DAT)h       ));}
    else if(xi>( h-lp) && xi<=( h+lp)){N = (h+lp-xi)*(h+lp-xi)*((DAT)1.0/((DAT)4.0*h*lp));}
    else                              {N = 0.0                                           ;}
    return(N);
}
__device__ DAT getdN(DAT xi,DAT lp,DAT h){
    DAT dN;
    if     (xi>(-h-lp) && xi<=(-h+lp)){dN = (h+lp+xi)*((DAT)1.0/((DAT)2.0*h*lp));}
    else if(xi>(-h+lp) && xi<=(  -lp)){dN =           ((DAT)1.0/((DAT)h       ));}
    else if(xi>(  -lp) && xi<=(   lp)){dN =-xi       *((DAT)1.0/((DAT)1.0*h*lp));}
    else if(xi>(   lp) && xi<=( h-lp)){dN =-1.0      *((DAT)1.0/((DAT)h       ));}
    else if(xi>( h-lp) && xi<=( h+lp)){dN =-(h+lp-xi)*((DAT)1.0/((DAT)2.0*h*lp));}
    else                              {dN = 0.0                                 ;}
    return(dN);
}
__global__ void basisD_TEST(DAT* N, DAT* dNx, DAT* dNy, DAT* dNz, int* p2n, int* p2e, int* e2N, DAT* xp, DAT* lp, DAT* xn, DAT dx, DAT dy, DAT dz, int nmp, int nn, int nex, int ney, int nez, int nel, int no, DAT xnmin, DAT ynmin, DAT znmin){
#define iD  ix+iy*nmp
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    int iy = blockIdx.y*blockDim.y + threadIdx.y; // thread ID, dimension y
    DAT Nxt,Nyt,Nzt,dNxt,dNyt,dNzt,edx,edy,edz;
    if ((ix<nmp)&&(iy<nn)){
        
        edx  = (xp[ix+0*nmp]-xnmin)                                            ;
        edy  = (xp[ix+1*nmp]-ynmin)                                            ;
        edz  = (xp[ix+2*nmp]-znmin)                                            ;
        
        p2e[ix] = (int)((floor((edz)*(1.0/dz)))+(nez)*floor(edx*(1.0/dx))+nez*nex*(floor(edy*(1.0/dy))));
        // Find connectivity
        p2n[iD] = e2N[p2e[ix]+iy*nel];
        // Compute basis functions and derivatives
        Nxt  = getN (xp[ix+0*nmp]-xn[p2n[iD]+0*no],lp[ix+0*nmp],dx);
        dNxt = getdN(xp[ix+0*nmp]-xn[p2n[iD]+0*no],lp[ix+0*nmp],dx);
        Nyt  = getN (xp[ix+1*nmp]-xn[p2n[iD]+1*no],lp[ix+1*nmp],dy);
        dNyt = getdN(xp[ix+1*nmp]-xn[p2n[iD]+1*no],lp[ix+1*nmp],dy); 
        Nzt  = getN (xp[ix+2*nmp]-xn[p2n[iD]+2*no],lp[ix+2*nmp],dz);
        dNzt = getdN(xp[ix+2*nmp]-xn[p2n[iD]+2*no],lp[ix+2*nmp],dz); 
        // Convolution of basis functions and derivatives
        N  [iD] =  Nxt* Nyt * Nzt;
        dNx[iD] = dNxt* Nyt * Nzt;
        dNy[iD] =  Nxt*dNyt * Nzt;
        dNz[iD] =  Nxt* Nyt *dNzt;
    }
#undef iD
}
//-----------------------------------------------------------------------//
__global__ void DPPlastD(DAT* sig, DAT* cohp, DAT* phip, DAT* epII, DAT Hp, DAT cohr, DAT Kc, DAT Gc, DAT psi, int nmp){
    int p = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    DAT c,Pr,tensile,sxx,syy,szz,sxy,syz,sxz,devxx,devyy,devzz,devxy,devyz,devxz,J2,tau,eta,etaB,xi,sigm,fs,ft,tauP,alpP,h,dlam,tauN,PrN,dep;
    tensile = (DAT)0.0;
    if(p<nmp){
        c  = cohp[p]+(Hp*epII[p]);
        if(c<cohr){c = cohr;}
        sxx   = sig[p+0*nmp] ;
        syy   = sig[p+1*nmp] ;
        szz   = sig[p+2*nmp] ;
        sxy   = sig[p+3*nmp] ;
        syz   = sig[p+4*nmp] ;
        sxz   = sig[p+5*nmp] ;
        Pr    = (sxx+syy+szz)*((DAT)1.0/(DAT)3.0);
        devxx = sxx - Pr     ; 
        devyy = syy - Pr     ;
        devzz = szz - Pr     ;
        devxy = sxy          ;
        devyz = syz          ;
        devxz = sxz          ;
        tau   = sqrt((DAT)0.5*(devxx*devxx+devyy*devyy+devzz*devzz)+devxy*devxy+devyz*devyz+devxz*devxz);
        eta   = ((DAT)6.0*  sin(phip[p]))*(DAT)1.0/((DAT)sqrt((DAT)3.0)*((DAT)3.0+sin(phip[p])));
        etaB  = ((DAT)6.0*  sin(psi    ))*(DAT)1.0/((DAT)sqrt((DAT)3.0)*((DAT)3.0+sin(psi    )));
        xi    = ((DAT)6.0*c*cos(phip[p]))*(DAT)1.0/((DAT)sqrt((DAT)3.0)*((DAT)3.0+sin(phip[p])));        
    
        
        //sigm  = min(tensile,xi*(DAT)1.0/((DAT)eta));
        sigm  = xi*(DAT)1.0/((DAT)eta);
        fs    = tau+eta*Pr-xi;
        ft    = Pr-sigm         ;
        tauP  = xi-eta*sigm  ;
        alpP  = sqrt((DAT)1.0+eta*eta)-eta;
        h     = tau-tauP-alpP*(Pr-sigm);
        if((fs>(DAT)0.0 && Pr<sigm)||(h>(DAT)0.0 && Pr>=sigm)){
            dlam          = fs*(DAT)1.0/((DAT)Gc+Kc*eta*etaB)           ;
            PrN           = Pr-Kc*etaB*dlam                             ;
            tauN          = xi-eta*PrN                                  ;
            sig [p+0*nmp] = devxx*((DAT)tauN/((DAT)tau))+PrN            ;
            sig [p+1*nmp] = devyy*((DAT)tauN/((DAT)tau))+PrN            ;
            sig [p+2*nmp] = devzz*((DAT)tauN/((DAT)tau))+PrN            ;
            sig [p+3*nmp] = devxy*((DAT)tauN/((DAT)tau))                ;
            sig [p+4*nmp] = devyz*((DAT)tauN/((DAT)tau))                ;
            sig [p+5*nmp] = devxz*((DAT)tauN/((DAT)tau))                ;
            dep           = (dlam*sqrt((DAT)0.33333333+(DAT)0.22222222*etaB*etaB));
            epII[p      ]+= dep                                         ;
        }
        if((h<=0.0)&&(Pr>=sigm)){
            dlam          = (Pr-sigm)*((DAT)1.0/((DAT)Kc))        ;
            sig[p+0*nmp] += (sigm-Pr)                             ;
            sig[p+1*nmp] += (sigm-Pr)                             ;
            sig[p+2*nmp] += (sigm-Pr)                             ;
            sig[p+3*nmp] += (DAT)0.0                              ;
            sig[p+4*nmp] += (DAT)0.0                              ;
            sig[p+5*nmp] += (DAT)0.0                              ;
            dep           = (sqrt(2.0)*dlam*((DAT)1.0/((DAT)3.0)));
            epII[p      ]+= dep                                   ;
        }  
    }
}
//-----------------------------------------------------------------------//