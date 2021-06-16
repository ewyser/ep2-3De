// header gpu_kernels.h : set of all functions call by the device
//-----------------------------------------------------------------------//
//SOLVER-----------------------------------------------------------------//
//-----------------------------------------------------------------------//
__global__ void initD(DAT* mn, DAT* pn, DAT* fen, DAT* fin, DAT* pel, int no, int nel){
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    if(ix<no){
        mn [ix      ] = 0.0; 
        pn [ix+0*no ] = 0.0;
        pn [ix+1*no ] = 0.0;
        fen[ix+1*no ] = 0.0;
        fin[ix+0*no ] = 0.0;
        fin[ix+1*no ] = 0.0;
    }
    if(ix<nel){        
        pel[ix+0*nel] = 0.0;
        pel[ix+1*nel] = 0.0;
    }
}// WRITE: 6*no+2*nel || READ: - ||  TOTAL IO: 6*no+2*nel
//-----------------------------------------------------------------------//
__global__ void basisD(DAT* N, DAT* dNx, DAT* dNy, int* p2n, int* p2e, int* e2N, DAT* xp, DAT* lp, DAT* xn, DAT dx, DAT dy, int nmp, int nn, int ney, int nel, int no, DAT xnmin, DAT ynmin){
#define iD  ix+iy*nmp
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    int iy = blockIdx.y*blockDim.y + threadIdx.y; // thread ID, dimension y
    DAT Nxt,Nyt,dNxt,dNyt,xi,et;
    if ((ix<nmp)&&(iy<nn)){
        // Find connectivity
        p2e[ix] = (int)(floor((xp[ix+1*nmp]-ynmin)*(1.0/dx))+ney*floor((xp[ix+0*nmp]-xnmin)*(1.0/dy)));
        p2n[iD] = e2N[p2e[ix]+iy*nel];
        // Calculate xi and eta
        xi = xp[ix+0*nmp] - xn[p2n[iD]+0*no];
        et = xp[ix+1*nmp] - xn[p2n[iD]+1*no];
        // Initialize basis functions and derivative
        Nxt = Nyt = dNxt = dNyt = 0.0;
        // Compute basis functions and derivatives
        if( fabs(xi)< (   lp[ix+0*nmp])                            ){
            Nxt =1.0-(4.0*(xi*xi)+((2.0*lp[ix+0*nmp])*(2.0*lp[ix+0*nmp])))*(1.0/(8.0*dx*lp[ix+0*nmp]));
            dNxt=(-8.0*xi)*(1.0/(8.0*dx*lp[ix+0*nmp]));
        }
        if((fabs(xi)>=(   lp[ix+0*nmp])) && (fabs(xi)<(dx-lp[ix+0*nmp]))){
            Nxt=1.0-(fabs(xi)*(1.0/dx));
            if(xi<0.0){dNxt=( 1.0/dx);}
            if(xi>0.0){dNxt=(-1.0/dx);}
        }
        if((fabs(xi)>=(dx-lp[ix+0*nmp])) && (fabs(xi)<(dx+lp[ix+0*nmp]))){
            Nxt=((dx+lp[ix+0*nmp]-fabs(xi))*(dx+lp[ix+0*nmp]-fabs(xi)))*(1.0/(4.0*dx*lp[ix+0*nmp]));
            if(xi<0.0){dNxt=(  dx+lp[ix+0*nmp]-fabs(xi))*(1.0/(2.0*dx*lp[ix+0*nmp]));}
            if(xi>0.0){dNxt=(-(dx+lp[ix+0*nmp]-fabs(xi))*(1.0/(2.0*dx*lp[ix+0*nmp])));}
        }
        
        
        
        
        
        if( fabs(et)< (   lp[ix+1*nmp])                            ){
            Nyt =1.0-(4.0*(et*et)+((2.0*lp[ix+1*nmp])*(2.0*lp[ix+1*nmp])))*(1.0/(8.0*dy*lp[ix+1*nmp]));
            dNyt=(-8.0*et)*(1.0/(8.0*dy*lp[ix+1*nmp]));
        }
        if((fabs(et)>=(   lp[ix+1*nmp])) && (fabs(et)<(dy-lp[ix+1*nmp]))){
            Nyt=1.0-(fabs(et)*(1.0/dy));
            if(et<0.0){dNyt=( 1.0/dy);}
            if(et>0.0){dNyt=(-1.0/dy);}
        }
        if((fabs(et)>=(dy-lp[ix+1*nmp])) && (fabs(et)<(dy+lp[ix+1*nmp]))){
            Nyt=((dy+lp[ix+1*nmp]-fabs(et))*(dy+lp[ix+1*nmp]-fabs(et)))*(1.0/(4.0*dy*lp[ix+1*nmp]));
            if(et<0.0){dNyt=(  dy+lp[ix+1*nmp]-fabs(et))*(1.0/(2.0*dy*lp[ix+1*nmp]));}
            if(et>0.0){dNyt=(-(dy+lp[ix+1*nmp]-fabs(et))*(1.0/(2.0*dy*lp[ix+1*nmp])));}
        }
        // Convolution of basis functions and derivatives
        N  [iD] =  Nxt* Nyt;
        dNx[iD] = dNxt* Nyt;
        dNy[iD] =  Nxt*dNyt;
    }
#undef iD
}// WRITE: 3*nn*nmp || READ: 4*nmp+nn*nmp  ||  TOTAL IO: 4*nmp + 4*nn*nmp
//-----------------------------------------------------------------------//
__global__ void accumD(DAT* mn, DAT* pn, DAT* fen, DAT* fin, DAT* N, DAT* dNx, DAT* dNy, DAT* sig, DAT* vp, DAT* mp, DAT* vol, int* p2n, DAT g, int nmp, int nn, int no){
#define iD  p2n[ix+iy*nmp]
#define iDx iD+0*no
#define iDy iD+1*no
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    int iy = blockIdx.y*blockDim.y + threadIdx.y; // thread ID, dimension y
    int id;
    if ((ix<nmp)&&(iy<nn)){
        id = ix+iy*nmp                                                           ;
        atomicAdd(&mn [iD ],N[id]*mp[ix]                                       ); // map material point mass on nodes
        atomicAdd(&pn [iDx],N[id]*mp[ix]*vp[ix+0*nmp]                           );
        atomicAdd(&pn [iDy],N[id]*mp[ix]*vp[ix+1*nmp]                           );
        atomicAdd(&fen[iDy],N[id]*mp[ix]*(-g)                                  );
        atomicAdd(&fin[iDx],vol[ix]*(dNx[id]*sig[ix+0*nmp]+dNy[id]*sig[ix+3*nmp]));
        atomicAdd(&fin[iDy],vol[ix]*(dNx[id]*sig[ix+3*nmp]+dNy[id]*sig[ix+1*nmp]));
    }
#undef iD
#undef iDx
#undef iDy
}// WRITE: 6*no || READ: 4*nmp+3*nn*nmp+2*nmp*nstr || TOTAL IO: 6*no+nmp*(4+3*nn+2*nstr)
//-----------------------------------------------------------------------//
__global__ void solveD(DAT* fen, DAT* fin, DAT* mn, DAT* an, DAT* pn, DAT* vn, DAT* un, DAT dt, int no, int* bc){       
#define iDx ix+0*no  
#define iDy ix+1*no  
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    DAT dmp,fx,fy,px,py,vx,vy;
    // solve momentum equation on the background mesh
    if(ix<no){
        an[iDx]=0.0;
        vn[iDx]=0.0;
        un[iDx]=0.0;
        an[iDy]=0.0;
        vn[iDy]=0.0;
        un[iDy]=0.0;
        if(mn[ix]>0.0){
            fx  = fen[iDx]-fin[iDx]               ;
            fy  = fen[iDy]-fin[iDy]               ;
            px  = pn[iDx]                         ;
            py  = pn[iDy]                         ;
            vx  = px*(DAT)1.0/(DAT)mn[ix]         ;
            vy  = py*(DAT)1.0/(DAT)mn[ix]         ;
            dmp = sqrt(fx*fx+fy*fy)               ; 
            fx  = fx-D*dmp*((DAT)vx/(DAT)fabs(vx));
            if(fabs(vx)<(DAT)1E-3){
                fx = fen[iDx]-fin[iDx];
            }
            fy  = fy-D*dmp*((DAT)vy/(DAT)fabs(vy));
            if(fabs(vy)<(DAT)1E-3){
                fy = fen[iDy]-fin[iDy];
            }
            dmp     = (DAT)1.0/(DAT)mn[ix]        ;
            an[iDx] = (fx        )*dmp*bc[iDx]    ;
            vn[iDx] = (px + dt*fx)*dmp*bc[iDx]    ;            
            an[iDy] = (fy        )*dmp*bc[iDy]    ;
            vn[iDy] = (py + dt*fy)*dmp*bc[iDy]    ;
        }
    }
#undef iDx
#undef iDy
}// WRITE: 3*2*no || READ: 3*2*no || TOTAL IO: 12*no
//-----------------------------------------------------------------------//
__global__ void projectD(DAT* an, DAT* vn, DAT* N, DAT* vp, DAT* xp, DAT* up, int* p2n, DAT dt, int nmp, int nn, int no){
#define iD  ix+k*nmp
#define iDx p2n[iD]+0*no
#define iDy p2n[iD]+1*no
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    DAT dax,day,dvx,dvy,vxp,vyp;
    // initialize cache
    dax=day=dvx=dvy=0.0;
    if ((ix<nmp)){
        // parallel reduction
        for(int k=0;k<nn;k++){
            dax+=N[iD]*an[iDx];
            day+=N[iD]*an[iDy];
            dvx+=N[iD]*vn[iDx];
            dvy+=N[iD]*vn[iDy];
        }
        // update mp's velocity and coordinates
        if(alpha>0.0){
            vxp           = vp[ix+0*nmp];
            vyp           = vp[ix+1*nmp];
            vp[ix+0*nmp] += (dt*dax)-alpha*(vxp-dvx);
            vp[ix+1*nmp] += (dt*day)-alpha*(vyp-dvy);
            xp[ix+0*nmp] += (dvx-0.5*alpha*(vxp-dvx))*dt+0.5*(dax*dt*dt);
            xp[ix+1*nmp] += (dvy-0.5*alpha*(vyp-dvy))*dt+0.5*(day*dt*dt);
        }
        else{
            vp[ix+0*nmp] += dt*dax;
            vp[ix+1*nmp] += dt*day;
            xp[ix+0*nmp] += dt*dvx;
            xp[ix+1*nmp] += dt*dvy;
            up[ix+0*nmp] += dt*dvx;
            up[ix+1*nmp] += dt*dvy;
        }
    }
#undef iD
#undef iDx
#undef iDy
}// WRITE: 6*nmp || READ: 3*nmp*nn || TOTAL IO: nmp*(6+3*nn)
//-----------------------------------------------------------------------//
__global__ void DMD(DAT* mn, DAT* un, DAT* N, DAT* vp, DAT* mp, int* p2n, DAT dt, int nmp, int nn, int* bc, int no){
#define iD  p2n[ix+iy*nmp]
#define iDx iD+0*no
#define iDy iD+1*no
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    int iy = blockIdx.y*blockDim.y + threadIdx.y; // thread ID, dimension y
    DAT prod;
    if ((ix<nmp)&&(iy<nn)&&(mn[iD]>0.0)){
        prod = dt*N[ix+iy*nmp]*mp[ix]*((DAT)1.0/(DAT)mn[iD]);
        atomicAdd(&un[iDx],prod*vp[ix+0*nmp]*bc[iDx]);
        atomicAdd(&un[iDy],prod*vp[ix+1*nmp]*bc[iDy]);
    }
#undef iD    
#undef iDx
#undef iDy
}// WRITE: 2*no || READ: 1*nmp*nn+2*nmp || TOTAL IO: 2*no+nmp*(nn+2)
//-----------------------------------------------------------------------//
__global__ void getdFD(DAT* un, DAT* dF, DAT* dNx, DAT* dNy, int* p2n, int nmp, int nn, int no){
#define iD  ip+k*nmp
#define iDx p2n[iD]+0*no
#define iDy p2n[iD]+1*no
    int ip = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    DAT dFxx,dFxy,dFyx,dFyy;
    // initialize cache
    dFxx=dFxy=dFyx=dFyy=0.0;
    if ((ip<nmp)){
        // parallel reduction
        for(int k=0;k<nn;k++){
            dFxx+=dNx[iD]*un[iDx];
            dFxy+=dNy[iD]*un[iDx];
            dFyx+=dNx[iD]*un[iDy];
            dFyy+=dNy[iD]*un[iDy];
        }
        // incremental deformation gradient
        dF[ip+0*nmp]=1.0+dFxx;
        dF[ip+1*nmp]=0.0+dFxy;
        dF[ip+2*nmp]=0.0+dFyx;
        dF[ip+3*nmp]=1.0+dFyy;
    }
#undef iD
#undef iDx
#undef iDy
}// WRITE: 4*nmp || READ: 4*nmp*nn || TOTAL IO: 4*nmp*(1+nn)
//-----------------------------------------------------------------------//
__global__ void elastD(DAT* dF, DAT* eps, DAT* sig, DAT* ome, DAT* lp, DAT* vol, int nmp, DAT* Del, DAT dt){
    int p = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    DAT dexx,deyy,dexy,doxy,sxx0,syy0,szz0,sxy0,J;
    if (p<nmp){
        // store previous stresses
        sxx0 = sig[p+0*nmp];
        syy0 = sig[p+1*nmp];
        szz0 = sig[p+2*nmp];
        sxy0 = sig[p+3*nmp];
        // compute strains 
        dexx = (dF[p+0*nmp]-1.0        )*(DAT)1.0/(DAT)dt;  
        deyy = (dF[p+3*nmp]-1.0        )*(DAT)1.0/(DAT)dt;
        dexy = (dF[p+1*nmp]+dF[p+2*nmp])*(DAT)1.0/(DAT)dt;
        // compute spin rate
        doxy = 0.5*(dF[p+1*nmp]-dF[p+2*nmp])*(DAT)1.0/(DAT)dt;   
        // update objective stress
        sig[p+0*nmp]+= (dt*(Del[0]*dexx+Del[4]*deyy+Del[12]*dexy)+dt*(2.0*sxy0   *doxy));
        sig[p+1*nmp]+= (dt*(Del[1]*dexx+Del[5]*deyy+Del[13]*dexy)-dt*(2.0*sxy0   *doxy));
        sig[p+2*nmp]+= (dt*(Del[2]*dexx+Del[6]*deyy+Del[14]*dexy)                      );
        sig[p+3*nmp]+= (dt*(Del[3]*dexx+Del[7]*deyy+Del[15]*dexy)+dt*((syy0-sxx0)*doxy));
        // update material point volume and domain lengths
        J           = (dF[p+0*nmp]*dF[p+3*nmp])-(dF[p+1*nmp]*dF[p+2*nmp]);
        vol[p]      = J*vol[p];
        lp[p+0*nmp] = sqrt(J)*lp[p+0*nmp];
        lp[p+1*nmp] = sqrt(J)*lp[p+1*nmp];
    } 
}// WRITE: 6*nmp || READ: 8*nmp || TOTAL IO: 14*nmp
//-----------------------------------------------------------------------//
__global__ void DPPlastD(DAT* sig, DAT* cohp, DAT* phip, DAT* epII, DAT Hp, DAT cohr, DAT Kc, DAT Gc, DAT psi, int nmp){
    int p = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    DAT c,Pr,tensile,sxx,syy,szz,sxy,devxx,devyy,devzz,devxy,J2,tau,eta,etaB,xi,sigm,fs,ft,tauP,alpP,h,dlam,tauN,PrN;
    tensile = 0.0;
    if(p<nmp){
        c  = cohp[p]+(Hp*epII[p]);
        if(c<cohr){c = cohr;}
        sxx   = sig[p+0*nmp] ;
        syy   = sig[p+1*nmp] ;
        szz   = sig[p+2*nmp] ;
        sxy   = sig[p+3*nmp] ;    
        Pr    = (sxx+syy+szz)*((DAT)1.0/(DAT)3.0);
        devxx = sxx - Pr     ; 
        devyy = syy - Pr     ;
        devzz = szz - Pr     ;
        devxy = sxy          ;
        tau   = sqrt(0.5*(devxx*devxx+devyy*devyy+devzz*devzz)+devxy*devxy);       
        eta   = (3.0*  tan(phip[p]))*(DAT)1.0/((DAT)sqrt(9.0+12.0*tan(phip[p])*tan(phip[p])));
        etaB  = (3.0*  tan(psi    ))*(DAT)1.0/((DAT)sqrt(9.0+12.0*tan(psi    )*tan(psi    )));
        xi    = (3.0*c             )*(DAT)1.0/((DAT)sqrt(9.0+12.0*tan(phip[p])*tan(phip[p])));
        
        
        sigm  = min(tensile,xi*(DAT)1.0/((DAT)eta));
        sigm  = xi*(DAT)1.0/((DAT)eta)*0.0;
        fs    = tau+eta*Pr-xi;
        ft    = Pr-sigm         ;
        tauP  = xi-eta*sigm  ;
        alpP  = sqrt(1.0+eta*eta)-eta;
        h     = tau-tauP-alpP*(Pr-sigm);
        if((fs>0.0 && Pr<sigm)||(h>0.0 && Pr>=sigm)){
            dlam          = fs*(DAT)1.0/((DAT)Gc+Kc*eta*etaB)          ;
            PrN           = Pr-Kc*etaB*dlam                             ;
            tauN          = xi-eta*PrN                               ;
            sig [p+0*nmp] = devxx*((DAT)tauN/((DAT)tau))+PrN            ;
            sig [p+1*nmp] = devyy*((DAT)tauN/((DAT)tau))+PrN            ;
            sig [p+2*nmp] = devzz*((DAT)tauN/((DAT)tau))+PrN            ;
            sig [p+3*nmp] = devxy*((DAT)tauN/((DAT)tau))                ;
            epII[p      ]+= (dlam*sqrt(0.33333333+0.22222222*etaB*etaB));
        }
        if((h<=0.0)&&(Pr>=sigm)){
            dlam          = (Pr-sigm)*((DAT)1.0/((DAT)Kc))        ;
            sig[p+0*nmp] += (sigm-Pr)                             ;
            sig[p+1*nmp] += (sigm-Pr)                             ;
            sig[p+2*nmp] += (sigm-Pr)                             ;
            sig[p+3*nmp] += 0.0                                   ;
            epII[p      ]+= (sqrt(2.0)*dlam*((DAT)1.0/((DAT)3.0)));
        }  
    }
}
//-----------------------------------------------------------------------//
__global__ void MCPlastD(DAT* sig, DAT* cohp, DAT* phip, DAT* epII, DAT* iDel, DAT Hp, DAT cohr, int nmp){
    int p = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    DAT coh,phi,ds,tau,sigma,f,beta,sxx,syy,sxy,sxxn,syyn,sxyn,dsxx,dsyy,dsxy,dexx,deyy,dezz,dexy,temp;
    if(p<nmp){
        coh  = cohp[p]+(Hp*epII[p]);
        phi  = phip[p];
        sxx  = sig[p+0*nmp];
        syy  = sig[p+1*nmp];
        sxy  = sig[p+3*nmp];
        if(coh<cohr){coh = cohr;}
        ds    = sxx-syy;
        tau   = sqrt(0.25*(ds*ds)+(sxy*sxy));
        sigma = 0.5*(sxx+syy);
        f     = tau+sigma*sin(phi)-coh*cos(phi);
        if(f>0.0){
            temp = coh*((DAT)1.0/((DAT)tan(phi)));
            if(sigma<=temp){
                beta = fabs(coh*cos(phi)-sin(phi)*sigma)*((DAT)1.0/((DAT)tau));
                sxxn = sigma + 0.5*beta*ds;
                syyn = sigma - 0.5*beta*ds;
                sxyn = beta*sxy;
            }
            else if(sigma>temp){
                sxxn = temp;
                syyn = temp;
                sxyn = 0.0;
            }
            dsxx = sxxn - sxx;
            dsyy = syyn - syy;
            dsxy = sxyn - sxy;
            
            dexx = iDel[0]*dsxx+iDel[4]*dsyy+iDel[12]*dsxy;
            deyy = iDel[1]*dsxx+iDel[5]*dsyy+iDel[13]*dsxy;
            dexy = iDel[3]*dsxx+iDel[7]*dsyy+iDel[15]*dsxy;
            
            epII[p] += sqrt(0.6667*(dexx*dexx+deyy*deyy+2.0*dexy*dexy));
            
            sig[p+0*nmp] = sxxn;
            sig[p+1*nmp] = syyn;
            sig[p+3*nmp] = sxyn;
        }
    }
}
//-----------------------------------------------------------------------//
__global__ void volLockD1(DAT* pel,DAT* sig, DAT* dev, DAT* vol, int* p2e, int nmp, int nel){
#define iD  p2e[ix]
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    DAT pr,sxx,syy,szz,sxy;
    // accumulate material point pressure on element
    if (ix<nmp){
        sxx          =sig[ix+0*nmp]                     ;
        syy          =sig[ix+1*nmp]                     ;
        szz          =sig[ix+2*nmp]                     ;
        sxy          =sig[ix+3*nmp]                     ;
        pr           =-(sxx+syy+szz)*((DAT)1.0/(DAT)3.0);
        dev[ix+0*nmp]=sxx+pr                            ;
        dev[ix+1*nmp]=syy+pr                            ;
        dev[ix+2*nmp]=szz+pr                            ;
        dev[ix+3*nmp]=sxy                               ;
        atomicAdd(&pel[iD+0*nel],pr*vol[ix])            ; 
        atomicAdd(&pel[iD+1*nel],   vol[ix])            ;
    }
}// WRITE: 6*nmp || READ: 5*nmp || TOTAL IO: 11*nmp
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
    }
#undef iD
}// WRITE: 4*nmp || READ: 5*nmp || TOTAL IO: 9*nmp
//-----------------------------------------------------------------------//
//UTILITIES--------------------------------------------------------------//
DAT checkGpuMem(){
    float free_m,total_m,used_m;
    size_t free_t,total_t;
    cudaMemGetInfo(&free_t,&total_t);
    free_m =(double)free_t/1E9 ;
    total_m=(double)total_t/1E9; 
    printf("\n-----------------------------------");
    if(sizeof(DAT)==8){
        printf("\n-----------------------------------");
        printf("\nGPU memory: double arithmetic");
    }
    if(sizeof(DAT)==4){
        printf("\nGPU memory: float arithmetic");
    }
    printf("\n-----------------------------------");
    printf("\n  total: %.2f GB --> 100 %%",total_m);
    printf("\n  free : %.2f GB --> %.2f %%",free_m        ,free_m/total_m*100.0          );
    printf("\n  used : %.2f GB --> %.2f %%",total_m-free_m,(total_m-free_m)/total_m*100.0);
    return((DAT)((total_m-free_m)/total_m));   
}
//-----------------------------------------------------------------------//
void checkCudaErrors(){
    cudaError_t err = cudaGetLastError();
    if ( err != cudaSuccess )
    {
        printf("\n !! CUDA Error: %s !!\n", cudaGetErrorString(err));
        exit(-1);
    }
}
//-----------------------------------------------------------------------//
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
__global__ void basisD_TEST(DAT* N, DAT* dNx, DAT* dNy, int* p2n, int* p2e, int* e2N, DAT* xp, DAT* lp, DAT* xn, DAT dx, DAT dy, int nmp, int nn, int ney, int nel, int no, DAT xnmin, DAT ynmin){
#define iD  ix+iy*nmp
    int ix = blockIdx.x*blockDim.x + threadIdx.x; // thread ID, dimension x
    int iy = blockIdx.y*blockDim.y + threadIdx.y; // thread ID, dimension y
    DAT Nxt,Nyt,dNxt,dNyt;
    if ((ix<nmp)&&(iy<nn)){
        // Find connectivity
        p2e[ix] = (int)(floor((xp[ix+1*nmp]-ynmin)*(1.0/dx))+ney*floor((xp[ix+0*nmp]-xnmin)*(1.0/dy)));
        p2n[iD] = e2N[p2e[ix]+iy*nel];
        // Compute basis functions and derivatives
        Nxt  = getN (xp[ix+0*nmp]-xn[p2n[iD]+0*no],lp[ix+0*nmp],dx);
        dNxt = getdN(xp[ix+0*nmp]-xn[p2n[iD]+0*no],lp[ix+0*nmp],dx);
        Nyt  = getN (xp[ix+1*nmp]-xn[p2n[iD]+1*no],lp[ix+1*nmp],dy);
        dNyt = getdN(xp[ix+1*nmp]-xn[p2n[iD]+1*no],lp[ix+1*nmp],dy);   
        // Convolution of basis functions and derivatives
        N  [iD] =  Nxt* Nyt;
        dNx[iD] = dNxt* Nyt;
        dNy[iD] =  Nxt*dNyt;
    }
#undef iD
}
