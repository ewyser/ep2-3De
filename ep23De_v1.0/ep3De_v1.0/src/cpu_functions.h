// header fH.h : set of all functions call by the host
//-----------------------------------------------------------------------//
__host__ DAT CFL(DAT* vp, DAT dx, DAT dy, DAT dz, DAT yd, DAT tg, DAT tw, int nmp){
    DAT vxpmax,vypmax,vzpmax,dt;
    // determine suitable dt from CFL condition
    vxpmax = fabs(vp[0+0*nmp]);
    vypmax = fabs(vp[0+1*nmp]);
    vzpmax = fabs(vp[0+2*nmp]);
    for(int p=1;p<nmp;p++){
        if(vxpmax<fabs(vp[p+0*nmp])){
            vxpmax=fabs(vp[p+0*nmp]);
        }
        if(vypmax<fabs(vp[p+1*nmp])){
            vypmax=fabs(vp[p+1*nmp]);
        }
        if(vzpmax<fabs(vp[p+2*nmp])){
            vzpmax=fabs(vp[p+2*nmp]);
        }
    }
    dt = 0.5*min((DAT)dx/(DAT)(yd + vxpmax),(DAT)dy/(DAT)(yd + vypmax))      ;
    return dt;
}
//-----------------------------------------------------------------------//
__host__ DAT getG(DAT tw, DAT tg){
    DAT g = 0.0;
    if(tw<=tg){
        g = 9.81*tw*((DAT)1.0/(DAT)tg);
    }
    else{
        g = 9.81;
    }
    return g;
}
//-----------------------------------------------------------------------//
__host__ void topol(DAT* xp, int* p2e, int* p2n, int* e2n, DAT xnmin, DAT ynmin, DAT znmin, DAT dx, DAT dy, DAT dz, int nmp, int nn, int nex, int ney, int nez, int nel){
    for(int i=0;i<nmp;i++){
        p2e[i] = (int)((floor((xp[i+2*nmp]-znmin)*(1.0/dz)))+(nez)*floor((xp[i+0*nmp]-xnmin)*(1.0/dx))+nez*nex*(floor((xp[i+1*nmp]-ynmin)*(1.0/dy))));
        
        for(int j=0;j<nn;j++){
            p2n[i+j*nmp] = e2n[p2e[i]+j*nel];
        }
    }
}
//-----------------------------------------------------------------------//
__host__ void getNdN(DAT* N_dN,DAT xi,DAT lp,DAT dx){
    N_dN[0] = N_dN[1] = 0.0;
    if( fabs(xi)< (   lp)                     ){N_dN[0]=1.0-(4.0*(xi*xi)+((2.0*lp)*(2.0*lp)))*((DAT)1.0/((DAT)8.0*dx*lp));
                                                N_dN[1]=(-8.0*xi)*((DAT)1.0/((DAT)8.0*dx*lp));
    }
    if((fabs(xi)>=(   lp))&&(fabs(xi)<(dx-lp))){N_dN[0]=1.0-(fabs(xi)*((DAT)1.0/(DAT)dx));
                                                if(xi<0.0){N_dN[1]=( (DAT)1.0/(DAT)dx);}
                                                if(xi>0.0){N_dN[1]=(-(DAT)1.0/(DAT)dx);}
    }
    if((fabs(xi)>=(dx-lp))&&(fabs(xi)<(dx+lp))){N_dN[0]=((dx+lp-fabs(xi))*(dx+lp-fabs(xi)))*((DAT)1.0/((DAT)4.0*dx*lp));
                                                if(xi<0.0){N_dN[1]=(  dx+lp-fabs(xi))*((DAT)1.0/((DAT)2.0*dx*lp)) ;}
                                                if(xi>0.0){N_dN[1]=(-(dx+lp-fabs(xi))*((DAT)1.0/((DAT)2.0*dx*lp)));}
    }
}
//-----------------------------------------------------------------------//
__host__ void basis(DAT* xp, DAT* xn, DAT* N, DAT* dNx, DAT* dNy, DAT* dNz, int* p2n, DAT* lp, DAT dx, DAT dy, DAT dz, int nmp, int nn, int no){
    DAT Nx_dNx[2],Ny_dNy[2],Nz_dNz[2];
    int iD;
    for(int p=0;p<nmp;p++){
        for(int n=0;n<nn;n++){
            iD  = p+n*nmp;       
            // x-basis and x-derivative
            getNdN(Nx_dNx,(xp[p+0*nmp]-xn[p2n[iD]+0*no]),lp[p+0*nmp],dx);
            // y-basis and y-derivative   
            getNdN(Ny_dNy,(xp[p+1*nmp]-xn[p2n[iD]+1*no]),lp[p+1*nmp],dy);
            // y-basis and y-derivative   
            getNdN(Nz_dNz,(xp[p+2*nmp]-xn[p2n[iD]+2*no]),lp[p+2*nmp],dz);
            // convolution of basis
            N[iD]   = Nx_dNx[0]*Ny_dNy[0]*Nz_dNz[0];
            dNx[iD] = Nx_dNx[1]*Ny_dNy[0]*Nz_dNz[0];
            dNy[iD] = Nx_dNx[0]*Ny_dNy[1]*Nz_dNz[0];
            dNz[iD] = Nx_dNx[0]*Ny_dNy[0]*Nz_dNz[1];
        }
    }
}
//-----------------------------------------------------------------------//
__host__ void accum(DAT* mn, DAT* pn, DAT* fen, DAT* fin, DAT* N, DAT* dNx, DAT* dNy, DAT* dNz, DAT* mp, DAT* vp, DAT* sig, DAT* vol, int* p2n, DAT g, int nmp, int nn, int no){
    int iD, iDx, iDy, iDz;
    // initialize
    for(int n=0;n<(3*no);n++){
        if(n<no){
            mn[n] = 0.0;
        }
        pn [n] = fen[n] = fin[n] = 0.0;
    }
    // accumulate material point contributions on nodes
    for(int p=0;p<nmp;p++){
        for(int n=0;n<nn; n++){
            iD        = (p2n[p+n*nmp]);
            iDx       = iD+0*no;
            iDy       = iD+1*no;
            iDz       = iD+2*no;
            mn [iD ] += (N[p+n*nmp]*mp[p]);
            pn [iDx] += (N[p+n*nmp]*mp[p]*vp[p+0*nmp]);
            pn [iDy] += (N[p+n*nmp]*mp[p]*vp[p+1*nmp]);
            pn [iDz] += (N[p+n*nmp]*mp[p]*vp[p+2*nmp]);
            fen[iDz] -= (N[p+n*nmp]*mp[p]*g);
            fin[iDx] += vol[p]*(dNx[p+n*nmp]*sig[0+6*p]+dNy[p+n*nmp]*sig[3+6*p]+dNz[p+n*nmp]*sig[5+6*p]);
            fin[iDy] += vol[p]*(dNx[p+n*nmp]*sig[3+6*p]+dNy[p+n*nmp]*sig[1+6*p]+dNz[p+n*nmp]*sig[4+6*p]);
            fin[iDz] += vol[p]*(dNx[p+n*nmp]*sig[5+6*p]+dNy[p+n*nmp]*sig[4+6*p]+dNz[p+n*nmp]*sig[2+6*p]);
        }
    }
}
//-----------------------------------------------------------------------//
__host__ void solve(DAT* fn, DAT* fen, DAT* fin, DAT* mn, DAT* an, DAT* pn, DAT* vn, int* bc, DAT dt, int no){
    int iDx,iDy,iDz;
    DAT dmp,m,fx,fy,fz,px,py,pz,vx,vy,vz;
    // initialize
    for(int n=0;n<3*no;n++){
        fn[n] = an[n] = vn[n] = 0.0;
    }
    // solve momentum equation on the background mesh
    for(int n=0;n<no;n++){
        iDx = n+0*no;
        iDy = n+1*no;
        iDz = n+2*no;
        if(mn[n]>0.0){
            fx = (fen[iDx]-fin[iDx]);
            fy = (fen[iDy]-fin[iDy]);
            fz = (fen[iDz]-fin[iDz]);
            px = pn[iDx];
            py = pn[iDy];
            pz = pn[iDz];
            vx = px*DAT(1.0)/((DAT)mn[n]);
            vy = py*DAT(1.0)/((DAT)mn[n]);
            vz = pz*DAT(1.0)/((DAT)mn[n]);
            dmp= sqrt(fx*fx+fy*fy+fz*fz);
            fx = fx-D*dmp*((DAT)vx/(DAT)fabs(vx));
            if(fabs(vx)<1E-3){
                fx = (fen[iDx]-fin[iDx]);
            }
            fy = fy-D*dmp*((DAT)vy/(DAT)fabs(vy));
            if(fabs(vy)<1E-3){
                fy = (fen[iDy]-fin[iDy]);
            }
            fz = fz-D*dmp*((DAT)vz/(DAT)fabs(vz));
            if(fabs(vz)<1E-3){
                fz = (fen[iDz]-fin[iDz]);
            }
            m       = DAT(1.0)/((DAT)mn[n]);
            an[iDx] = (fx      )*m*bc[iDx];
            an[iDy] = (fy      )*m*bc[iDy];
            an[iDz] = (fz      )*m*bc[iDz];
            vn[iDx] = (px+dt*fx)*m*bc[iDx];
            vn[iDy] = (py+dt*fy)*m*bc[iDy];
            vn[iDz] = (pz+dt*fz)*m*bc[iDz];
        }
    }
}
//-----------------------------------------------------------------------//
__host__ void FLIP(DAT* an, DAT* vn, DAT* N, DAT* vp, DAT* xp, int* p2n, DAT dt, int nmp, int nn, int no){
    // update material point velocity
    int iDx,iDy,iDz;
    DAT dvpx,dvpy,dvpz,dxp,dyp,dzp;
    for(int p=0;p<nmp;p++){
        dvpx = dvpy = dvpz = dxp = dyp = dzp = (DAT)0.0;
        for(int n=0;n<nn;n++){
            iDx   = p2n[p+n*nmp]+0*no   ;
            iDy   = p2n[p+n*nmp]+1*no   ;
            iDz   = p2n[p+n*nmp]+2*no   ;
            dvpx += (N[p+n*nmp]*an[iDx]);
            dvpy += (N[p+n*nmp]*an[iDy]);
            dvpz += (N[p+n*nmp]*an[iDz]);
            dxp  += (N[p+n*nmp]*vn[iDx]);
            dyp  += (N[p+n*nmp]*vn[iDy]);
            dzp  += (N[p+n*nmp]*vn[iDz]);
        }
        vp[p+0*nmp] += (dt*dvpx);
        vp[p+1*nmp] += (dt*dvpy);
        vp[p+2*nmp] += (dt*dvpz);
        xp[p+0*nmp] += (dt*dxp );
        xp[p+1*nmp] += (dt*dyp );
        xp[p+2*nmp] += (dt*dzp );
    }   
}
//-----------------------------------------------------------------------//
__host__ void DM_BC(DAT* un, DAT* pn, DAT* mn, DAT* N, DAT* mp, DAT* vp, DAT* up, int* bc, DAT dt, int* p2n, int nmp, int nn, int no){
    int iDx, iDy, iDz;
    DAT m,duxp,duyp,duzp;
    // initialize nodal momentum
    for(int n=0;n<(3*no);n++){
        pn[n] = un[n] = 0.0;
    }
    // accumulate material point momentum
    for(int p=0;p<nmp;p++){
        for(int n=0;n<nn;n++){
            iDx      = p2n[p+n*nmp]+0*no             ;
            iDy      = p2n[p+n*nmp]+1*no             ;
            iDz      = p2n[p+n*nmp]+2*no             ;
            pn[iDx] += (N[p+n*nmp]*mp[p]*vp[p+0*nmp]);
            pn[iDy] += (N[p+n*nmp]*mp[p]*vp[p+1*nmp]);
            pn[iDz] += (N[p+n*nmp]*mp[p]*vp[p+2*nmp]);
        }
    }
    // solve for nodal incremental displacement
    for(int n=0;n<no;n++){
        iDx = n+0*no ;
        iDy = n+1*no ;
        iDz = n+2*no ;
        if(mn[n]>0.0){
            m       = ((DAT)1.0/(DAT)mn[n]);
            un[iDx] = dt*pn[iDx]*m*bc[iDx];
            un[iDy] = dt*pn[iDy]*m*bc[iDy];
            un[iDz] = dt*pn[iDz]*m*bc[iDz];
        }
    }
    // update material point displacement
    for(int p=0;p<nmp;p++){
        duxp = duyp = duzp = 0.0;
        for(int n=0;n<nn;n++){
            iDx   = p2n[p+n*nmp]+0*no   ;
            iDy   = p2n[p+n*nmp]+1*no   ;
            iDz   = p2n[p+n*nmp]+2*no   ;
            duxp += (N[p+n*nmp]*un[iDx]);
            duyp += (N[p+n*nmp]*un[iDy]);
            duzp += (N[p+n*nmp]*un[iDz]);
        }
        up[p+0*nmp] += (dt*duxp);
        up[p+1*nmp] += (dt*duyp);
        up[p+2*nmp] += (dt*duzp);
    } 
}
//-----------------------------------------------------------------------//   
__host__ void strains(DAT* un, DAT* dNx, DAT* dNy, DAT* dNz, DAT* dF, DAT* eps, DAT* ome, DAT* lp, DAT* vol, int* p2n, DAT dt, int nmp, int nn, int no){
    int iD,iDx,iDy,iDz;
    DAT J;
    // update material point state variables
    for(int p=0;p<nmp;p++){
        // compute incremental deformation gradient
        dF[0+p*9] = dF[4+p*9] = dF[8+p*9] = (DAT)1.0;
        dF[1+p*9] = dF[2+p*9] = dF[3+p*9] = dF[5+p*9] = dF[6+p*9] = dF[7+p*9] = (DAT)0.0;
        for(int n=0;n<nn;n++){
            iD         = p+n*nmp          ;
            iDx        = p2n[iD]+0*no     ;
            iDy        = p2n[iD]+1*no     ;
            iDz        = p2n[iD]+2*no     ;
            dF[0+p*9] += (dNx[iD]*un[iDx]);
            dF[1+p*9] += (dNy[iD]*un[iDx]);
            dF[2+p*9] += (dNz[iD]*un[iDx]);
            
            dF[3+p*9] += (dNx[iD]*un[iDy]);
            dF[4+p*9] += (dNy[iD]*un[iDy]);
            dF[5+p*9] += (dNz[iD]*un[iDy]);
            
            dF[6+p*9] += (dNx[iD]*un[iDz]);
            dF[7+p*9] += (dNy[iD]*un[iDz]);
            dF[8+p*9] += (dNz[iD]*un[iDz]);
        }
        // compute incremental strains
        eps[0+p*6] = dF[0+p*9]-1.0;
        eps[1+p*6] = dF[4+p*9]-1.0;
        eps[2+p*6] = dF[8+p*9]-1.0;
        eps[3+p*6] = dF[1+p*9]+dF[3+p*9];
        eps[4+p*6] = dF[7+p*9]+dF[5+p*9];
        eps[5+p*6] = dF[6+p*9]+dF[2+p*9];
        // compute incremental rotation
        ome[0+p*3]     = 0.5*(dF[3+p*9]-dF[1+p*9]);
        ome[1+p*3]     = 0.5*(dF[5+p*9]-dF[7+p*9]);
        ome[2+p*3]     = 0.5*(dF[2+p*9]-dF[6+p*9]);
        // update material point volume and domain lengths
        J          = (DAT)1.0+(eps[0+p*6]+eps[1+p*6]+eps[2+p*6]);
        vol[p]     = J*vol[p];
        J          = pow(J,(DAT)0.3333);
        lp[p+0*nmp] = J*lp[p+0*nmp];
        lp[p+1*nmp] = J*lp[p+1*nmp];
        lp[p+2*nmp] = J*lp[p+2*nmp];
    }
}
//-----------------------------------------------------------------------//
__host__ void elast(DAT* sig, DAT* eps, DAT* ome, DAT* Del, int nmp, DAT dt){
    DAT sxx0,syy0,szz0,sxy0,syz0,sxz0,dexx,deyy,dezz,dexy,deyz,dexz,doxy,doyz,doxz;
    for(int p=0;p<nmp;p++){
        // store previous stresses
        sxx0   = sig[0+6*p];
        syy0   = sig[1+6*p];
        szz0   = sig[2+6*p];
        sxy0   = sig[3+6*p];
        syz0   = sig[4+6*p];
        sxz0   = sig[5+6*p];
        // strain rate
        dexx   = eps[0+6*p]/(DAT)dt;
        deyy   = eps[1+6*p]/(DAT)dt;
        dezz   = eps[2+6*p]/(DAT)dt;
        dexy   = eps[3+6*p]/(DAT)dt;
        deyz   = eps[4+6*p]/(DAT)dt;
        dexz   = eps[5+6*p]/(DAT)dt;
        // spin
        doxy   = ome[0+3*p]/dt;
        doyz   = ome[1+3*p]/dt;
        doxz   = ome[2+3*p]/dt;
        // update objective stress        
        sig[0+6*p] +=  (DAT)2*dt*(sxy0*doxy+sxz0*doxz);
        sig[1+6*p] += -(DAT)2*dt*(sxy0*doxy-syz0*doyz);
        sig[2+6*p] += -(DAT)2*dt*(sxz0*doxz+syz0*doyz);
        sig[3+6*p] +=         dt*(doxy*(syy0-sxx0)+syz0*doxz+sxz0*doyz);
        sig[4+6*p] +=         dt*(doyz*(szz0-syy0)-sxy0*doxz-sxz0*doxy);
        sig[5+6*p] +=         dt*(doxz*(szz0-sxx0)+syz0*doxy-sxy0*doyz);
        // incremental strain
        sig[0+6*p] += dt*(Del[0]*dexx+Del[6 ]*deyy+Del[12]*dezz+Del[18]*dexy+Del[24]*deyz+Del[30]*dexz);
        sig[1+6*p] += dt*(Del[1]*dexx+Del[7 ]*deyy+Del[13]*dezz+Del[19]*dexy+Del[25]*deyz+Del[31]*dexz);
        sig[2+6*p] += dt*(Del[2]*dexx+Del[8 ]*deyy+Del[14]*dezz+Del[20]*dexy+Del[26]*deyz+Del[32]*dexz);
        sig[3+6*p] += dt*(Del[3]*dexx+Del[9 ]*deyy+Del[15]*dezz+Del[21]*dexy+Del[27]*deyz+Del[33]*dexz);
        sig[4+6*p] += dt*(Del[4]*dexx+Del[10]*deyy+Del[16]*dezz+Del[22]*dexy+Del[28]*deyz+Del[34]*dexz);
        sig[5+6*p] += dt*(Del[5]*dexx+Del[11]*deyy+Del[17]*dezz+Del[23]*dexy+Del[29]*deyz+Del[35]*dexz);   
    }
}
//-----------------------------------------------------------------------//
__host__ void DPPlast(DAT* sig, DAT* cohp, DAT* phip, DAT* epII, DAT Hp, DAT cohr, DAT Kc, DAT Gc, DAT psi, int nmp){
    DAT c,Pr,tensile,sxx,syy,szz,sxy,syz,sxz,devxx,devyy,devzz,devxy,devyz,devxz,tau,eta,etaB,xi,sigm,fs,ft,tauP,alpP,h,dlam,tauN,PrN,dep;
    tensile = (DAT)0.0;
    for(int p=0; p<nmp; p++){
        c  = cohp[p]+(Hp*epII[p]);
        if(c<cohr){c = cohr;}
        sxx   = sig[0+6*p] ;
        syy   = sig[1+6*p] ;
        szz   = sig[2+6*p] ;
        sxy   = sig[3+6*p] ;
        syz   = sig[4+6*p] ;
        sxz   = sig[5+6*p] ;
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
    
        
        sigm  = min(tensile,xi*(DAT)1.0/((DAT)eta));
        sigm  = xi*(DAT)1.0/((DAT)eta);
        fs    = tau+eta*Pr-xi;
        ft    = Pr-sigm         ;
        tauP  = xi-eta*sigm  ;
        alpP  = sqrt((DAT)1.0+eta*eta)-eta;
        h     = tau-tauP-alpP*(Pr-sigm);
        if((fs>(DAT)0.0 && Pr<sigm)||(h>(DAT)0.0 && Pr>=sigm)){
            dlam        = fs*(DAT)1.0/((DAT)Gc+Kc*eta*etaB)           ;
            PrN         = Pr-Kc*etaB*dlam                             ;
            tauN        = xi-eta*PrN                                  ;
            sig [0+6*p] = devxx*((DAT)tauN/((DAT)tau))+PrN            ;
            sig [1+6*p] = devyy*((DAT)tauN/((DAT)tau))+PrN            ;
            sig [2+6*p] = devzz*((DAT)tauN/((DAT)tau))+PrN            ;
            sig [3+6*p] = devxy*((DAT)tauN/((DAT)tau))                ;
            sig [4+6*p] = devyz*((DAT)tauN/((DAT)tau))                ;
            sig [5+6*p] = devxz*((DAT)tauN/((DAT)tau))                ;
            dep         = (dlam*sqrt((DAT)0.33333333+(DAT)0.22222222*etaB*etaB));
            epII[p    ]+= dep                                         ;
        }
        if((h<=0.0)&&(Pr>=sigm)){
            dlam        = (Pr-sigm)*((DAT)1.0/((DAT)Kc))        ;
            sig[0+6*p] += (sigm-Pr)                             ;
            sig[1+6*p] += (sigm-Pr)                             ;
            sig[2+6*p] += (sigm-Pr)                             ;
            sig[3+6*p] += (DAT)0.0                              ;
            sig[4+6*p] += (DAT)0.0                              ;
            sig[5+6*p] += (DAT)0.0                              ;
            dep         = (sqrt(2.0)*dlam*((DAT)1.0/((DAT)3.0)));
            epII[p    ]+= dep                                   ;
        }  
    }
}
//-----------------------------------------------------------------------//
__host__ void volLock(DAT* pel,DAT* sig, DAT* dev, DAT* vol, int* p2e, int nmp, int nel){
    int iD;
    DAT pr;
    // initialize element pressure
    for(unsigned int e=0;e<nel;e++){
        pel[0+2*e] = 0.0;
        pel[1+2*e] = 0.0;
    }
    // accumulate material point pressure on element
    for(unsigned int p=0;p<nmp;p++){
        iD          = p2e[p]                      ;
        pr          = -(sig[0+6*p]+sig[1+6*p]+sig[2+6*p])/(DAT)3.0;
        dev[0+6*p ] = sig[0+6*p]+pr               ;
        dev[1+6*p ] = sig[1+6*p]+pr               ;
        dev[2+6*p ] = sig[2+6*p]+pr               ;
        dev[3+6*p ] = sig[3+6*p]+0.0              ;
        dev[4+6*p ] = sig[4+6*p]+0.0              ;
        dev[5+6*p ] = sig[5+6*p]+0.0              ;
        pel[0+2*iD]+= pr*vol[p]                   ;
        pel[1+2*iD]+= vol[p]                      ;
    }
    // assign element pressure to material point
    for(unsigned int p=0;p<nmp;p++){
        iD          = p2e[p]                      ;
        pr          = pel[0+2*iD]*1.0/pel[1+2*iD] ;
        sig[0+6*p]  = dev[0+6*p]-pr               ;
        sig[1+6*p]  = dev[1+6*p]-pr               ;
        sig[2+6*p]  = dev[2+6*p]-pr               ;
        sig[3+6*p]  = dev[3+6*p]                  ;
        sig[4+6*p]  = dev[4+6*p]                  ;
        sig[5+6*p]  = dev[5+6*p]                  ;
    }
}
//-----------------------------------------------------------------------//