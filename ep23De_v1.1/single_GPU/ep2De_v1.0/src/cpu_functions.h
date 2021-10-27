// header fH.h : set of all functions call by the host
//-----------------------------------------------------------------------//
__host__ DAT CFL(DAT* vp, DAT dx, DAT dy, DAT yd, DAT tg, DAT tw, int nmp){
    DAT vxpmax,vypmax,dt;
    // determine suitable dt from CFL condition
    vxpmax = fabs(vp[0+0*nmp]);
    vypmax = fabs(vp[0+1*nmp]);
    for(int p=1;p<nmp;p++){
        if(vxpmax<fabs(vp[p+0*nmp])){
            vxpmax=fabs(vp[p+0*nmp]);
        }
        if(vypmax<fabs(vp[p+1*nmp])){
            vypmax=fabs(vp[p+1*nmp]);
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
__host__ void topol(DAT* xp, int* p2e, int* p2n, int* e2n, DAT xnmin, DAT ynmin, DAT dx, DAT dy, int nmp, int nn, int ney, int nel){
    for(int i=0;i<nmp;i++){
        p2e[i] = (int)((floor((xp[i+1*nmp]-ynmin)*((DAT)1.0/(DAT)dx)))+ney*floor((xp[i+0*nmp]-xnmin)*((DAT)1.0/(DAT)dy)));
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
__host__ void basis(DAT* xp, DAT* xn, DAT* N, DAT* dNx, DAT* dNy, int* p2n, DAT* lp, DAT dx, DAT dy, int nmp, int nn, int no){
    DAT Nx_dNx[2],Ny_dNy[2];
    int iD;
    for(int p=0;p<nmp;p++){
        for(int n=0;n<nn;n++){
            iD  = p+n*nmp;       
            // x-basis and x-derivative
            getNdN(Nx_dNx,(xp[p+0*nmp]-xn[p2n[iD]+0*no]),lp[p+0*nmp],dx);
            // y-basis and y-derivative   
            getNdN(Ny_dNy,(xp[p+1*nmp]-xn[p2n[iD]+1*no]),lp[p+1*nmp],dy);
            // convolution of basis
            N[iD]   = Nx_dNx[0]*Ny_dNy[0];
            dNx[iD] = Nx_dNx[1]*Ny_dNy[0];
            dNy[iD] = Nx_dNx[0]*Ny_dNy[1];
        }
    }
}
//-----------------------------------------------------------------------//
__host__ void accum(DAT* mn, DAT* pn, DAT* fen, DAT* fin, DAT* N, DAT* dNx, DAT* dNy, DAT* mp, DAT* vp, DAT* sig, DAT* vol, int* p2n, DAT g, int nmp, int nn, int no){
    int iD, iDx, iDy;
    // initialize
    for(int n=0;n<(2*no);n++){
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
            mn [iD ] += (N[p+n*nmp]*mp[p]);
            pn [iDx] += (N[p+n*nmp]*mp[p]*vp[p+0*nmp]);
            pn [iDy] += (N[p+n*nmp]*mp[p]*vp[p+1*nmp]);
            fen[iDy] -= (N[p+n*nmp]*mp[p]*g);
            fin[iDx] += (vol[p]*(dNx[p+n*nmp]*sig[0+p*4]+dNy[p+n*nmp]*sig[3+p*4]));
            fin[iDy] += (vol[p]*(dNx[p+n*nmp]*sig[3+p*4]+dNy[p+n*nmp]*sig[1+p*4]));
        }
    }
}
//-----------------------------------------------------------------------//
__host__ void solve(DAT* fn, DAT* fen, DAT* fin, DAT* mn, DAT* an, DAT* pn, DAT* vn, int* bc, DAT dt, int no){
    int iDx,iDy;
    DAT dmp,m,fx,fy,px,py,vx,vy;
    // initialize
    for(int n=0;n<2*no;n++){
        fn[n] = an[n] = vn[n] = 0.0;
    }
    // solve momentum equation on the background mesh
    for(int n=0;n<no;n++){
        iDx = n+0*no;
        iDy = n+1*no;
        if(mn[n]>0.0){
            fx = (fen[iDx]-fin[iDx]);
            fy = (fen[iDy]-fin[iDy]);
            px = pn[iDx];
            py = pn[iDy];
            vx = px*DAT(1.0)/((DAT)mn[n]);
            vy = py*DAT(1.0)/((DAT)mn[n]);
            dmp= sqrt(fx*fx+fy*fy);
            fx = fx-D*dmp*((DAT)vx/(DAT)fabs(vx));
            if(fabs(vx)<1E-3){
                fx = (fen[iDx]-fin[iDx]);
            }
            fy = fy-D*dmp*((DAT)vy/(DAT)fabs(vy));
            if(fabs(vy)<1E-3){
                fy = (fen[iDy]-fin[iDy]);
            }
            m       = DAT(1.0)/((DAT)mn[n]);
            an[iDx] = (fx      )*m*bc[iDx];
            an[iDy] = (fy      )*m*bc[iDy];
            vn[iDx] = (px+dt*fx)*m*bc[iDx];
            vn[iDy] = (py+dt*fy)*m*bc[iDy];
        }
    }
}
//-----------------------------------------------------------------------//
__host__ void FLIP(DAT* an, DAT* vn, DAT* N, DAT* vp, DAT* xp, int* p2n, DAT dt, int nmp, int nn, int no){
    // update material point velocity
    int iDx,iDy;
    DAT dvpx,dvpy,dxp,dyp;
    for(int p=0;p<nmp;p++){
        dvpx = dvpy = dxp  = dyp  = 0.0;
        for(int n=0;n<nn;n++){
            iDx   = p2n[p+n*nmp]+0*no   ;
            iDy   = p2n[p+n*nmp]+1*no   ;
            dvpx += (N[p+n*nmp]*an[iDx]);
            dvpy += (N[p+n*nmp]*an[iDy]);
            dxp  += (N[p+n*nmp]*vn[iDx]);
            dyp  += (N[p+n*nmp]*vn[iDy]);
        }
        vp[p+0*nmp] += (dt*dvpx);
        vp[p+1*nmp] += (dt*dvpy);
        xp[p+0*nmp] += (dt*dxp );
        xp[p+1*nmp] += (dt*dyp );
    }   
}
//-----------------------------------------------------------------------//
__host__ void PICFLIP2(DAT* an, DAT* vn, DAT* N, DAT* vp, DAT* xp, int* p2n, DAT dt, int nmp, int nn, int no){
    // update material point velocity
    int iDx,iDy;
    DAT dapx,dapy,dvpx,dvpy,vxp,vyp;
        for(int p=0;p<nmp;p++){
        dapx = dapy = dvpx = dvpy = 0.0;
        for(int n=0;n<nn;n++){
            iDx   = p2n[p+n*nmp]+0*no   ;
            iDy   = p2n[p+n*nmp]+1*no   ;
            dapx += (N[p+n*nmp]*an[iDx]);// FLIP update
            dapy += (N[p+n*nmp]*an[iDy]);// FLIP update
            dvpx += (N[p+n*nmp]*vn[iDx]);// PIC update
            dvpy += (N[p+n*nmp]*vn[iDy]);// PIC update
        }
        vxp          = vp[p+0*nmp];
        vyp          = vp[p+1*nmp];
        vp[p+0*nmp] += (dt*dapx)-alpha*(vxp-dvpx);
        vp[p+1*nmp] += (dt*dapy)-alpha*(vyp-dvpy);
        xp[p+0*nmp] += (dvpx-0.5*alpha*(vxp-dvpx))*dt+0.5*(dapx*dt*dt);
        xp[p+1*nmp] += (dvpy-0.5*alpha*(vyp-dvpy))*dt+0.5*(dapy*dt*dt);
    }  
}
//-----------------------------------------------------------------------//
__host__ void DM_BC(DAT* un, DAT* pn, DAT* mn, DAT* N, DAT* mp, DAT* vp, DAT* up, int* bc, DAT dt, int* p2n, int nmp, int nn, int no){
    int iDx, iDy;
    DAT m,duxp,duyp;
    // initialize nodal momentum
    for(int n=0;n<(2*no);n++){
        pn[n] = un[n] = 0.0;
    }
    // accumulate material point momentum
    for(int p=0;p<nmp;p++){
        for(int n=0;n<nn;n++){
            iDx      = p2n[p+n*nmp]+0*no             ;
            iDy      = p2n[p+n*nmp]+1*no             ;
            pn[iDx] += (N[p+n*nmp]*mp[p]*vp[p+0*nmp]);
            pn[iDy] += (N[p+n*nmp]*mp[p]*vp[p+1*nmp]);
        }
    }
    // solve for nodal incremental displacement
    for(int n=0;n<no;n++){
        iDx = n+0*no ;
        iDy = n+1*no ;
        if(mn[n]>0.0){
            m       = ((DAT)1.0/(DAT)mn[n]);
            un[iDx] = dt*pn[iDx]*m*bc[iDx];
            un[iDy] = dt*pn[iDy]*m*bc[iDy];
        }
    }
    // update material point displacement
    for(int p=0;p<nmp;p++){
        duxp = duyp = 0.0;
        for(int n=0;n<nn;n++){
            iDx   = p2n[p+n*nmp]+0*no   ;
            iDy   = p2n[p+n*nmp]+1*no   ;
            duxp += (N[p+n*nmp]*un[iDx]);
            duyp += (N[p+n*nmp]*un[iDy]);
        }
        up[p+0*nmp] += (dt*duxp);
        up[p+1*nmp] += (dt*duyp);
    } 
}
//-----------------------------------------------------------------------//   
__host__ void strains(DAT* un, DAT* dNx, DAT* dNy, DAT* dF, DAT* eps, DAT* ome, DAT* lp, DAT* vol, int* p2n, int nmp, int nn, int no){
    int iD,iDx,iDy;
    DAT J;
    // update material point state variables
    for(int p=0;p<nmp;p++){
        // compute incremental deformation gradient
        dF[0+p*4] = dF[3+p*4] = 1.0;
        dF[1+p*4] = dF[2+p*4] = 0.0;
        for(int n=0;n<nn;n++){
            iD         = p+n*nmp          ;
            iDx        = p2n[iD]+0*no     ;
            iDy        = p2n[iD]+1*no     ;
            dF[0+p*4] += (dNx[iD]*un[iDx]);
            dF[1+p*4] += (dNy[iD]*un[iDx]);
            dF[2+p*4] += (dNx[iD]*un[iDy]);
            dF[3+p*4] += (dNy[iD]*un[iDy]);
        }
        // compute incremental strains
        eps[0+p*3] = dF[0+p*4]-1.0;
        eps[1+p*3] = dF[3+p*4]-1.0;
        eps[2+p*3] = dF[1+p*4]+dF[2+p*4];
        // compute incremental rotation
        ome[p]     = 0.5*(dF[1+p*4]-dF[2+p*4]);
        // update material point volume and domain lengths
        J          = (dF[0+p*4]*dF[3+p*4])-(dF[1+p*4]*dF[2+p*4]);
        vol[p]     = J*vol[p];
        lp[p+0*nmp] = sqrt(J)*lp[p+0*nmp];
        lp[p+1*nmp] = sqrt(J)*lp[p+1*nmp];
    }
}
//-----------------------------------------------------------------------//
__host__ void elast(DAT* sig, DAT* eps, DAT* ome, DAT* Del, int nmp, DAT dt){
    DAT sxx,syy,sxy,dexx,deyy,dexy,doxy;
    for(int p=0;p<nmp;p++){
        // store previous stresses
        sxx    = sig[0+4*p];
        syy    = sig[1+4*p];
        sxy    = sig[3+4*p];
        // strain rate
        dexx   = eps[0+3*p]/dt;
        deyy   = eps[1+3*p]/dt;
        dexy   = eps[2+3*p]/dt;
        // spin
        doxy   = ome[p]/dt;
        // update stress
        sig[0+4*p]+= dt*(Del[0]*dexx+Del[4]*deyy+Del[12]*dexy)+dt*2.0*sxy*doxy  ;
        sig[1+4*p]+= dt*(Del[1]*dexx+Del[5]*deyy+Del[13]*dexy)-dt*2.0*sxy*doxy  ;
        sig[2+4*p]+= dt*(Del[2]*dexx+Del[6]*deyy+Del[14]*dexy)                  ;
        sig[3+4*p]+= dt*(Del[3]*dexx+Del[7]*deyy+Del[15]*dexy)+dt*(syy-sxx)*doxy;
    }
}
//-----------------------------------------------------------------------//
__host__ void plast(DAT* sig, DAT* cohp, DAT* phip, DAT* epII, DAT* iDel, DAT Hp, DAT cohr, int nmp){
    DAT coh,phi,ds,tau,sigma,f,beta,sxx,syy,sxy,sxxn,syyn,sxyn,dsxx,dsyy,dsxy,dexx,deyy,dexy,temp;
    for(int p=0;p<nmp;p++){
        coh  = cohp[p]+(Hp*epII[p]);
        phi  = phip[p];
        sxx  = sig[0+4*p];
        syy  = sig[1+4*p];
        sxy  = sig[3+4*p];        
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
                
                sig[0+4*p] = sxxn;
                sig[1+4*p] = syyn;
                sig[3+4*p] = sxyn;
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
        pr          = -(sig[0+4*p]+sig[1+4*p]+sig[2+4*p])/(DAT)3.0;
        dev[0+4*p ] = sig[0+4*p]+pr               ;
        dev[1+4*p ] = sig[1+4*p]+pr               ;
        dev[2+4*p ] = sig[2+4*p]+pr               ;
        dev[3+4*p ] = sig[3+4*p]+0.0              ;
        pel[0+2*iD]+= pr*vol[p]                   ;
        pel[1+2*iD]+= vol[p]                      ;
    }
    // assign element pressure to material point
    for(unsigned int p=0;p<nmp;p++){
        iD          = p2e[p]                      ;
        pr          = pel[0+2*iD]*1.0/pel[1+2*iD] ;
        sig[0+4*p]  = dev[0+4*p]-pr               ;
        sig[1+4*p]  = dev[1+4*p]-pr               ;
        sig[2+4*p]  = dev[2+4*p]-pr               ;
        sig[3+4*p]  = dev[3+4*p]                  ;
    }
}
//-----------------------------------------------------------------------//



__host__ void elast_TEST(DAT* sig, DAT* eps, DAT* ome, DAT* Del, int nmp, DAT dt){
    DAT sxx,syy,sxy,cos2,sin2,sincos,sum;
    for(int p=0;p<nmp;p++){
        cos2   = cos(-ome[p])*cos(-ome[p]);
        sin2   = sin(-ome[p])*sin(-ome[p]);
        sincos = sin(-ome[p])*cos(-ome[p]);
        sxx    = sig[0+3*p];
        syy    = sig[1+3*p];
        sxy    = sig[2+3*p];
        // objective stress
        sig[0+3*p] = sxx*cos2+syy*sin2-2*sxy*sincos;
        sig[1+3*p] = sxx*sin2+syy*cos2+2*sxy*sincos;
        sig[2+3*p] = (sxx-syy)*sincos+sxy*(cos2-sin2);
        // update stress
        for(int N=0; N<3; N++){
            sum = 0.0;
            for(int n=0; n<3; n++){
                sum+=(Del[N+3*n]*eps[n+3*p]);
            }
            sig[N+3*p]+= sum;
        }
//         sig[0+3*p]+= (Del[0]*eps[0+3*p]+Del[3]*eps[1+3*p]+Del[6]*eps[2+3*p]);
//         sig[1+3*p]+= (Del[1]*eps[0+3*p]+Del[4]*eps[1+3*p]+Del[7]*eps[2+3*p]);
//         sig[2+3*p]+= (Del[2]*eps[0+3*p]+Del[5]*eps[1+3*p]+Del[8]*eps[2+3*p]);
    }
}