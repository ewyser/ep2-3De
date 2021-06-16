#include "stdio.h"
#include "cuda.h"
#include "time.h"
#include "macros.h"
#include "cpu_functions.h"

int main(){
    // Set up GPU
    int  gpu_id=-1;
    dim3 Grid, Block; 
    cudaSetDevice(gpu_id); cudaGetDevice(&gpu_id);
    cudaDeviceReset(); cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);  // set L1 to prefered
    // Timer       
    clock_t start, end;
    double CPUinfo[2];
    // Import
    #include "I_O.h"
    load(bcs,3*no,1,"bcs.dat",int);
    // Solver 
    start = clock();
    while(tw<t){
        // get adaptative dt
        dt = CFL(vp_h,dx,dy,dz,yd,tg,tw,nmp);
        // linear gravitational increase
        g  = getG(tw,tg);
        // 
        topol(xp_h,p2e_h,p2n_h,e2n_h,xnmin,ynmin,znmin,dx,dy,dz,nmp,nn,nex,ney,nez,nel);
        basis(xp_h,xn_h,N_h,dNx_h,dNy_h,dNz_h,p2n_h,lp_h,dx,dy,dz,nmp,nn,no);
        accum(mn_h,pn_h,fen_h,fin_h,N_h,dNx_h,dNy_h,dNz_h,mp_h,vp_h,sig_h,vol_h,p2n_h,g,nmp,nn,no);
        solve(fn_h,fen_h,fin_h,mn_h,an_h,pn_h,vn_h,bcs_h,dt,no);
        FLIP(an_h,vn_h,N_h,vp_h,xp_h,p2n_h,dt,nmp,nn,no);
        DM_BC(un_h,pn_h,mn_h,N_h,mp_h,vp_h,up_h,bcs_h,dt,p2n_h,nmp,nn,no);
        strains(un_h,dNx_h,dNy_h,dNz_h,dF_h,eps_h,ome_h,lp_h,vol_h,p2n_h,dt,nmp,nn,no);
        elast(sig_h,eps_h,ome_h,Del_h,nmp,dt);
        if(tw>te){
            DPPlast(sig_h,cohp_h,phip_h,epII_h,Hp,cohr,Kc,Gc,psi0,nmp);
        }
        volLock(pel_h,sig_h,dev_h,vol_h,p2e_h,nmp,nel);
        // update time & iteration
        tw+=dt;
        it++;
    }
    end = clock();
    CPUinfo[0] = (double)(((double)(end-start))/CLOCKS_PER_SEC);
    CPUinfo[1] = it/CPUinfo[0];
    CPUinfo[2] = (IO*it*sizeof(DAT)/(1024*1024*1024*CPUinfo[0]));    
    printf("\n-----------------------------------");
    printf("\nCPU summary: MTPeff = %.2f [GB/s]",CPUinfo[2]);
    printf("\n-----------------------------------");  
    printf("\n  CPU time is %.2f s\n  after %d iterations \n  i.e., %.2f it/s\n",CPUinfo[0],it,CPUinfo[1]);
    // save data
    save(epII,1*nmp,1,"C_epII.dat");
    save(xp  ,3*nmp,1,"C_xp.dat"  );
    save(lp  ,3*nmp,1,"C_lp.dat"  );
    save(up  ,3*nmp,1,"C_up.dat"  );
    save(sig ,6*nmp,1,"C_sig.dat" );

FILE* fidw=fopen("CPUinfo.dat", "wb");
fwrite(CPUinfo, 4*sizeof(DAT), 1, fidw);
fclose(fidw);

    cudaDeviceReset();
}
