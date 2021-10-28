#include "stdio.h"
#include "cuda.h"
#include "time.h"
#include "macros.h"
#include "cpu_functions.h"
#include "gpu_kernels.h"
int main(){
    // Set up GPU
    int  gpu_id=0;
    dim3 Grid, Block;
    dim3 grid, block;
    cudaSetDevice(gpu_id); cudaGetDevice(&gpu_id);
    cudaDeviceReset(); cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);  // set L1 to prefered
    // Timer
    int saveData=0;
    cudaEvent_t startD, stopD;
    float milliseconds = (DAT)0.0;
    // Get initial GPU memory load
    DAT GPUinfo[4];        
    double total_m,free_m,used_m;
    size_t free_t,total_t;
    cudaMemGetInfo(&free_t,&total_t);
    free_m = free_t          /(double)(1024.0*1024.0*1024.0);
    total_m= total_t         /(double)(1024.0*1024.0*1024.0);
    used_m = (total_t-free_t)/(double)(1024.0*1024.0*1024.0);
    GPUinfo[3]=(DAT)((total_m-free_m)/total_m);       
    // Import
    #include "I_O.h"
    load(bcs,3*no,1,"bcs.dat",int);//load(bcN,3*no,1,"bcN.dat",int);
    // CUDA
    Block.x = BCK1;
    Block.y = nn;
    Grid.x  = ceil(max(nmp,no)/Block.x)+1;
    Grid.y  = 1;
    block.x = BCK2;
    block.y = 1;
    grid.x  = ceil(max(nmp,no)/block.x)+1;
    grid.y  = 1;
    // Solver
    int flag = 0;
    DAT DT   = 0.0;
    printf("\nGPU set up:");
    printf("\n-----------------------------------");
    printf("\n  Grid = %d, Block = %d",Grid.x*Grid.y,Block.x*Block.y);
    printf("\n  grid = %d, block = %d",grid.x*grid.y,block.x*block.y);        
    cudaEventCreate(&startD);
    cudaEventCreate(&stopD);
    cudaEventRecord(startD);
    // Check Mem Usage
    GPUinfo[3]=checkGpuMem()-GPUinfo[3];
    //checkCudaErrors();
    while(tw<t){
        if(tw>te && flag == 0){
            H2D(up ,3,nmp,DAT);
            D2H(xp ,3,nmp,DAT);
            save(xp,3*nmp,1,"xp0.dat");
            //swapBCs<<<grid,block>>>(bcN_d,bcs_d,no);
            flag++;
        }
        // get adaptative dt
        D2H(vp ,3,nmp,DAT);
        dt = CFL(vp_h,dx,dy,dz,yd,tg,tw,nmp);
        // linear gravitational increase
        g  = getG(tw,tg);
        //              
        initD<<<grid,block>>>(mn_d,pn_d,fen_d,fin_d,pel_d,no,nel);
        basisD<<<Grid,Block>>>(N_d,dNx_d,dNy_d,dNz_d,p2n_d,p2e_d,e2n_d,xp_d,lp_d,xn_d,dx,dy,dz,nmp,nn,nex,ney,nez,nel,no,xnmin,ynmin,znmin);
        cudaDeviceSynchronize();
        accumD<<<Grid,Block>>>(mn_d,pn_d,fen_d,fin_d,N_d,dNx_d,dNy_d,dNz_d,sig_d,vp_d,mp_d,vol_d,p2n_d,g,nmp,nn,no);
        cudaDeviceSynchronize();
        solveD<<<grid,block>>>(fen_d,fin_d,mn_d,an_d,pn_d,vn_d,un_d,dt,no,bcs_d);
        cudaDeviceSynchronize();
        projectD<<<grid,block>>>(an_d,vn_d,N_d,vp_d,xp_d,up_d,p2n_d,dt,nmp,nn,no);
        cudaDeviceSynchronize();
        DMD<<<Grid,Block>>>(mn_d,un_d,N_d,vp_d,mp_d,p2n_d,dt,nmp,nn,bcs_d,no);
        cudaDeviceSynchronize();
        getdFD<<<grid,block>>>(un_d,dF_d,dNx_d,dNy_d,dNz_d,p2n_d,nmp,nn,no);
        cudaDeviceSynchronize();
        elastD<<<grid,block>>>(dF_d,eps_d,sig_d,ome_d,lp_d,vol_d,nmp,Del_d,dt);
        cudaDeviceSynchronize();
        if(tw>te){
            DPPlastD<<<grid,block>>>(sig_d,cohp_d,phip_d,epII_d,Hp,cohr,Kc,Gc,psi0,nmp);
            cudaDeviceSynchronize();        
        }
        volLockD1<<<grid,block>>>(pel_d,sig_d,dev_d,vol_d,p2e_d,nmp,nel);
        cudaDeviceSynchronize();
        volLockD2<<<grid,block>>>(pel_d,sig_d,dev_d,p2e_d,nmp,nel);
        cudaDeviceSynchronize();
        // update time & iteration
        tw+=dt;
        DT+=dt;
        it++;
        // GPU workload & save data
        if(((int)FPS>(int)0)&&(DT>=((DAT)1.0/(DAT)FPS))){
            // save data
            if(tw>te){
                saveData++;
                D2HD(xp  ,3,nmp,"xp"  ,DAT,saveData,SIM);
                D2HD(up  ,3,nmp,"up"  ,DAT,saveData,SIM);
                D2HD(lp  ,3,nmp,"lp"  ,DAT,saveData,SIM);
                D2HD(sig ,6,nmp,"sig" ,DAT,saveData,SIM);
                D2HD(epII,1,nmp,"epII",DAT,saveData,SIM);  
            }  
            // display workload
            printf("\n  workload = %.2f %%",100*tw/t);
            // reset time interval counter DT
            DT = 0.0;
        }
    }
    cudaEventRecord(stopD);
    cudaEventSynchronize(stopD);
    milliseconds = (DAT)0.0;
    cudaEventElapsedTime(&milliseconds, startD, stopD);
    
    GPUinfo[0] = (DAT)1E-3*milliseconds;
    GPUinfo[1] = it/GPUinfo[0];
    GPUinfo[2] = (IO*it*sizeof(DAT)/(1024*1024*1024*GPUinfo[0]));
    printf("\n-----------------------------------");
    printf("\nGPU summary: MTPeff = %.2f [GB/s]",GPUinfo[2]);
    printf("\n-----------------------------------");
    printf("\n  time is %.2f s\n  after %d iterations \n  i.e., %.2f it/s\n",GPUinfo[0],it,GPUinfo[1]);
    // save data
    D2H(epII,1,nmp,DAT);
    D2H(xp  ,3,nmp,DAT);
    D2H(up  ,3,nmp,DAT);
    D2H(lp  ,3,nmp,DAT);
    D2H(sig ,6,nmp,DAT);
    D2H(un  ,3,no ,DAT);
    save(epII,nmp  ,1,"epII.dat");
    save(xp  ,3*nmp,1,"xp.dat"  );
    save(lp  ,3*nmp,1,"lp.dat"  );
    save(up  ,3*nmp,1,"up.dat"  );
    save(sig ,6*nmp,1,"sig.dat" );
    FILE* fidw=fopen("GPUinfo.dat", "wb");
    fwrite(GPUinfo, 4*sizeof(DAT), 1, fidw);
    fclose(fidw);
    FILE *fidw2=fopen("numSaved.dat","wb");
    fwrite(&saveData,1,sizeof(saveData),fidw2);
    fclose(fidw2);
    // clear host memory & clear device memory
    //free_all(pel);
    //free_all(mn);
    //free_all(pn);
    //free_all(fen);
    //free_all(fin);
    //free_all(an);
    //free_all(vn);
    //free_all(un);
    //free_all(p2e);
    //free_all(p2n);
    //free_all(sig);
    //free_all(epII);
    //free_all(N);
    //free_all(dNx);
    //free_all(dNy);
    //free_all(dNz);
    //free_all(vp);
    //free_all(up);
    //free_all(dF);
    //free_all(eps);
    //free_all(ome);
    //free_all(dev);
    //free_all(bcs);
    cudaDeviceReset();
}
