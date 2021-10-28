#include "stdio.h"
#include "stdlib.h"
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
    DAT GPUinfo[3];
    cudaEvent_t startD, stopD;
    float milliseconds = 0;
    // Import
    #include "I_O.h"
    load(bcs,2*no,1,"bcs.dat",int);
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
    GPUinfo[3]=checkGpuMem();
    while(tw<t){
        if(tw>te && flag == 0){H2D(up ,2,nmp,DAT);flag++;}
        // get adaptative dt
        D2H(vp ,2,nmp,DAT);
        dt = CFL(vp_h,dx,dy,yd,tg,tw,nmp);
        // linear gravitational increase
        g  = getG(tw,tg);
        //              
        initD<<<grid,block>>>(mn_d,pn_d,fen_d,fin_d,pel_d,no,nel);
        basisD<<<Grid,Block>>>(N_d,dNx_d,dNy_d,p2n_d,p2e_d,e2n_d,xp_d,lp_d,xn_d,dx,dy,nmp,nn,ney,nel,no,xnmin,ynmin);

        accumD<<<Grid,Block>>>(mn_d,pn_d,fen_d,fin_d,N_d,dNx_d,dNy_d,sig_d,vp_d,mp_d,vol_d,p2n_d,g,nmp,nn,no);
        cudaDeviceSynchronize();

        solveD<<<grid,block>>>(fen_d,fin_d,mn_d,an_d,pn_d,vn_d,un_d,dt,no,bcs_d);
        cudaDeviceSynchronize();
        projectD<<<grid,block>>>(an_d,vn_d,N_d,vp_d,xp_d,up_d,p2n_d,dt,nmp,nn,no);
        cudaDeviceSynchronize();
        DMD<<<Grid,Block>>>(mn_d,un_d,N_d,vp_d,mp_d,p2n_d,dt,nmp,nn,bcs_d,no);
        cudaDeviceSynchronize();
        getdFD<<<grid,block>>>(un_d,dF_d,dNx_d,dNy_d,p2n_d,nmp,nn,no);
        cudaDeviceSynchronize();
        elastD<<<grid,block>>>(dF_d,eps_d,sig_d,ome_d,lp_d,vol_d,nmp,Del_d,dt);
        cudaDeviceSynchronize();
        if(tw>te){
            MCPlastD<<<grid,block>>>(sig_d,cohp_d,phip_d,epII_d,iDel_d,Hp,cohr,nmp);
            cudaDeviceSynchronize();
//             DPPlastD<<<grid,block>>>(sig_d,cohp_d,phip_d,epII_d,Hp,cohr,Kc,Gc,psi0,nmp);
//             cudaDeviceSynchronize();        
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
        if(FPS > 0 && DT>=((DAT)1.0/(DAT)FPS)){
            // save data
            if(tw>te){
                saveData++;
                D2HD(xp  ,2,nmp,"xp"  ,DAT,saveData,SIM);
                D2HD(lp  ,2,nmp,"lp"  ,DAT,saveData,SIM);
                D2HD(sig ,4,nmp,"sig" ,DAT,saveData,SIM);
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
    milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, startD, stopD);
    
    GPUinfo[0] = 1E-3*milliseconds;
    GPUinfo[1] = it/GPUinfo[0];
    GPUinfo[2] = (IO*it*sizeof(DAT)/(1024*1024*1024*GPUinfo[0]));
    printf("\n-----------------------------------");
    printf("\nGPU summary: MTPeff = %.2f [GB/s]",GPUinfo[2]);
    printf("\n-----------------------------------");
    printf("\n  time is %.2f s\n  after %d iterations \n  i.e., %.2f it/s\n",GPUinfo[0],it,GPUinfo[1]);
    // save data
    D2H(epII,1,nmp,DAT);
    D2H(xp ,2,nmp,DAT);
    D2H(up ,2,nmp,DAT);
    D2H(lp ,2,nmp,DAT);
    D2H(sig,4,nmp,DAT);
    save(epII,nmp  ,1,"CUDA_epII.dat");
    save(xp  ,2*nmp,1,"CUDA_xp.dat"  );
    save(lp  ,2*nmp,1,"CUDA_lp.dat"  );
    save(up  ,2*nmp,1,"CUDA_up.dat"  );
    save(sig ,4*nmp,1,"CUDA_sig.dat" );

    FILE* fidw1=fopen("GPUinfo.dat", "wb");
    fwrite(GPUinfo, 2*sizeof(DAT), 1, fidw1);
    fclose(fidw1);
    FILE *fidw2=fopen("numSaved.dat","wb");
    fwrite(&saveData,1,sizeof(saveData),fidw2);
    fclose(fidw2);

    cudaDeviceReset();
}




