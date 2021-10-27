#include "stdio.h"
#include "cuda.h"
#include "time.h"
#include "mpi_macros.h"
#include "cpu_functions.h"
#include "gpu_kernels.h"
#include "mpi_gpu.h"
#include "cuda_profiler_api.h"
int main(int argc, char *argv[]){
    // Set up GPUs
    set_up_parallelisation();
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
    load_mpi(param,13,1,"param_rank",me,DAT)
    int nmp  = (int)param_h[0];
    int no  = (int)param_h[1];
    int nel = (int)param_h[2];
    int nex = (int)param_h[3];
    int ney = (int)param_h[4];
    int nez = (int)param_h[5];
    int nn  = (int)param_h[6];
    DAT dx = param_h[7];
    DAT dy = param_h[8];
    DAT dz = param_h[9];
    DAT xnmin = param_h[10];
    DAT ynmin = param_h[11];
    DAT znmin = param_h[12];
    int nnx = nex+1;
    int nny = ney+1; 
    int nnz = nez+1;
    //if(me==0){
        printf("\n np=%d, nel=%d, no=%d, nn=%d, xnmin = %f, ynmin = %f, dx = %f, dy = %f, ney = %d for rank %d \n",nmp,nel,no,nn,xnmin,ynmin,dx,dy,ney,me);
    //}
    load(Del ,36        ,1,"Del.dat" ,DAT); 
    load_mpi(e2n,nel,nn,"e2n_rank",me,int)
    load_mpi(mp,nmp,1,"mp_rank",me,DAT)
    load_mpi(lp,nmp,3,"lp_rank",me,DAT)
    load_mpi(xp,nmp,3,"xp_rank",me,DAT)
    load_mpi(xn,no,3,"xn_rank",me,DAT)
    load_mpi(bcs,no,3,"bcs_rank",me,int)
    load_mpi(vol,nmp,1,"vol_rank",me,DAT)
    load_mpi(cohp,nmp,1,"cohp_rank",me,DAT)
    load_mpi(phip,nmp,1,"phip_rank",me,DAT)
    // initialize arrays
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
    // MPI send & recv buffers
    init_sides(mn ,1*(BFS)*nnz*nnx);
    init_sides(fen ,3*(BFS*nnz*nnx));
    init_sides(fin ,3*(BFS*nnz*nnx));
    init_sides(pn ,3*(BFS*nnz*nnx));
    init_sides(un ,3*(BFS*nnz*nnx));
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
    // CUDA
    dim3 Grid, Block;
    dim3 grid, block;
    Block.x = BCK1;
    Block.y = nn;
    Grid.x  = ceil(max(nmp,no)/Block.x)+1;
    Grid.y  = 1;
    block.x = BCK2;
    block.y = 1;
    grid.x  = ceil(max(nmp,no)/block.x)+1;
    grid.y  = 1;
    dim3 grid_mpi0, block_mpi0;
    block_mpi0.x = 1; grid_mpi0.x = 3;      
    block_mpi0.y = 512; grid_mpi0.y =ceil(BFS*(nnz*nnx)/block_mpi0.y)+1; 

    // Solver
    int flag = 0;
    DAT DT   = 0.0;
    if(me==0){
        printf("\nGPU set up:");
        printf("\n-----------------------------------");
        printf("\n  Grid = %d, Block = %d",Grid.x*Grid.y,Block.x*Block.y);
        printf("\n  grid = %d, block = %d",grid.x*grid.y,block.x*block.y);        
    }
    cudaEventCreate(&startD);
    cudaEventCreate(&stopD);
    cudaEventRecord(startD);




    // Check Mem Usage
    GPUinfo[3]=checkGpuMem(me)-GPUinfo[3];
    checkCudaErrors();
    dt = (DAT)0.5*(DAT)dx/(DAT)yd;
    if(me==0){
        printf("\n-----------------------------------");
        printf("\nWorkload: Profiling");
        printf("\n-----------------------------------");
    }
    while(tw<t){
                if (it==90 && (me==0)){ cudaProfilerStart(); }
        if (it==100 && (me==0)){ cudaProfilerStop();  }
        if(me==0 & it%100==0){
        printf("\n t/T = %.2f %%",tw/t*100);    
        }
//         if (it==90 && (me==0)){ cudaProfilerStart(); }
//         if (it==180 && (me==0)){ cudaProfilerStop();  }
        MPI_Barrier(topo_comm);
        // get adaptative dt
        
        // linear gravitational increase
        g  = getG(tw,tg);
        //              
        initD<<<grid,block>>>(mn_d,pn_d,fen_d,fin_d,pel_d,no,nel);
        cudaDeviceSynchronize();        
        basisD<<<Grid,Block>>>(N_d,dNx_d,dNy_d,dNz_d,p2n_d,p2e_d,e2n_d,xp_d,lp_d,xn_d,dx,dy,dz,nmp,nn,nex,ney,nez,nel,no,xnmin,ynmin,znmin);
        cudaDeviceSynchronize();
        accumD<<<Grid,Block>>>(mn_d,pn_d,fen_d,fin_d,N_d,dNx_d,dNy_d,dNz_d,sig_d,vp_d,mp_d,vol_d,p2n_d,g,nmp,nn,no);
        cudaDeviceSynchronize();
            MPI_Barrier(topo_comm);
            update_sides(mn,BFS*nnx*nnz,1);
            cudaDeviceSynchronize();
            update_sides(fen,BFS*nnx*nnz,3);
            cudaDeviceSynchronize();
            update_sides(fin,BFS*nnx*nnz,3);
            cudaDeviceSynchronize();
            update_sides(pn,BFS*nnx*nnz,3);
            cudaDeviceSynchronize();
            MPI_Barrier(topo_comm);
        solveD<<<grid,block>>>(fen_d,fin_d,mn_d,an_d,pn_d,vn_d,un_d,dt,no,bcs_d);
        cudaDeviceSynchronize();
        projectD<<<grid,block>>>(an_d,vn_d,N_d,vp_d,xp_d,up_d,p2n_d,dt,nmp,nn,no);
        cudaDeviceSynchronize();
        DMD<<<Grid,Block>>>(mn_d,un_d,N_d,vp_d,mp_d,p2n_d,dt,nmp,nn,bcs_d,no);
        cudaDeviceSynchronize();
// // 
// //             DMD1<<<Grid,Block>>>(un_d,N_d,vp_d,mp_d,p2n_d,nmp,nn,no);
// //             cudaDeviceSynchronize();
// //             DMD2<<<grid,block>>>(mn_d,un_d,dt,bcs_d,no);
// //             cudaDeviceSynchronize();
            MPI_Barrier(topo_comm);
            update_sides(un,BFS*nnx*nnz,3);
            cudaDeviceSynchronize();
            MPI_Barrier(topo_comm);
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
     }
    cudaEventRecord(stopD);
    cudaEventSynchronize(stopD);
    milliseconds = (DAT)0.0;
    cudaEventElapsedTime(&milliseconds, startD, stopD);
    
    GPUinfo[0] = (DAT)1E-3*milliseconds;
    GPUinfo[1] = it/GPUinfo[0];
    GPUinfo[2] = (IO*it*sizeof(DAT)/(1024*1024*1024*GPUinfo[0]));
    if(me==0){
        printf("\n-----------------------------------");
        printf("\nGPU summary: MTPeff = %.2f [GB/s]",GPUinfo[2]);
        printf("\n-----------------------------------");
        printf("\n  time is %.2f s\n  after %d iterations \n  i.e., %.2f it/s\n",GPUinfo[0],it,GPUinfo[1]);
    }
    // save data
save_mpi(mn,no,1,"mn_rank",me,DAT)
save_mpi(fin,no,3,"un_rank",me,DAT)
save_mpi(up,nmp,3,"up_rank",me,DAT)
save_mpi(xp,nmp,3,"xp_rank",me,DAT)
save_mpi(sig,nmp,6,"sig_rank",me,DAT)
save_mpi(epII,nmp,1,"epII_rank",me,DAT)


getname("GPUinfo",me)
    FILE* fidw=fopen(filename, "wb");
    fwrite(GPUinfo, 4*sizeof(DAT), 1, fidw);
    fclose(fidw);



    MPI_Finalize();
    cudaDeviceReset();
}