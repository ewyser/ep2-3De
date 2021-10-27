#include "mpi.h"

#define MPI_DAT  MPI_REAL

#define NDIMS  1
#define zeros_d(A,nxy)      DAT *A##_d; cudaMalloc(&A##_d,(nxy)*sizeof(DAT));
#define NREQS               (2*2*NDIMS)
#define neighbours(dim,nr)  __neighbours[nr + dim*2]
int dims[3]={DIMS_X,0,0};
int coords[3]={0,0,0};
int* coords_d=NULL;
int nprocs=-1, me=-1, me_loc=-1, gpu_id=-1;
int __neighbours[2*NDIMS]={-1}; // DEBUG neighbours(DIM.nr) macro in my case oposite to Sam's
int reqnr=0, tag=0;
int periods[NDIMS]={0};
int reorder=1;
MPI_Comm    topo_comm=MPI_COMM_NULL;
MPI_Request req[NREQS]={MPI_REQUEST_NULL};

// // // // // // #define set_up_grid()  dim3 grid, block, grid_mpi0, block_mpi0, grid_mpi1, block_mpi1; \
// // // // // //     block.x      = BLOCK_X; grid.x      = GRID_X; \
// // // // // //     block.y      = BLOCK_Y; grid.y      = GRID_Y; \
// // // // // //     block_mpi0.x = 1;       grid_mpi0.x = 1;      \
// // // // // //     block_mpi0.y = BLOCK_Y; grid_mpi0.y = GRID_Y; \
// // // // // //     block_mpi1.x = BLOCK_X; grid_mpi1.x = GRID_X; \
// // // // // //     block_mpi1.y = 1;       grid_mpi1.y = 1;

void __set_up_parallelisation(int argc, char *argv[]){
    // GPU STUFF
    cudaSetDeviceFlags(cudaDeviceMapHost); // DEBUG: needs to be set before context creation !
    const char* me_str     = getenv("OMPI_COMM_WORLD_RANK");
    const char* me_loc_str = getenv("OMPI_COMM_WORLD_LOCAL_RANK");
    me     = atoi(me_str);
    me_loc = atoi(me_loc_str);
    gpu_id = me_loc;
    cudaSetDevice(gpu_id); cudaGetDevice(&gpu_id);
    cudaDeviceReset(); cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);  // set L1 to prefered
    // MPI STUFF
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Dims_create(nprocs, NDIMS, dims);
    MPI_Cart_create(MPI_COMM_WORLD, NDIMS, dims, periods, reorder, &topo_comm);
    MPI_Comm_rank(topo_comm, &me);
    MPI_Cart_coords(topo_comm, me, NDIMS, coords);
    cudaMalloc(&coords_d,3*sizeof(int)); cudaMemcpy(coords_d ,coords,3*sizeof(int),cudaMemcpyHostToDevice);
    for (int i=0; i<NDIMS; i++){ MPI_Cart_shift(topo_comm, i, 1, &(neighbours(i,0)), &(neighbours(i,1))); }
    if (me==0){ printf("nprocs=%d,dims(1)=%d,dims(2)=%d \n", nprocs,dims[0],dims[1]); 
    printf("gpu id = %d with lett and right folks %d and %d \n", me,neighbours(0,0),neighbours(0,1)); 
    }
}
#define set_up_parallelisation()  __set_up_parallelisation(argc, argv);
// MPI buffer init
// SEND
__global__ void write_to_mpi_sendbuffer_00(DAT* A_send_00,DAT* A,  const int no_A,  const int no, const int nny, int ncomp){
    int iy = blockIdx.y*blockDim.y + threadIdx.y; // thread ID, dimension y
    int ix = blockIdx.x;    
    if (iy<(no_A) && ix<ncomp){
        A_send_00[iy+ix*no_A] = A[iy+ix*no];
    }
} // left buffer
__global__ void write_to_mpi_sendbuffer_01(DAT* A_send_01,DAT* A, const int no_A, const int no, const int nny, int ncomp){
    int iy = blockIdx.y*blockDim.y + threadIdx.y; // thread ID, dimension y
    int ix = blockIdx.x;
    if (iy<(no_A) && ix<ncomp){
        A_send_01[no_A-1-iy+ix*no_A] = A[(no-1-iy)+ix*no]; 
    }
} // right buffer
// RECEIVE
__global__ void read_from_mpi_recvbuffer_00(DAT* A,DAT* A_recv_00, const int no_A, const int no, const int nny, int ncomp){
    int iy = blockIdx.y*blockDim.y + threadIdx.y; // thread ID, dimension y
    int ix = blockIdx.x;
    if (iy<(no_A) && ix<ncomp){
        A[iy+ix*no] = A[iy+ix*no]+A_recv_00[iy+ix*no_A]; 
    }
} // left buffer
__global__ void read_from_mpi_recvbuffer_01(DAT* A,DAT* A_recv_01, const int no_A, const int no, const int nny, int ncomp){
    int iy = blockIdx.y*blockDim.y + threadIdx.y; // thread ID, dimension y
    int ix = blockIdx.x;
    if (iy<(no_A) && ix<ncomp){
        A[((no-1)-iy)+ix*no] = A[((no-1)-iy)+ix*no] + A_recv_01[(no_A-1-iy)+ix*no_A]; 
    }
} // right buffer
#define update_sides(A,no_A,ncomp) cudaDeviceSynchronize(); \
                                  if (neighbours(0,1) != MPI_PROC_NULL)    write_to_mpi_sendbuffer_01<<<grid_mpi0,block_mpi0>>>(A##_send_01_d, A##_d, no_A, no, nnx*nnz,ncomp);cudaDeviceSynchronize(); \
                                  if (neighbours(0,0) != MPI_PROC_NULL){  MPI_Irecv(A##_recv_00_d, (no_A*ncomp), MPI_DAT, neighbours(0,0), tag, topo_comm, &(req[reqnr]));  reqnr++;  } \
                                  if (neighbours(0,1) != MPI_PROC_NULL){  MPI_Isend(A##_send_01_d, (no_A*ncomp), MPI_DAT, neighbours(0,1), tag, topo_comm, &(req[reqnr]));  reqnr++;  } \
                                   cudaDeviceSynchronize(); \
                                  if (neighbours(0,0) != MPI_PROC_NULL)    write_to_mpi_sendbuffer_00<<<grid_mpi0,block_mpi0>>>(A##_send_00_d, A##_d, no_A, no, nnx*nnz,ncomp); cudaDeviceSynchronize(); \
                                  if (neighbours(0,1) != MPI_PROC_NULL){  MPI_Irecv(A##_recv_01_d, (no_A*ncomp), MPI_DAT, neighbours(0,1), tag, topo_comm, &(req[reqnr]));  reqnr++;  } \
                                  if (neighbours(0,0) != MPI_PROC_NULL){  MPI_Isend(A##_send_00_d, (no_A*ncomp), MPI_DAT, neighbours(0,0), tag, topo_comm, &(req[reqnr]));  reqnr++;  } \
                                   MPI_Waitall(reqnr,req,MPI_STATUSES_IGNORE);  reqnr=0;  for (int j=0; j<NREQS; j++){ req[j]=MPI_REQUEST_NULL; }; \
                                   cudaDeviceSynchronize(); \
                                  if (neighbours(0,0) != MPI_PROC_NULL)   read_from_mpi_recvbuffer_00<<<grid_mpi0,block_mpi0>>>(A##_d, A##_recv_00_d, no_A, no, nnx*nnz,ncomp); cudaDeviceSynchronize(); \
                                  if (neighbours(0,1) != MPI_PROC_NULL)   read_from_mpi_recvbuffer_01<<<grid_mpi0,block_mpi0>>>(A##_d, A##_recv_01_d, no_A, no, nnx*nnz,ncomp); cudaDeviceSynchronize(); \
                                  cudaDeviceSynchronize(); \

#define init_sides(A,no_A)   zeros_d(A##_send_00 ,no_A); zeros_d(A##_send_01 ,no_A); zeros_d(A##_recv_00 ,no_A); zeros_d(A##_recv_01 ,no_A);
#define free_sides(A)             cudaFree(A##_send_00_d); cudaFree(A##_send_01_d); cudaFree(A##_send_10_d); cudaFree(A##_send_11_d); cudaFree(A##_recv_00_d); cudaFree(A##_recv_01_d); cudaFree(A##_recv_10_d); cudaFree(A##_recv_11_d);
