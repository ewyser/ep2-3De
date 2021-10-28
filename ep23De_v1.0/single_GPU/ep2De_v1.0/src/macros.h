// header macros.h
#define PI (DAT)3.1415926535897931
#define alpha (DAT)0.0
#define IO (DAT)(26*no+nmp*(54+41*nn+2*4))
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define zeros(A,nx,ny,type)            type *A##_h,*A##_d                                                                 ;\
                                       A##_h = (type*)malloc(nx*ny*sizeof(type))                                          ;\
                                       cudaMalloc(&A##_d,(nx)*(ny)*sizeof(type))                                          ;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define load(A,nx,ny,Aname,type)       type *A##_h,*A##_d                                                                 ;\
                                       A##_h = (type*)malloc(nx*ny*sizeof(type))                                          ;\
                                       FILE* A##fid=fopen(Aname, "rb")                                                    ;\
                                       fread(A##_h, sizeof(type), nx*ny, A##fid)                                          ;\
                                       fclose(A##fid)                                                                     ;\
                                       cudaMalloc(&A##_d,nx*ny*sizeof(type))                                              ;\
                                       cudaMemcpy(A##_d,A##_h,nx*ny*sizeof(type),cudaMemcpyHostToDevice)                  ;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define save(A,nx,ny,Aname)            FILE* A##fidw=fopen(Aname, "wb")                                                   ;\
                                       fwrite(A##_h, sizeof(DAT), ((nx)*(ny)), A##fidw)                                   ;\
                                       fclose(A##fidw)                                                                    ;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define H2D(A,nx,ny,type)              cudaMemcpy(A##_d,A##_h,nx*ny*sizeof(type),cudaMemcpyHostToDevice)                  ;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define D2H(A,nx,ny,type)              cudaMemcpy(A##_h,A##_d,nx*ny*sizeof(type),cudaMemcpyDeviceToHost)                  ; 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
#define free_all(A)              free(A##_h); cudaFree(A##_d)                                                         ;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
#define D2HD(A,nx,ny,name,type,it,sim) char str##A[100],sim##A[100],fname##A[100] = name                                  ;\
                                       sprintf(str##A,"_%d",it)                                                           ;\
                                       sprintf(sim##A,"_%d",sim)                                                          ;\
                                       strcat(fname##A,str##A)                                                            ;\
                                       strcat(fname##A,sim##A)                                                            ;\
                                       strcat(fname##A,".dat")                                                            ;\
                                       D2H (A,nx,ny,type)                                                                 ;\
                                       save(A,nx,ny,fname##A)                                                             ;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////                                            