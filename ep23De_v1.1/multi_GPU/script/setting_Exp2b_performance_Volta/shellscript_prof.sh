#!/bin/bash

#------------------------------------------------------------------
# load what should be loaded
module load cuda/11.0
module load openmpi/215_g65_c10
#------------------------------------------------------------------
# remove what should be removed
rm *.out
#------------------------------------------------------------------
mpirun=$(which mpirun)
nvprof=$(which nvprof)
nprocs=$7
#------------------------------------------------------------------
#compile GPU code
#nvcc --compiler-bindir mpic++ --compiler-options -O3 ../../src_mpi/gpu_main.cu -o gpu.out -arch=sm_70 -O3 -DSIM=$1 -DDAT=$2 -DFPS=$3 -DD=$4 -DBCK1=$5 -DBCK2=$6 -DDIMS_X=$7 -DBFS=$8
#run_cmd="-np $nprocs -rf gpu_rankfile_volta --mca btl_openib_if_exclude mlx4_0,mlx4_1 --mca btl_openib_ignore_locality 1"
nvcc ../../src_mpi/gpu_main_prof.cu -o gpu.out -arch=sm_52 --compiler-options -O3 --compiler-bindir mpic++ -DSIM=$1 -DDAT=$2 -DFPS=$3 -DD=$4 -DBCK1=$5 -DBCK2=$6 -DDIMS_X=$7 -DBFS=$8
run_cmd="-np $nprocs -rf gpu_rankfile_titan --mca btl_openib_if_include mlx4_0,mlx4_1 --mca btl_openib_ignore_locality 1"

echo $mpirun $run_cmd
# mpi run
$mpirun $run_cmd ./gpu.out
# nvidia profiler
$mpirun $run_cmd $nvprof --profile-from-start off --export-profile mpm_3D.%q{OMPI_COMM_WORLD_RANK}.prof -f ./gpu.out