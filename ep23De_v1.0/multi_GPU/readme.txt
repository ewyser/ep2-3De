______________________________________________________________
--------------------------------------------------------------
DISCLAIMER
--------------------------------------------------------------
The multi-GPU implementation uses the openMPI API and, it was only tested on a CentOS system. We do not guarantee the current implementation is stable on any other computing architectures. However, we believe it can still provides a general overview of a multi-GPU implementation of the material point method. 
______________________________________________________________
--------------------------------------------------------------
GPU ARCHITECTURE: few examples (Maxwell to Ampere architectures) 
--------------------------------------------------------------
We provide only few architecture flags, i.e.,
	gpu-architecture=sm_xx:
	xx=52; Maxwell architecture, Geforce 900 series
	xx=60; Pascal architecture , Geforce 10 series
				     Tesla P100
				     Titan Xp
	xx=70; Volta arhitecture   , Tesla V100 
				     Titan V
	xx=75; Turing architecture , Geforce 20 series
				     Geforce 16 series
	xx=80; Ampere architecture , Geforce 30 series
				     A100
For other architectures, check on the Nvidia website http://Nvidia.com and/or
https://patrickorcl.medium.com/compile-with-nvcc-3566fbdfdbf for a more exhaustive 
listing of architecture flags.
______________________________________________________________