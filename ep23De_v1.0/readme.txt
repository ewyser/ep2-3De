______________________________________________________________________________________
--------------------------------------------------------------------------------------
QUICK START: example for Experience 1a
--------------------------------------------------------------------------------------
	1) In ep23DSuite, select the 3D solver folder (ep3De_v1.0) and navigate in:
		.../ep3De_v1.0/scripts/setting_Exp1a/
	2) Open the MATLAB routine compileExp1a.m and:
		a) select the actual GPU architecture installed on the computer, i.e., 
		--gpu_architecture=sm_xx,
		b) define the number of threads per block BCK1 and BCK2,
		c) set the precision arithmetic typeD = 'float' or typeD = 'double',
		d) (opt.) define FPS (frame per seconds) and D (local damping).
	2) Open the MATLAB routine Exp1a.m, run.
	3) Post-process the results using the MATLAB routine postprocessing.m 
______________________________________________________________________________________
--------------------------------------------------------------------------------------
GPU ARCHITECTURE: few examples (Maxwell to Ampere architectures) 
--------------------------------------------------------------------------------------
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
______________________________________________________________________________________