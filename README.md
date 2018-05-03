# PPC Final Project - Sagar Bhatt, Scott Peters #
----------------------------------------------------------------
## Grain Growth simulation using Hybrid Potts-Phase field model.

The code is presented in two folders. The folder 'gg_bgq' contains modifications for BG/Q especially. Any execution of the code for BG/Q should be done from that folder only. The folder 'ggOriginal' contains the code in its entierity for all non-BG/Q systems. To change the grid size, open gg.cpp in either folder and change the value of `#define L`. 

### To compile and run for Blue Gene:
1) include all this code in a folder in BGQ. 
2) edit the `srun_MC.sh` script to correctly represent both your account name, and the location of this folder
	1) Change the email from "peters9@rpi.edu" to your preferred email if you like. Otherwise please remove this line!
	2) The folder my code was run from was at `.../PCP7/PCP7ptrs/gg2/`, change lines that have this to match your location
3) don't forget to type `module load xl` if you haven't done so already!
4) compile the code by typing `make bgqmc`
5) submit the code to be run by typing `sbatch srun_MC.sh`

### To Compile and run for non-BG/Q systems:
1) `cd ggOriginal/`
2) To make use: `make parallel`
3) To Run : 
	1) Create a grid: `mpirun -n <np> ./parallel --example 2 <file_name> `		
	2) Run from the previously created grid: `mpirun -n <np> ./parallel <file_name> <number_of_iteration> <step increment to write files>`
