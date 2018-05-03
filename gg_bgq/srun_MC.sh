#!/bin/bash
#
# USAGE: /full/path/to/./q_MC.out [--help]
#                                 [--init dimension [outfile]]
#                                 [--nonstop dimension outfile steps [increment]]
#                                 [infile [outfile] steps [increment]]
#
#SBATCH --job-name=PPC_Project
#SBATCH --mail-type=END
#SBATCH --mail-user=peters9@rpi.edu
#SBATCH --partition medium
#SBATCH -t 00:05:00
#SBATCH -N 256
#SBATCH -n 16384
#SBATCH --overcommit
#SBATCH -D /gpfs/u/home/PCP7/PCP7ptrs/gg2
#SBATCH -o /gpfs/u/home/PCP7/PCP7ptrs/gg2/Proj.log


srun --runjob-opts="--mapping TEDCBA" /gpfs/u/home/PCP7/PCP7ptrs/gg2/q_MC.out bgstart 1000
