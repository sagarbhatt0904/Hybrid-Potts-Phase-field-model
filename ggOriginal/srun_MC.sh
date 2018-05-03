#!/bin/bash
#
# USAGE: /full/path/to/./q_MC.out [--help]
#                                 [--init dimension [outfile]]
#                                 [--nonstop dimension outfile steps [increment]]
#                                 [infile [outfile] steps [increment]]
#
#SBATCH --job-name=PPC_Project
#SBATCH --mail-type=END
#SBATCH --mail-user=bhatts8@rpi.edu
#SBATCH --partition debug 
#SBATCH -t 00:10:00
#SBATCH -N 1 
#SBATCH -n 64
#SBATCH --overcommit
#SBATCH -D /gpfs/u/home/PCP7/PCP7sgrb/scratch/gg/
#SBATCH -o /gpfs/u/home/PCP7/PCP7sgrb/scratch/gg/Proj.log


srun --runjob-opts="--mapping TEDCBA" /gpfs/u/home/PCP7/PCP7sgrb/scratch/gg/q_MC.out bgstart 1000
