#!/bin/bash
#
# USAGE: /full/path/to/./q_MC.out [--help]
#                                 [--init dimension [outfile]]
#                                 [--nonstop dimension outfile steps [increment]]
#                                 [infile [outfile] steps [increment]]
#
#SBATCH --job-name=Ti64Datafit
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhatts8@rpi.edu
#SBATCH --partition small 
#SBATCH -t 03:00:00
#SBATCH -N 32 
#SBATCH -n 2048
#SBATCH --overcommit
#SBATCH -D /gpfs/u/home/ACME/ACMEsgrb/scratch/Ti6Al4V_betaPhase/16
#SBATCH -o /gpfs/u/home/ACME/ACMEsgrb/scratch/Ti6Al4V_betaPhase/16/Ti.log


srun --runjob-opts="--mapping TEDCBA" /gpfs/u/home/ACME/ACMEsgrb/scratch/Ti6Al4V_betaPhase/q_MC.out --nonstop 2 Ti.000000.dat 250000 5000 0 4 1361 1361 > result.txt
