#!/bin/bash


#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mail-user=yfs2333333@gmail.com  # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err


module load gcc/10.2.0-zuvaafu
module load r/4.0.4-py3-4khjixy
module load cmake
Rscript Three_arm.R --r $1 > NMA_$SLURM_JOBID.Rout


