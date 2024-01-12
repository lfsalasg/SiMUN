#!/bin/sh
#SBATCH --partition=LocalQ
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=maxTime
#SBATCH --job-name=filename
#SBATCH -o resultqiime_%N_%j.out
#SBATCH -e resultqiime_%N_%j.err

$RASPA_DIR/bin/simulate $1
