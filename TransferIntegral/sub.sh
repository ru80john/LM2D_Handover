#!/bin/bash
#SBATCH -J TI
#SBATCH -o ./TI.stdout%j
#SBATCH -e ./TI.stderr%j
##SBATCH --mail-type=end
##SBATCH --mail-type=fail
#SBATCH --mem=256000
#SBATCH --ntasks-per-node=48
#SBATCH --reservation=PIYAN


python calc_TI_parallel.py
