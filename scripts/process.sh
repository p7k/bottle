#!/bin/bash
#SBATCH -A p30041
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8GB
#SBATCH -t 00:10:00
#SBATCH --job-name="process"
#SBATCH -o ../logs/out/process_expansions
#SBATCH -e ../logs/error/e_tmp_1
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu
ulimit -c 0
module load python/anaconda3.6
source activate casp
python -u process_all.py