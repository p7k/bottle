#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --mem=0
#SBATCH -t 3:00:00
#SBATCH --job-name="map_iv3"
#SBATCH -o ../logs/out/map_iv3
#SBATCH -e ../logs/error/e_tmp
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu
ulimit -c 0
module load python/anaconda3.6
source activate mine
python -u map_operators.py ../data/rules/JN3604IMT_rules.tsv ../data/sprhea/sprhea_240310_v3.json ../artifacts/operator_mapping_sprhea_v3_imt_ops.tsv