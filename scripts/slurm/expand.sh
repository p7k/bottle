#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --mem=0
#SBATCH -t 4:00:00
#SBATCH --job-name="exp"
#SBATCH --output=/path/to/outlogs/%A
#SBATCH --error=/path/to/error/logs/%A
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu

# Args
scripts_dir=/home/spn1560/bottle/scripts
starters=bottle_targets_24
generations=2
processes=50

# Commands
ulimit -c 0
module purge
source /home/spn1560/.cache/pypoetry/virtualenvs/bottle-OcfYl-TY-py3.12/bin/activate
python $scripts_dir/run_pickaxe.py $starters $generations -p $processes