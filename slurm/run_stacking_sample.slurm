#!/bin/bash
#SBATCH -J cwa_bayes_runs
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p skx-dev
#SBATCH -o logs/cwa_bayes_runs_test.o
#SBATCH -e logs/cwa_byaes_runs_test.e
#SBATCH -t 2:00:00
#SBATCH -A tagmap
#SBATCH --mail-user=Christopher.R.Peterson@utexas.edu
#SBATCH --mail-type=fail
#------------------------------------------------------

# N is number of nodes you want to use
# n is number of total jobs; should be N * 12
# p is the queue; use skx-dev for small jobs/testing and skx-normal for larger ones

module load launcher
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_hurdle.job
$LAUNCHER_DIR/paramrun