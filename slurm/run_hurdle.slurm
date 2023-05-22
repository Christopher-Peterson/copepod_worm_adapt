#!/bin/bash
#SBATCH -J cwa_bayes_runs
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p skx-dev
#SBATCH -o logs/cwa_bayes_runs_test.o
#SBATCH -e logs/cwa_byaes_runs_test.e
#SBATCH -t 2:00:00
#------------------------------------------------------

# N is number of nodes you want to use
# n is number of total jobs; should be N * 12
# p is the queue; use skx-dev for small jobs/testing and skx-normal for larger ones

module load launcher
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=jobs/run_hurdle.job

# Create the job file & model specifications
docker_stan R/prep_runs.r $LAUNCHER_JOB_FILE
# Run them
$LAUNCHER_DIR/paramrun