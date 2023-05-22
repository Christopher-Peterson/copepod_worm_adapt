#!/bin/bash
#SBATCH -J cwa_stack
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -p skx-normal
#SBATCH -o logs/cwa_stack.o
#SBATCH -e logs/cwa_stack.e
#SBATCH -t 2:00:00
#------------------------------------------------------

# N is number of nodes you want to use
# n is number of total jobs; should be N * 12
# p is the queue; use skx-dev for small jobs/testing and skx-normal for larger ones

module load launcher
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=jobs/model_stacking

# Create the job file & model specifications
# arguments; feel free to adjust
seed=1234
cores=48
reps=60
chunks=4

docker_stan R/weights_1_pbma.r $seed $LAUNCHER_JOB_FILE $cores $reps $chunks

# Divide the data into 4 chunks & run stacking on each
$LAUNCHER_DIR/paramrun

wait
# Combine the four chunks, filter out the models that non-zero performance, and re-run stacking on them.
docker_stan R/weights_3_stacking_calculate.r
