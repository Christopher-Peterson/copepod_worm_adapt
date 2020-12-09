
<!-- README.md is generated from README.Rmd. Please edit that file -->

# copepod\_worm\_adapt

Scripts for the Bayesian hurdle models in this manuscript (add details
on submission).

## Instructions

### TACC setup

Part of this uses TACC, though it could probably be run on a different
HPC system. This is configured for **Stampede2**.

1.  Create a directory in `$SCRATCH` for the project, then copy the R &
    data folders to it. Make a folder called “logs” with `mkdir`.
2.  Download the docker image.  
    a. Make the folder `$WORK/singularity` and go there.  
    b. Open an `idev` session c. Run `ml tacc-singularity` d. Run
    `singularity pull docker://crpeters/docker-stan:21.1.2` e. This will
    download a docker image that has all the analysis packages in it.
    f. Run `ls *.sif` and copy the name of the sif file.
3.  Configure docker\_stan shell script a. Copy the docker\_stan file
    into your local bin folder (could be `$HOME/bin` or `$WORK/apps/bin`
    or something like that). b. Edit the file with nano or something;
    replace the sif file it references with the one you coppied in the
    previous step, then save. c. Make the file executable.
4.  In an idev session, go to the project directory and run
    `docker-stan R/prep_runs.r`. Make sure that the .job and .rds file
    are created.
5.  Create a SLURM script to run the .job file and submit it. The
    sample.slurm file is a starting place.
