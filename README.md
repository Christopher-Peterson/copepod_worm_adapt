
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Data and analyses for “Local adaptation and host specificity to copepod intermediate hosts by the *Schistocephalus solidus* tapeworm”

These include the data and scripts for the Bayesian hurdle model
ensemble analysis.

### Datasets

There are 3 datasets in here: “chapter2.copepods.cleaned.Aug2020.csv”,
“chapter2.copepods.cleaned.csv”, and “chapter_2\_copepod_for_bayes.csv.”
They are identical except that the Aug2020 one had a few extra
space/tabs for a few individuals in the numb.worm column added by
mistake during data collection. The “for_bayes” version of the dataset
was modified from the Aug2020 version by a script in the repository, and
is the one used in the actual analysis.

## Re-analysis instructions

### TACC setup

The hurdle models were originally run on the Texas Advanced Computing
Center’s (TACC) Stampede2 supercomputer. It could probably be adapted to
other HPC systems.

1.  Create a directory in `$SCRATCH` for the project, then copy the R &
    data folders to it. Make a folder called “logs” with `mkdir`.
2.  Download the docker image.  
    a. Make the folder `$WORK/singularity` and go there.  
    b. Open an `idev` session c. Run `ml tacc-singularity` d. Run
    `singularity pull docker://crpeters/docker-stan:21.1.2` e. This will
    download a docker image that has all the analysis packages in it.
    f. Run `ls *.sif` and copy the name of the sif file.
3.  Configure docker_stan shell script a. Copy the docker_stan file into
    your local bin folder (could be `$HOME/bin` or `$WORK/apps/bin` or
    something like that). b. Edit the file with nano or something;
    replace the sif file it references with the one you coppied in the
    previous step, then save. c. Make the file executable.

### Run the Hurdle models

This will create a list of all possible hurdle model candidates and
submit them to TACC .  

1.  Submit the hurdle model batch job with
    `sbatch slurm/run_hurdle.slurm`. You may need to edit the slurm
    file’s parameters (cores, nodes, etc) to suite your HPC system.

2.  If the all of the tasks aren’t completed, run
    `sbatch slurm/continue_hurdle.slurm` to continue the run. This
    should be done repeatedly until all models have been run

For each model, there should be two output files: the posterior
distribution and the loo log predictive density.

### Model Stacking

Model stacking is a three-step process, which is fully handled by the
`slurm/model_stacking.slurm` script. If it times out, the easiest option
is to give the job more time and then re-run it.

### Post-processing

Several other scripts are useful for analyzing the stacking posterior:

- `get_ev.r` provides expected values
- `get_effect_sizes.R` calculates effect sizes
- `make_effect_size_plots.r` visualizes the effect sizes
- `copepod_corrrelation_plot.r` creates figure S1 in the manuscript. All
  of these are probably best run on the HPC system used for the other
  analyses.
