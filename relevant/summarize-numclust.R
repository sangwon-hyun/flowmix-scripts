## Synopsis: simple script to summarize all 80-89 simulation jobs.

## Load.
source("blockcv-helpers.R")
parse_args(args = commandArgs(trailingOnly = TRUE), verbose=TRUE)

## Setup
library(devtools)
library(parallel)
codedir = "/scratch/sangwonh/repos/flowcy/flowcy"
outputdir = "/scratch/sangwonh/output"
datadir = "/scratch/sangwonh/data"
scriptdir = "/scratch/sangwonh/scripts"
source(file.path(scriptdir,"blockcv-helpers.R"))
load_all(codedir)

## ## Figure out how many cores
## mynodename = strsplit(system("scontrol show hostname $SLURM_NODELIST | paste -d, -s", intern=TRUE),
##                            ",")[[1]]
## mycommand = "sinfo -o '%40N %c' --Node --long | grep "
## mystring = system(paste0(mycommand, mynodename), intern=TRUE)
## mc.cores = as.numeric(strsplit(mystring, "\\s+")[[1]][2])
## print("Number of cores:")
## print(mc.cores)
## ## mc.cores = 8
## ## Also try the shell command: nproc --all
## ## mc.cores = 24


## if(is.na(mc.cores)){
  mynodename = strsplit(system("scontrol show hostname $SLURM_NODELIST | paste -d, -s", intern=TRUE),
                             ",")[[1]]
  mycommand = "sinfo -o '%40N %c' --Node | grep "
  ## mycommand = "sinfo -o '%40N %c' --Node --long | grep "
  mystring = system(paste0(mycommand, mynodename), intern=TRUE)
  mc.cores = as.numeric(strsplit(mystring, "\\s+")[[1]][2])
## }

print("Number of cores:")
print(mc.cores)


## Measure the CV
blockcv_summary_sim2(nsim = nsim,
                     ## Simulation setup
                     blocktype = 2,
                     datatype = 9,
                     numclust = numclust,
                     nfold = 5,
                     nrep = 10,
                     cv_gridsize = 7,
                     ## Other settings.
                     outputdir = outputdir,
                     datadir = datadir,
                     mc.cores = mc.cores)
