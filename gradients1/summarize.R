# Synopsis: simple script to summarize all 80-89 simulation jobs.

## Number of simmulations.
## nsim = 100

## Adit
source("blockcv-helpers.R")
parse_args(args = commandArgs(trailingOnly = TRUE), verbose=TRUE)

## Setup
outputdir = "~/stagedir/output"
datadir = "~/stagedir/data"
codedir = "~/stagedir/repos/flowcy/flowcy"
load_all(codedir)

## Figure out how many cores
mynodename = strsplit(system("scontrol show hostname $SLURM_NODELIST | paste -d, -s", intern=TRUE),
                           ",")[[1]]
mycommand = "sinfo -o '%40N %c' --Node --long | grep "
mystring = system(paste0(mycommand, mynodename), intern=TRUE)
mc.cores = as.numeric(strsplit(mystring, "\\s+")[[1]][2])
print("Number of cores:")
print(mc.cores)
## mc.cores = 8
## Also try the shell command: nproc --all
## mc.cores = 24

## Measure the errors and fits
## datatypes = c(82:83)
## for(datatype in datatypes){
## print(datatype)

## Sample setup
## datatype = 83
## nsim = 30


## Measure the CV
la(codedir)
blockcv_summary_sim(nsim = nsim,
                     ## Simulation setup
                     blocktype = 2, ##2,
                     datatype = datatype,
                     numclust = 2,
                     cv_gridsize = 7,
                     ## Other settings.
                     outputdir = outputdir,
                     datadir = datadir,
                     mc.cores = mc.cores)

## }
