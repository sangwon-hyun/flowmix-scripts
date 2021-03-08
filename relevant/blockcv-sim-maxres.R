## Synopsis: testing out the blocked CV approach
## codedir = "/home/rcf-proj3/sh7/repos/flowcy/flowcy"
codedir = "/staging/sh7/sangwonh/repos/flowcy/flowcy"
## outputdir = "/home/rcf-proj3/sh7/output"
outputdir = "/staging/sh7/sangwonh/output"
## datadir = "/home/rcf-proj3/sh7/data"
datadir = "/staging/sh7/sangwonh/data"
source("blockcv-helpers.R")
load_all(codedir)

###########
## Setup ##
###########
maxdev = NULL ## 1.5 ## Slightly larger than the truth, which is about (2.66-(-0.84)) * 0.3 = 1.05
append = NULL
blocksize = 20 ## Size for hourlong blocks (block type 2).
fitonly = FALSE ## Only do refitting
random_order = FALSE ## Make iilist into random order
maxres_once = FALSE
test = FALSE
manual_index = FALSE ## Also input ialpha and ibeta for this to work.
sigma_fac = 1

## Things related to the additive noise simulations
sim = FALSE ## In normal times, this is set to false, and only one round is run.
nsim = 1

############################################
### Read in some command line arguments ####
############################################
parse_args(args = commandArgs(trailingOnly = TRUE), verbose=TRUE)
stopifnot(!is.na(blocktype) & !is.na(datatype))
master_node = (arraynum == 1)
if(test) append="test"

## Test code:
## arraynum_max=10; blocktype=1;  datatype=6;  numclust=10; cv_gridsize=5; nfold=5;nrep=5
## arraynum_max=10; blocktype=1;  datatype=100;  numclust=4; cv_gridsize=5; nfold=5;nrep=5

################################################
### See how many cores are available to use ####
################################################
mynodename = strsplit(system("scontrol show hostname $SLURM_NODELIST | paste -d, -s", intern=TRUE),
                           ",")[[1]]
mycommand = "sinfo -o '%40N %c' --Node --long | grep "
mystring = system(paste0(mycommand, mynodename), intern=TRUE)
mc.cores = as.numeric(strsplit(mystring, "\\s+")[[1]][2])
print("Number of cores:")
print(mc.cores)

#############################################
#### Create directory to save results in ####
#############################################
destin = build_destin(outputdir, blocktype, datatype, numclust, append)
if(!dir.exists(destin)){
  dir.create(destin, recursive = TRUE)
  cat("Creating destin: ", destin, fill=TRUE)
} else {
  cat("All output goes out to destin: ", destin, fill = TRUE)
}

#######################################
## Load data once, from isim==1  ######
#######################################
isim = 1
datobj = loadRDS(file = file.path(destin,
                                  paste0("data-", datatype, "-", isim, ".Rds")))

################################
## Calculate maximum lambda ####
################################
if(maxres_once){
  print("maxres starting")
  maxres = get_max_lambda(destin, "maxres.Rdata",
                          ylist = ylist,
                          countslist = countslist,
                          X = X,
                          numclust = numclust,
                          maxdev = maxdev,
                          ## max_lambda_alpha = 10,
                          ## max_lambda_beta =  10)
                          max_lambda_alpha = 1,
                          max_lambda_beta =  3)
  print("maxres ending")
  q()
  if(fitonly)q()
} else {
  stop("maxres only!")
}
