## Synopsis: generate new simulated data (datatype=80:89, isim=1:nsim)

#######################
## General settings ###
#######################
blocktype = 2
numclust = 2
sigmalist = seq(from = 3 * 0.1, to = 3 * 0.9, length = 9)
nsim = 100

## codedir = "/home/rcf-proj3/sh7/repos/flowcy/flowcy"
codedir = "/staging/sh7/sangwonh/repos/flowcy/flowcy"
## outputdir = "/home/rcf-proj3/sh7/output"
outputdir = "/staging/sh7/sangwonh/output"
## datadir = "/home/rcf-proj3/sh7/data"
datadir = "/staging/sh7/sangwonh/data"
source("blockcv-helpers.R")
load_all(codedir)

################################################
## Obtain the number of cores in this machine ##
################################################
mynodename = strsplit(system("scontrol show hostname $SLURM_NODELIST | paste -d, -s", intern=TRUE),
                           ",")[[1]]
mycommand = "sinfo -o '%40N %c' --Node --long | grep "
mystring = system(paste0(mycommand, mynodename), intern=TRUE)
mc.cores = as.numeric(strsplit(mystring, "\\s+")[[1]][2])
print("Number of cores:")
print(mc.cores)

#########################################
## Generate some new simulated data. ####
#########################################
for(datatype in 80:89){
  print("Datatype")
  printprogress(datatype, 80:89, fill=TRUE)
  destin = file.path(outputdir,
                     paste0("blockcv-", blocktype, "-", datatype, "-", numclust))

  ## for(isim in 1:nsim){
  mclapply(1:nsim, function(isim){
    print("Simulation")
    printprogress(isim, nsim, fill=TRUE)

    ## Generate data once
    set.seed(datatype * 1000 + isim)
    source("blockcv-generate-data.R")
    Xorig = X

    ## Add noise to the sunlight and noise covariates once
    ii = datatype %% 80
    if( ii > 0 ){
      X[,"par"] = Xorig[,"par"] + rnorm(TT, 0, sigmalist[ii])
      for(icol in 3:ncol(X)){
        X[,icol] = Xorig[,icol] + rnorm(TT, 0, sigmalist[ii])
      }
    }
    ## Save the data, bundled into |datobj|
    datobj = list(ylist = ylist, countslist = countslist,
                  X = X, Xorig = Xorig)
    saveRDS(datobj, file=file.path(destin,
                                   paste0("data-", datatype, "-", isim, ".Rds")))
  }, mc.cores = mc.cores)
}
