## Synopsis: testing out the blocked CV approach
## codedir = "/scratch/sangwonh/repos/flowmix/flowmix"
outputdir = "/scratch/sangwonh/output"
datadir = "/scratch/sangwonh/data"
## datadir = "~/Dropbox/research/usc/flow-cytometry/data"
library(flowmix)
library(dplyr)

###################
## Default Setup ##
###################
maxdev_inv = 2
append = NULL
fitonly = FALSE ## Only do refitting
random_order = FALSE ## Make iilist into random order
maxres = FALSE
test = FALSE
sigma_fac = 1
numclust = 5
cv_gridsize = 10
nfold = 5
nrep = 10
blocksize = 20 ## Size for hourlong blocks (block type 2).
limit_y = 0

############################################
### Read in some command line arguments ####
############################################
parse_args(args = commandArgs(trailingOnly = TRUE), verbose=TRUE)
master_node = (arraynum == 1)

################################################
### See how many cores are available to use ####
################################################
mc.cores = system("echo $SLURM_CPUS_PER_TASK", intern=TRUE) %>% as.numeric()
print("mc.cores")
print(mc.cores)
mc.cores = 1

#############################################
#### Create directory to save results in ####
#############################################
## destin = file.path(outputdir, paste0("gradients1-numclust", numclust, "-blocksize", blocksize, "-maxdevinv", maxdev_inv))##, "-limity", limit_y)) ##blocktype, datatype, numclust, append)
## destin = file.path(outputdir, paste0("gradients1-numclust", numclust, "-blocksize", blocksize, "-maxdevinv", maxdev_inv, "-limity", limit_y)) ##blocktype, datatype, numclust, append)
if(!dir.exists(destin)){
  dir.create(destin, recursive = TRUE)
  cat("Creating destin: ", destin, fill=TRUE)
}
cat("All output goes out to destin: ", destin, fill = TRUE)

##########################################
## Load data (ylist, countslist and X)####
##########################################
## load(file = file.path(datadir, "KOK1606.Rdata"))
## load(file = file.path(datadir, "KOK1606-binned.Rdata"), verbose = TRUE)
## load(file = file.path(datadir, paste0("KOK1606-binned-1d-diam.Rdata")), verbose=TRUE)
load(file = file.path(datadir, paste0("MGL1704-binned-hourly.Rdata")), verbose=TRUE) ## Fix this.

## Do whatever processing was done before, from generate-data.R
X = X[-(1:12),] %>% as.matrix()
ylist = ylist[-(1:12)]
## countslist = qclist[-(1:12)]
countslist = countslist[-(1:12)]


## TODO: just load 20 different bootstrap replicates pre-made and saved, instead of doing it on the fly
load(file = file.path(datadir, paste0("MGL1704-binned-hourly-bootstrap1-1.Rdata")), verbose=TRUE) ## Fix this.
load(file = file.path(datadir, paste0("MGL1704-binned-hourly-bootstrap1-20.Rdata")), verbose=TRUE) ## Fix this.
bootrep = 15
load(file = file.path(datadir, paste0("MGL1704-binned-hourly-bootstrap1-", bootrep, ".Rdata")), verbose=TRUE) ## Fix this.


#########################################################
## Calculate maximum lambda (save as "maxres.Rdata") ####
#########################################################
if(maxres){
  maxres = get_max_lambda(destin, "maxres.Rdata", ylist = ylist,
                          countslist = countslist, X = X, numclust = numclust,
                          maxdev = maxdev, max_prob_lambda = 1,
                          max_mean_lambda = 1, verbose = TRUE)
  q()
  if(fitonly) q()
} else {
  load(file.path(destin, "maxres.Rdata"), verbose=TRUE)
}
## TODO: simply load the maxres from before

prob_lambdas =  logspace(min = 0.0001, max = maxres$alpha, length = cv_gridsize)
mean_lambdas = logspace(min = 0.0001, max = maxres$beta,  length = cv_gridsize)

######################
## Create CV folds ###
######################
folds = make_cv_folds(ylist, nfold, verbose = TRUE, blocksize = blocksize)

####################################################
### Save meta data once PER bootstrap replicate ####
####################################################
if(master_node){
  save(folds,
       nfold,
       nrep, ## Added recently
       cv_gridsize,
       mean_lambdas,
       prob_lambdas,
       ylist, countslist, X,
       ## Save the file
       file = file.path(destin, 'meta.Rdata'))
  print(paste0("wrote meta data to ", file.path(destin, 'meta.Rdata')))
}

###################################################
## Run ONE CV job (ialpha, ibeta, ifold, irep) ####
###################################################
if(!fitonly){

  ## Chunk the jobs up.
  iimat = make_iimat(cv_gridsize, nfold, nrep)
  iimax = nrow(iimat)
  ends = round(seq(from = 0, to = iimax, length = arraynum_max + 1))
  iilist = Map(function(a,b){ (a+1):b}, ends[-length(ends)], ends[-1])

  ## Run one chunk of jobs in one machine, in parallel.
  parallel::mclapply(iilist[[arraynum]], function(ii){

    ## Indices defining the jobs to run
    ialpha = iimat[ii,"ialpha"]
    ibeta = iimat[ii,"ibeta"]
    ifold = iimat[ii,"ifold"]
    irep = iimat[ii,"irep"]

    cat('(ialpha, ibeta, ifold, irep)=', c(ialpha, ibeta, ifold, irep), fill=TRUE)

    if(min(mean_lambdas) == 1E-5) stop("mean lambdas are not as expected") ## temporary

    one_job(ialpha = ialpha,
            ibeta = ibeta,
            ifold = ifold,
            irep = irep,
            folds = folds,
            destin = destin,
            mean_lambdas = mean_lambdas,
            prob_lambdas = prob_lambdas,
            ## Arguments for flowmix()
            ylist = ylist, countslist = countslist, X = X,
            sigma_fac = sigma_fac,
            ## Additional arguments for flowmix()
            numclust = numclust,
            maxdev = maxdev,
            ## Temporary
            verbose = TRUE)


    return(NULL)
  }, mc.cores = mc.cores)
}

if(fitonly){

  ## Chunk the jobs up
  iimat = make_iimat_small(cv_gridsize)
  iimax = nrow(iimat)
  ends = round(seq(from = 0, to = nrow(iimat), length = arraynum_max + 1))
  iilist = Map(function(a,b){ (a+1):b}, ends[-length(ends)], ends[-1])


  ## Run one chunk of jobs in one machine, in parallel.
  mclapply(iilist[[arraynum]], function(ii){

    ialpha = iimat[ii, "ialpha"]
    ibeta = iimat[ii, "ibeta"]

    print(paste0('refitting (ialpha, ibeta)=(', ialpha, ",", ibeta, ")"))

    one_job_refit(ialpha = ialpha, ibeta = ibeta, destin = destin,
                  mean_lambdas = mean_lambdas, prob_lambdas = prob_lambdas,
                  ## Arguments to covarem()
                  ylist = ylist, countslist = countslist, X = X,
                  sigma_fac = sigma_fac,
                  numclust = numclust,
                  maxdev = maxdev,
                  nrep = nrep)

    return(NULL)
  }, mc.cores=round(mc.cores*0.75))
}
