## Synopsis: testing out the blocked CV approach
## codedir = "/home/rcf-proj3/sh7/repos/flowcy/flowcy"
## codedir = "/staging/sh7/sangwonh/repos/flowcy/flowcy"
codedir = "/scratch/sangwonh/repos/flowmix/flowmix"
## outputdir = "/home/rcf-proj3/sh7/output"
## outputdir = "/staging/sh7/sangwonh/output"

outputdir = "/scratch/sangwonh/output"
## datadir = "/home/rcf-proj3/sh7/data"
## datadir = "/staging/sh7/sangwonh/data"
datadir = "/scratch/sangwonh/data"
source("blockcv-helpers.R")
load_all(codedir)

###########
## Setup ##
###########
maxdev = 0.5
append = NULL
blocksize = 20 ## Size for hourlong blocks (block type 2).
fitonly = FALSE ## Only do refitting
random_order = FALSE ## Make iilist into random order
maxres_once = FALSE
test = FALSE
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

################################################
### See how many cores are available to use ####
################################################
mynodename = strsplit(system("scontrol show hostname $SLURM_NODELIST | paste -d, -s", intern=TRUE),
                           ",")[[1]]
mycommand = "sinfo -o '%40N %c' --Node | grep "
## mycommand = "sinfo -o '%40N %c' --Node --long | grep "
mystring = system(paste0(mycommand, mynodename), intern=TRUE)
mc.cores = as.numeric(strsplit(mystring, "\\s+")[[1]][2])
print("Number of cores:")
print(mc.cores)

#############################################
#### Create directory to save results in ####
#############################################
destin = file.path(outputdir, paste0("experiment-type-", experiment_type))
if(!dir.exists(destin)){
  dir.create(destin, recursive = TRUE)
  cat("Creating destin: ", destin, fill=TRUE)
} else {
  cat("All output goes out to destin: ", destin, fill = TRUE)
}


################################
## Generate data and folds ####
################################
source("blockcv-generate-data.R")

## Instead, just load the data here.

if(experiment_type == 1){

  ## Load the 3d data for a 10-cluster analysis
  dat = readRDS(file = file.path(datadir, "MGL1704-hourly-paper.RDS"))
  list2env(dat, globalenv())

} else if (experiment_type == 2){

  ## Load the 1d data for a 5-cluster analysis
  dat = readRDS(file = file.path(datadir, "MGL1704-hourly-paper-1d-diam.RDS"))
  list2env(dat, globalenv())

} else if (experiment_type == 3){

  ## Generate realistic 5-cluster 1d data, in a few steps

  ## 1. Load data.
  obj = flowmix::generate_data_1d_pseudoreal_from_cv(datadir = datadir) ## loads a file called "1d-cvres.rds"
  ylist = obj$ylist
  X = obj$X
  countslist = obj$countslist

  ## 2. Bin the data
  dat.grid = flowmix::make_grid(ylist, gridsize = 40)
  obj = flowmix::bin_many_cytograms(ylist, dat.grid, mc.cores = 1, verbose = TRUE)
  ylist = lapply(obj$ybin_list, cbind)
  countslist = obj$counts_list

  ## Other setup
  maxdev = 0.1155245
  final_range = lapply(ylist, range) %>% unlist %>% range
  sigma_fac = (final_range[2] - final_range[1])/numclust

} else if (experiment_type == 4){

  ## Generate artificial 2-cluster 1d simulation data, in a few steps

  ## 1. Generate fake 1d data
  obj = generate_data_1d_pseudoreal(bin = FALSE, datadir = datadir, nt1 = 200, beta_par = 0.3, p = 10)
  X = obj$X
  ylist = obj$ylist
  countslist = obj$countslist
  TT = length(ylist)

  ## 2. Bin with just counts
  dat.grid = flowcy::make_grid(ylist, gridsize = 40)
  obj = flowcy::bin_many_cytograms(ylist, dat.grid, mc.cores = 1, verbose = TRUE)
  ylist = lapply(obj$ybin_list, cbind)
  countslist = obj$counts_list

  ## Check about noise level
  sigmalist = seq(from = 3 * 0.1, to = 3 * 0.9, length = 9)
  stopifnot(noise_ii %in% c(1:9))
}

## For simulations, X needs to be saved before modification.
Xorig = X

################################
## Calculate maximum lambda ####
################################
if(maxres_once){
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
  q()
} else {
  load(file.path(destin, "maxres.Rdata"), verbose=TRUE)
}

pie_lambdas =  logspace(min = 0.0001, max=maxres$alpha, length=cv_gridsize)
mean_lambdas = logspace(min = 0.0001, max=maxres$beta, length=cv_gridsize)


######################
## Create CV folds ###
######################
folds = blockcv_hourlong_make_folds(ylist, nfold, verbose = TRUE, blocksize = blocksize)


############################
### Save meta data once ####
############################
if(master_node){
  save(folds,
       nfold,
       nrep, ## Added recently
       cv_gridsize,
       mean_lambdas,
       pie_lambdas,
       ylist, countslist, X,
       ## Save the file
       file = file.path(destin, 'meta.Rdata'))
  print(paste0("wrote meta data to ", file.path(destin, 'meta.Rdata')))
}

#####################################
## Temporary, for simulations: ######
#####################################


###################################################
## Run ONE CV job (ialpha, ibeta, ifold, irep) ####
###################################################
if(!fitonly){
  iimat = make_iimat(cv_gridsize, nfold, nrep)

  ## Temporary:
  iimat = iimat[which(iimat[,"irep"] > 5),]

  ## Chunk the jobs up.
  iimax = nrow(iimat)
  ends = round(seq(from=0, to=iimax, length=arraynum_max+1))
  iilist = Map(function(a,b){ (a+1):b}, ends[-length(ends)], ends[-1])

  ## If applicable, randomize the order
  if(random_order){
    ## stop("Don't use random order!!") ## In fact, retire this soon.
    iirand = sample(1:iimax, iimax)
    iilist = lapply(iilist, function(ii){iirand[ii]})
  }

  cat("maxdev is:", maxdev, fill=TRUE)

  ## Run the jobs in parallel, in one machine.
  mclapply(iilist[[arraynum]], function(ii){

    ## Indices defining the jobs to run
    ialpha = iimat[ii,"ialpha"]
    ibeta = iimat[ii,"ibeta"]
    ifold = iimat[ii,"ifold"]
    irep = iimat[ii,"irep"]

    cat('(ialpha, ibeta, ifold, irep)=', c(ialpha, ibeta, ifold, irep), fill=TRUE)

    if(min(mean_lambdas)==1E-5) stop("mean lambdas are not as expected") ## temporary

    for(isim in 1:nsim){

      ## Add noise to X, if applicable
      ## if(sim & (80 <= datatype & datatype < 90)){
      if(experiment_type == 4){

        ## Add noise to the sunlight covariate
        if( ii > 0 ){
          X[,"par"] = Xorig[,"par"] + rnorm(TT, 0, sigmalist[noise_ii])

          ## Add noise to the spurious covariates by the same amount as well.
          for(icol in 3:ncol(X)){
            X[,icol] = Xorig[,icol] + rnorm(TT, 0, sigmalist[noise_ii])
          }
        }
      }

      one_job(ialpha = ialpha,
              ibeta = ibeta,
              ifold = ifold,
              irep = irep,
              folds = folds,
              destin = destin,
              mean_lambdas = mean_lambdas,
              pie_lambdas = pie_lambdas,
              ## Arguments for covarem()
              ylist = ylist, countslist = countslist, X = X,
              sigma_fac = sigma_fac,
              ## Additional arguments for covarem(), for ellipsis.
              numclust = numclust,
              maxdev = maxdev,
              ## Simulations? Yes or no
              sim = sim,
              isim = isim,
              verbose=TRUE)
    }
    return(NULL)
  }, mc.cores=round(mc.cores*0.75))
}
if(fitonly){

  ## Chunk the jobs up
  iimat = make_iimat_small(cv_gridsize)
  iimax = nrow(iimat)
  ends = round(seq(from=0, to=nrow(iimat), length=arraynum_max+1))
  iilist = Map(function(a,b){ (a+1):b}, ends[-length(ends)], ends[-1])

  ## If applicable, randomize the order
  if(random_order){
    iirand = sample(1:iimax, iimax)
    iilist = lapply(iilist, function(ii){iirand[ii]})
  }

  ## parallel::parLapplyLB(cl, 1:iimax, function(ii){
  mclapply(iilist[[arraynum]], function(ii){

    ialpha = iimat[ii, "ialpha"]
    ibeta = iimat[ii, "ibeta"]

    print(paste0('refitting (ialpha, ibeta)=(', ialpha, ",", ibeta, ")"))
    for(isim in 1:nsim){

      ## Add noise to X, if applicable
      if(sim & (80 <= datatype & datatype < 90)){
        ## sigmalist = seq(from = 0.1, to = 0.9, length = 9)
        sigmalist = seq(from = 3 * 0.1, to = 3 * 0.9, length = 9)
        ii = datatype %% 80  ## ii = 9
        if(ii > 0){
          X[,"par"] = Xorig[,"par"] + rnorm(TT, 0, sigmalist[ii])
        }
      }

      one_job_refit(ialpha = ialpha, ibeta = ibeta, destin = destin,
                    mean_lambdas = mean_lambdas, pie_lambdas = pie_lambdas,
                    ## Arguments to covarem()
                    ylist = ylist, countslist = countslist, X = X,
                    sigma_fac = sigma_fac,
                    numclust = numclust,
                    maxdev = maxdev,
                    nrep = nrep,
                    sim = sim,
                    isim = isim)
    }
    return(NULL)
  }, mc.cores=round(mc.cores*0.75))
}
