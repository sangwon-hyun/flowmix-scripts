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
maxdev = 1.5 ## NULL ## 1.5 ## Slightly larger than the truth, which is about (2.66-(-0.84)) * 0.3 = 1.05
blocksize = 20 ## Size for hourlong blocks (block type 2).
fitonly = FALSE ## Only do refitting
random_order = FALSE ## Make iilist into random order
maxres_once = FALSE
manual_index = FALSE ## If true, also set values for |ialpha| & |ibeta|.
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
destin = build_destin(outputdir, blocktype, datatype, numclust)
if(!dir.exists(destin)){
  dir.create(destin, recursive = TRUE)
  cat("Creating destin: ", destin, fill=TRUE)
} else {
  cat("All output goes out to destin: ", destin, fill = TRUE)
}

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
  load(file.path(destin, "maxres.Rdata"))
}

pie_lambdas =  logspace(min = 0.0001, max=maxres$alpha, length=cv_gridsize)
mean_lambdas = logspace(min = 0.0001, max=maxres$beta, length=cv_gridsize)

###################################################
## Run ONE CV job (ialpha, ibeta, ifold, irep) ####
###################################################
if(!fitonly){
  iimat = make_iimat(cv_gridsize, nfold, nrep)

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

      ## Load the data
      datobj = readRDS(file=file.path(destin,
                                      paste0("data-",
                                             datatype, "-", isim, ".Rds")))
      ylist = datobj$ylist
      countslist = datobj$countslist
      X = datobj$X

      ######################
      ## Create CV folds ###
      ######################
      if(blocktype==1) folds = blockcv_make_folds(ylist, nfold, verbose = TRUE)
      if(blocktype==2) folds = blockcv_hourlong_make_folds(ylist, nfold, verbose = TRUE, blocksize = blocksize)


      ############################
      ### Save meta data once ####
      ############################
      if(master_node & isim == 1 & ii==1){
        ## args = list(...)
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
              isim = isim)
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

  ## A temporary thing: run only the refit jobs that I want.
  if(manual_index){
    if(arraynum>2){
      cat("Only arraynum == 1 needed for manual indexing; quitting this node")
      q()
    }
    iirow = rbind(c(ialpha, ibeta))
    colnames(iirow) =  c("ialpha", "ibeta")
    iilist = list(iirow)
  }

  ## parallel::parLapplyLB(cl, 1:iimax, function(ii){
  mclapply(iilist[[arraynum]], function(ii){

    ialpha = iimat[ii, "ialpha"]
    ibeta = iimat[ii, "ibeta"]

    print(paste0('refitting (ialpha, ibeta)=(', ialpha, ",", ibeta, ")"))
    for(isim in 1:nsim){

      ## Load the data
      datobj = readRDS(file=file.path(destin,
                                      paste0("data-",
                                             datatype, "-", isim, ".Rds")))
      ylist = datobj$ylist
      countslist = datobj$countslist
      X = datobj$X

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
