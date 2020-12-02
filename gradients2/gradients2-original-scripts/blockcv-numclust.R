## Synopsis: testing out the blocked CV approach, ONLY FOR THE NUMCLUST simulation
codedir = "/scratch/sangwonh/repos/flowcy/flowcy"
outputdir = "/scratch/sangwonh/output"
datadir = "/scratch/sangwonh/data"
scriptdir = "/scratch/sangwonh/scripts"
source(file.path(scriptdir,"blockcv-helpers.R"))
load_all(codedir)

###########
## Setup ##
###########
maxdev = 0.5
append = NULL
blocksize = 20 ## Size for hourlong blocks (block type 2).
random_order = FALSE ## Make iilist into random order
maxres_once = FALSE
test = FALSE
sigma_fac = 1

## Things related to the simulations
sim = FALSE ## In normal times, this is set to false, and only one round is run.
nsim = 1

############################################
### Read in some command line arguments ####
############################################
parse_args(args = commandArgs(trailingOnly = TRUE), verbose=TRUE)
stopifnot(!is.na(blocktype) & !is.na(datatype))
master_node = (arraynum == 1)


## Some specific things for the numclust simulation
maxdev = 0.1155245

## ## ## ## ## Test settings ############
## arraynum_max=20
## blocktype=2
## datatype=9
## numclust=7
## ## numclust = 2
## cv_gridsize=7
## nfold=5
## nrep=5
## maxres_once=0
## sim=1
## nsim=40
## ## ## ## ## End of Test settings ####

################################################
### See how many cores are available to use ####
################################################

## First check if the job has a #CPU request, from the --cpus-per-task option.
## mycommand = "echo $SLURM_CPUS_PER_TASK"
## mystring = system(paste0(mycommand), intern=TRUE)
## mc.cores = as.numeric(mystring)## extract

## ## If the #CPUs is not specified, we have exclusive access to this node.
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

########################################
## Load data once, only from isim=1 ####
########################################
isim = 1
load(file = file.path(destin, paste0("sim-", isim), "ylist.Rdata"))
load(file = file.path(destin, paste0("sim-", isim), "X.Rdata"))
load(file = file.path(destin, paste0("sim-", isim), "countslist.Rdata"))

################################
## Calculate maximum lambda ####
################################
if(maxres_once){
  print("maxres starting")
  maxres = get_max_lambda(destin, "maxres.rdata",
                          ylist = ylist,
                          countslist = countslist,
                          X = X,
                          numclust = numclust,
                          maxdev = maxdev,
                          max_lambda_alpha = 40,
                          max_lambda_beta =  2)
  print("maxres ending"); q();
} else {
  load(file.path(destin, "maxres.Rdata"), verbose=TRUE)
}
pie_lambdas =  logspace(min = 0.0001, max=maxres$alpha, length=cv_gridsize)
mean_lambdas = logspace(min = 0.0001, max=maxres$beta, length=cv_gridsize)


######################
## Create CV folds ###
######################
if(blocktype==1) folds = blockcv_make_folds(ylist, nfold, verbose = TRUE)
if(blocktype==2) folds = blockcv_hourlong_make_folds(ylist, nfold, verbose = TRUE, blocksize = blocksize)

###################################################
## Run ONE CV job (ialpha, ibeta, ifold, irep) ####
###################################################
iimat = make_iimat(cv_gridsize, nfold, nrep)

## Chunk the jobs up.
iimax = nrow(iimat)
ends = round(seq(from=0, to=iimax, length=arraynum_max+1))
iilist = Map(function(a,b){ (a+1):b}, ends[-length(ends)], ends[-1])
cat("maxdev is:", maxdev, fill=TRUE)

## Run the jobs in parallel, in one machine.
mclapply(iilist[[arraynum]], function(ii){

  ## Indices defining the jobs to run
  ialpha = iimat[ii,"ialpha"]
  ibeta = iimat[ii,"ibeta"]
  ifold = iimat[ii,"ifold"]
  irep = iimat[ii,"irep"]

  for(isim in 1:nsim){

    ## Print progress
    cat('(isim, ialpha, ibeta, ifold, irep)=', c(isim, ialpha, ibeta, ifold, irep), fill=TRUE)

    ## Load data and get range
    rm(ylist)
    rm(countslist)
    rm(X)
    load(file = file.path(destin, paste0("sim-", isim), "ylist.Rdata"))
    load(file = file.path(destin, paste0("sim-", isim), "countslist.Rdata"))
    load(file = file.path(destin, paste0("sim-", isim), "X.Rdata"))
    final_range = lapply(ylist, range) %>% unlist %>% range
    sigma_fac = (final_range[2] - final_range[1])/numclust

    ## Run the model
    ## obj = profvis::profvis({
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
    ## })
    ## save(obj, file="~/profvis-glmnet-discovery.Rdata")
  }
  return(NULL)
}, mc.cores = round(mc.cores * 0.75), mc.preschedule = FALSE)


###############################################
## Run ONE REFIT job (ialpha, ibeta, irep) ####
###############################################
## Chunk the jobs up
iimat = make_iimat_small(cv_gridsize)
iimax = nrow(iimat)
ends = round(seq(from=0, to=nrow(iimat), length=arraynum_max+1))
iilist = Map(function(a,b){ (a+1):b}, ends[-length(ends)], ends[-1])

## parallel::parLapplyLB(cl, 1:iimax, function(ii){
mclapply(iilist[[arraynum]], function(ii){

  ialpha = iimat[ii, "ialpha"]
  ibeta = iimat[ii, "ibeta"]

  for(isim in 1:nsim){

    ## Load data and get range
    rm(ylist)
    rm(countslist)
    rm(X)
    load(file = file.path(destin, paste0("sim-", isim), "ylist.Rdata"))
    load(file = file.path(destin, paste0("sim-", isim), "countslist.Rdata"))
    load(file = file.path(destin, paste0("sim-", isim), "X.Rdata"))
    final_range = lapply(ylist, range) %>% unlist %>% range
    sigma_fac = (final_range[2] - final_range[1])/numclust


    cat('(isim, ialpha, ibeta)=', c(isim, ialpha, ibeta), fill=TRUE)
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
}, mc.cores = round(mc.cores*0.75), mc.preschedule = FALSE)
