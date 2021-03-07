## All directories
outputdir = "/scratch/sangwonh/output"
datadir = "/scratch/sangwonh/data"
scriptdir = "/scratch/sangwonh/scripts/gradients2" ## the current script's location
library(flowmix)

## ## Test settings
iboot = 1
## ## end of test settings



### Read in some command line arguments
summarize = FALSE
parse_args(args = commandArgs(trailingOnly = TRUE), verbose=TRUE)

## Form destination folder
folder1 = "subsample-simulation"
folder2 = "all-simulations"
folder3 = paste0("iboot-", iboot)
folder4 = paste0("orig")
destin = file.path(outputdir, folder1, folder2, folder3, folder4)
create_destin(destin)



## Basic setup
maxdev = 1.5
numclust = 2
cv_gridsize = 10
nfold = 5
nrep = 10

## Summarize
if(summarize){
  cv_summary(destin = destin,
             nfold = nfold,
             nrep = nrep,
             save = TRUE,
             filename = "summary.RDS")
  ## Also save in "subsample-summaries" under destin
  obj = readRDS(file = file.path(destin, "summary.RDS"))
  create_destin(file.path(outputdir, folder1, folder2, folder3, "subsample-summaries"))
  saveRDS(obj, file = file.path(outputdir,folder1, folder2, folder3, "subsample-summaries", paste0("summary.RDS")))
  create_destin(file.path(outputdir, folder1, folder2, "all-summaries"))
  saveRDS(obj, file = file.path(outputdir,folder1, folder2, "all-summaries", paste0("orig-summary-",iboot, ".RDS")))
  q()
}


## 1. Generate fake 1d data
set.seed(100000 * iboot)
obj = generate_data_1d_pseudoreal(bin = FALSE,
                                  datadir = datadir,
                                  nt = 200,
                                  beta_par = 0.3,
                                  p = 10)
set.seed(NULL)
X = obj$X
ylist = obj$ylist
countslist = obj$countslist
TT = length(ylist)

## 2. Bin with just counts
dat.grid = flowmix::make_grid(ylist, gridsize = 40)
obj = flowmix::bin_many_cytograms(ylist, dat.grid, mc.cores = 1, verbose = TRUE)
ylist = lapply(obj$ybin_list, cbind)
countslist = obj$counts_list

## 3. Obtain the maximum regularization parameters
if(file.exists(file.path(destin, "maxres.Rdata"))){
  load(file.path(destin, "maxres.Rdata"), verbose = TRUE)
} else {
  maxres = get_max_lambda(destin, "maxres.Rdata",
                          ylist = ylist,
                          countslist = countslist,
                          X = X,
                          numclust = numclust,
                          maxdev = maxdev,
                          max_prob_lambda = 1,
                          max_mean_lambda = 3,
                          verbose = TRUE)
  q()
}
prob_lambdas =  logspace(min = 0.0001, max=maxres$alpha, length=cv_gridsize)
mean_lambdas = logspace(min = 0.0001, max=maxres$beta, length=cv_gridsize)


## 4. Estimate model
for(refit in c(FALSE, TRUE)){

  ## Make CV job index matrices
  if(!refit){
    iimat = flowmix::make_iimat(cv_gridsize, nfold, nrep)
  } else {
    iimat = flowmix::make_iimat_small(cv_gridsize)
  }

  ## Divide them up into |arraynum_max| parts
  iilist = make_iilist(arraynum_max, iimat) ## TODO: add flowmix:: here
  if(arraynum > length(iilist)) return(NULL)
  iimat = iimat[iilist[[arraynum]],]

  ## Run the jobs in parallel, in one machine.
  cv.flowmix(
      ## Data
      ylist,
      countslist,
      X,
      ## Define the locations to save the CV.
      destin = destin,
      ## Regularization parameter values
      mean_lambdas,
      prob_lambdas,
      iimat,
      ## Other settings
      maxdev,
      numclust,
      nfold,
      nrep,
      refit = refit,
      save_meta = (arraynum == 1),
      mc.cores = 1)
}

