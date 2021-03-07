## All directories
outputdir = "/scratch/sangwonh/output"
datadir = "/scratch/sangwonh/data"
scriptdir = "/scratch/sangwonh/scripts/gradients2" ## the current script's location
library(flowmix)


## ## Temporary settings
## DF = 5
## isim = 1
## nsim = 5
## noisetype_ii=2
## summarize = 0
## ## End of temporary settings

### Read in some command line arguments
parse_args(args = commandArgs(trailingOnly = TRUE), verbose=TRUE)

## Basic setup
maxdev = 1.5
numclust = 2
cv_gridsize = 10
nfold = 5
nrep = 10

skew_alpha_list = alphas = seq(from=0, to =2, by = 0.5)

## Basic setup
noisetype = c("gaussian", "heavytail", "skewed")[noisetype_ii]

folder1 = paste0("heavytail-", noisetype)

if(noisetype == "heavytail"){
  folder1 = paste0(folder1, "-df-", DF)
}
if(noisetype == "skewed"){
  skew_alpha = skew_alpha_list[ialpha]
  folder1 = paste0(folder1, "-ialpha-", ialpha)
}

## Form destination folder
folder2 = paste0("sim-", isim)
destin = file.path(outputdir, folder1, folder2)
create_destin(destin)


## Summarize
if(run_summarize){
  cv_summary(destin = destin,
             nfold = nfold,
             nrep = nrep,
             save = TRUE,
             filename = "summary.RDS")
  ## Also save in "subsample-summaries" under destin
  obj = readRDS(file = file.path(destin, "summary.RDS"))
  create_destin(file.path(outputdir, folder1,"summaries"))
  saveRDS(obj, file = file.path(outputdir,folder1, "summaries", paste0("summary-", isim, ".RDS")))
  q()
}


## 1. Generate fake 1d data
set.seed(isim)
obj = generate_data_1d_pseudoreal(bin = FALSE,
                                  datadir = datadir,
                                  nt = 200,
                                  beta_par = 0.3,
                                  p = 10,
                                  noisetype = noisetype,
                                  df = DF,
                                  skew_alpha = skew_alpha)
set.seed(NULL)
X = obj$X
ylist = obj$ylist
countslist = obj$countslist
TT = length(ylist)

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

## Summarize
if(FALSE){
 cv_summary(destin = destin,
             nfold = 5,
             nrep = 10,
             save = TRUE,
             filename = "summary.RDS")
  q()
}
