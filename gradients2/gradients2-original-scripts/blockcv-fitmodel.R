## Synopsis: testing out the blocked CV approach
codedir = "/home/rcf-proj3/sh7/repos/flowcy/flowcy"
outputdir = "/home/rcf-proj3/sh7/output"
datadir = "/home/rcf-proj3/sh7/data"
load_all(codedir)
source("blockcv-helpers.R")

###########
## Setup ##
###########
nfold = 5
cv_gridsize = 10
maxdev = 0.5
append = NULL
coarse = FALSE
iimin = iimax = NULL
blocksize = 20 ## For hourlong blocks.
## admm_rho = 0.1
## admm_niter = 1000

############################################
### Read in some command line arguments ####
############################################
parse_args(args = commandArgs(trailingOnly = TRUE), verbose=TRUE)
stopifnot(!is.na(blocktype) & !is.na(datatype))
if(is.null(iimin) & is.null(iimax)){
  iirange= NULL
} else {
  iirange = iimin:iimax
}

#############################################
#### Create directory to save results in ####
#############################################
## append = "singlecore"
destin = build_destin(outputdir, blocktype, datatype, append)
if(!dir.exists(destin)){
  stop("Directory needs to exist already, for refitting!")
} else {
  cat("All output goes out to destin: ", destin, fill = TRUE)
}

################################
## Generate data and folds ####
################################
source("blockcv-generate-data.R")



################################
## Calculate maximum lambda ####
################################
maxres = get_max_lambda(destin, "maxres.Rdata",
                        ylist = ylist,
                        countslist = countslist,
                        X = X,
                        numclust = numclust,
                        maxdev = maxdev,
                        max_lambda_alpha = 10,
                        max_lambda_beta = 10)
pie_lambdas = logspace(min = 0.00001, max=maxres$alpha, length=cv_gridsize)
mean_lambdas = logspace(min = 0.00001, max=maxres$beta, length=cv_gridsize)

######################
## Parallel setup ####
######################
cl = NULL
if(parallel){
  cl = get_cl()
  parallel::clusterExport(cl, ls())## ls_without_ylist)
  ## ls_without_ylist = ls()[-which(ls() %in% "ylist")]
  ## print(ls())
  parallel::clusterCall(cl, function(){devtools::load_all(codedir)})
}

cat("Block CV refitting has started")
for(irep in 1:nrep){
  cat("Refitting for (ialpha, ibeta, irep) = (", ialpha, ibeta, irep, ") in progress.", fill=TRUE)
  filename = paste0(ialpha, "-", ibeta, "-", irep, "-fit.Rdata")
  res = covarem_once(ylist = ylist,
                     countslist = countslist,
                     mean_lambda = mean_lambda[ind_mean_lambda],
                     pie_lambda = pie_lambda[ind_pie_lambda],
                     X = X,
                     maxdev = maxdev,
                     numclust = numclust)
  save(res, file=file.path(destin, filename))
}

## blockcv_fitmodel(cl = cl,
##                  parallel = parallel,
##                  mean_lambdas = mean_lambdas,##rep(0, cv_gridsize),
##                  pie_lambdas = mean_lambdas,##rep(0, cv_gridsize),
##                  maxdev = maxdev,
##                  ylist = ylist,
##                  countslist = countslist,
##                  X = X,
##                  numclust = numclust,
##                  nrep = nrep,
##                  destin = destin)



## The rest of this script will be useful later.

## ##############################
## ## Aggregate the results  ####
## ##############################
## destin = "~/repos/flowcy/output/blockcv-blocktype-1-datatype-6"
## cv.gridsize = 10
## nfold = 5
## la("flowcy")
## obj = blockcv_aggregate(destin = destin, cv.gridsize = cv.gridsize, nfold = nfold)
## cvscoremat_orig = obj$cvscore.mat
## install.packages('pheatmap') # if not installed already
## library(pheatmap)
## cvscoremat = cvscoremat_orig - min(cvscoremat_orig, na.rm=TRUE)
## pdf(file=file.path("~/Desktop/blockcv-1-6-errors.pdf"), width=7, height=7)
## pheatmap(cvscoremat_orig,
##          display_numbers = T,
##          cluster_cols = FALSE,
##          cluster_rows=FALSE)
## graphics.off()
## cvscoremat_orig

## ##  What were the lambdas?
## ialpha = 9
## ibeta = 9
## igrid = 1
## filename = paste0(ialpha, "-", ibeta, "-", igrid, "-cvscore.Rdata")
## print(filename)
## load(file.path(destin, filename))
## mean_lambdas
## plot(pie_lambdas, type='l')
## cvscore
## min(cvscoremat, na.rm=TRUE)

## ## Fit the model for the most regularized solution.

## ## (Load data from  blockcv-speed.R in dropbox)
## set.seed(10)
## la('flowcy')
## obj = covarem(ylist = ylist,
##               X = X,
##               countslist = countslist,
##               numclust = 10,
##               mean_lambda = 1000,
##               pie_lambda = 1000,
##               verbose = TRUE,
##               plot = TRUE,
##               plotdir = "~/Desktop/blockcv-test-figures")

## ## filled.contour(1:(n-1), 1:(n-1), clockwise90(BBt), color=heat.colors)
