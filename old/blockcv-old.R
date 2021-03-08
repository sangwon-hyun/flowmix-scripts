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
fitonly = FALSE

############################################
### Read in some command line arguments ####
############################################
parse_args(args = commandArgs(trailingOnly = TRUE), verbose=TRUE)
stopifnot(!is.na(blocktype) & !is.na(datatype))
if(is.null(iimin) & is.null(iimax)){
  iirange= NULL
} else {
  iirange = (iimin:iimax)
}

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

min_mean_lambdas = maxres$beta*0.9
max_mean_lambdas = maxres$beta*1.2

## mean_lambdas = logspace(min = 0.00001, max=maxres$beta, length=cv_gridsize)
mean_lambdas = seq(from = min_mean_lambdas,
                   to = max_mean_lambdas,
                   length = cv_gridsize)

## ## # Temporary: break here for now
## if(datatype == 6){
##   print("Stopping after maxres, for now.")
##   q()
## }
## ## End of temporary


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

######################
## Create CV folds ###
######################
## fold_filename = file.path(destin, "folds.Rdata")
## if(file.exists(fold_filename)){
##   print(folds)
##   load(fold_filename)
## } else {
  if(blocktype==1) folds = blockcv_make_folds(ylist, nfold, verbose = TRUE)
  if(blocktype==2) folds = blockcv_hourlong_make_folds(ylist, nfold, verbose = TRUE, blocksize = blocksize)
##   save(folds, file = fold_filename)
## }



##################
## Run the CV ####
##################
cat("Block CV has started")
if(!fitonly){
blockcv(cl = cl,
        parallel = parallel,
        folds = folds,
        iirange = iirange,
        cv_gridsize = cv_gridsize,
        mean_lambdas = mean_lambdas, ##seq(from=0, to=1, length = cv_gridsize),
        pie_lambdas = pie_lambdas, ##seq(from=0, to=1, length = cv_gridsize),
        maxdev = maxdev,
        ylist = ylist,
        countslist = countslist,
        X = X,
        numclust = numclust,
        nrep = nrep,
        destin = destin)
}



cat("Block CV refitting has started")
blockcv_fitmodel(cl = cl,
                 destin = destin,
                 cv_gridsize = cv_gridsize,
                 mean_lambdas = mean_lambdas,##rep(0, cv_gridsize),
                 pie_lambdas = mean_lambdas,##rep(0, cv_gridsize),
                 ## Arguments to covarem
                 maxdev = maxdev,
                 ylist = ylist,
                 countslist = countslist,
                 X = X,
                 numclust = numclust,
                 nrep = nrep)



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
