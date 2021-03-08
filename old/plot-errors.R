## Synopsis: Plot estimation errors away from the edge.
codedir = "/home/rcf-proj3/sh7/repos/flowcy/flowcy"
outputdir = "/home/rcf-proj3/sh7/output"
datadir = "/home/rcf-proj3/sh7/data"


## ## For the laptop
## datadir = "~/repos/cruisedat/export"
## ## codedir = "~/repos/flowcy/flowcy"
la(codedir)

## Generate data
parse_args(args = commandArgs(trailingOnly = TRUE))
source("blockcv-generate-data.R")

## Model settings (some arbitrary lambda values, for now)
numclust = 10
mean_lambda = 0.001
pie_lambda = 0.001
nfold = 5

## Make CV folds
if(blocktype == 1) folds = blockcv_make_folds(ylist, nfold, verbose=TRUE)
if(blocktype == 2) folds = blockcv_hourlong_make_folds(ylist, nfold, verbose=TRUE)

## Define destination to save to
destin = file.path(outputdir, paste0("plot-error-blocktype-", blocktype, "-datatype-", datatype))
if(!dir.exists(destin)){
  dir.create(destin, recursive = TRUE)
  cat("Creating destin: ", destin, fill=TRUE)
}

## Span cores on available clusters.
cl = get_cl()
parallel::clusterCall(cl, function(){
  devtools::load_all("/home/rcf-proj3/sh7/repos/flowcy/flowcy")})
parallel::clusterExport(cl, ls())


## For each fold, calculate the test score (test likelihood).
parallel::parLapplyLB(cl, 1:nfold, function(ifold){

## 4
## lapply(1:nfold, function(ifold){
    printprogress(ifold, nfold, fill=TRUE)

    ## Train a model on the first half.
    test.ind = unlist(folds[ifold])
    ytrain = ylist[-test.ind]
    Xtrain = X[-test.ind,]
    ytest = ylist[test.ind]
    Xtest = X[test.ind,]
    counttrain = countslist[-test.ind]
    counttest = countslist[test.ind]

    ## Fit the algorithm once
    res = covarem(ylist = ytrain,
                  countslist = counttrain,
                  X = Xtrain,
                  pie_lambda = pie_lambda,
                  mean_lambda = mean_lambda,
                  numclust = numclust,
                  verbose = TRUE,
                  nrep = 1)
    ## niter = 2, ## temporary
    ## nrep = 1) ## temporary

    ## Make test predictions.
    pred = predict(res, newx = Xtest)

    ## Calculate test likelihoods.
    test.loglik = sapply(1:length(test.ind), function(ii){
      objective_subset(ii,
                       mu = pred$newmn,
                       pie = pred$newpie,
                       sigma = pred$sigma,
                       ylist = ytest,
                       countslist = counttest)
    })



    ## return(list(res = res, pred = pred, test.scores = test.scores))
    save(res = res, pred = pred, test.loglik = test.loglik,
         test.ind, ##forgot to do this! darn!
         file = file.path(destin, paste0("fitted-model-ifold", ifold, ".Rdata")))
})



## ##############################################
## ##  Load and plot results for type 1 block  ##
## ##############################################
## la('flowcy')
## ## mean_lambda = 0.001
## ## pie_lambda = 0.001
## mean_lambda = 0.0001
## pie_lambda = 0.0001
## numclust = 4
## ## par(mfrow=c(2, 4))
## par(mfrow=c(4, 4))

## set.seed(0)
## obj_orig = profvis::profvis({
## res_orig = covarem(ylist = ylist,
##               countslist = countslist,
##               X = X,
##               pie_lambda = pie_lambda,
##               mean_lambda = mean_lambda,
##               numclust = numclust,
##               verbose = TRUE,
##               ## admm_rho = 100,
##               ## admm_niter = 3000,
##               ## admm_local_adapt_niter = 1,
##               admm_rho = 0.01,
##               admm_niter = 1000,
##               admm_local_adapt_niter = 10,
##               admm_err_abs = 1E-4,
##               ## niter = 100,
##               ## niter = 50,
##               niter = 10,
##               nrep = 1,
##               always_first_iter=TRUE)
## })


## set.seed(0)
## obj_warmstart = profvis::profvis({
## res_warmstart = covarem(ylist = ylist,
##                         countslist = countslist,
##                         X = X,
##                         pie_lambda = pie_lambda,
##                         mean_lambda = mean_lambda,
##                         numclust = numclust,
##                         verbose = TRUE,
##                         ## admm_rho = 100,
##                         ## admm_niter = 3000,
##                         ## admm_local_adapt_niter = 1,
##                         admm_rho = 0.01,
##                         admm_niter = 1000,
##                         admm_local_adapt_niter = 10,
##                         admm_err_abs = 1E-4,
##                         ## niter = 100,
##                         ## niter = 50,
##                         niter = 10,
##                         nrep = 1,
##                         always_first_iter=FALSE)
## })

## summary(res_orig$beta[[1]] - res_warmstart$beta[[1]])
## summary(res_orig$beta[[2]] - res_warmstart$beta[[2]])
## summary(res_orig$beta[[3]] - res_warmstart$beta[[3]])
## plot(x=res_orig$mn, y= res_warmstart$mn)

## mtext(side=1,outer=TRUE, 'ADMM iteration')
## mtext(side=2,outer=TRUE, 'ADMM objective')


## dev.print(pdf, '~/Pictures/em-beta-mstep-warmstart-new.pdf', width=30, height=30)

## Do the fitted coefficients match?


## Currently, 30 seconds per iteration




## ############################################
## ## Load and plot results for type 1 block ##
## ############################################
blocktype = 1
nfold = 5
## datatype = 21
for(datatype in c(21,24,25)){
datadir = "~/repos/cruisedat/export"
##
if(blocktype == 1) folds = blockcv_make_folds(ylist, nfold, verbose=TRUE)
if(blocktype == 2) folds = blockcv_hourlong_make_folds(ylist, nfold, verbose=TRUE)
plot(NA, xlim = c(0, length(ylist)), ylim = c(1,3))## range(test.loglik)

loglik.list = lapply(c(1,2,3,4,5), function(ifold){
  destin = paste0("~/repos/flowcy/output/plot-error-blocktype-", blocktype,"-datatype-", datatype)
  load(file = file.path(destin, paste0("fitted-model-ifold", ifold, ".Rdata")))
  nn = length(test.loglik)
  if(ifold %in% c(1,5)){
    if(ifold==1){
      list(test.loglik[nn:((nn/2) + 1)])
    } else {
      list(test.loglik[1:(nn/2)])
    }
  } else {
    list(test.loglik[1:(nn/2)], test.loglik[nn:((nn/2) + 1)])
  }
})
loglik.list = unlist(loglik.list, recursive = FALSE)


## doing some fuddling because not all lengths are equal
lens = sapply(loglik.list, length)
loglik.list = lapply(loglik.list, function(a){if(length(a)==49) return(c(a,NA)) else return(a)})
loglik.mat = do.call(cbind, loglik.list)

## Make plot
## datatype = 21
## blocktype = 1
destin = paste0("~/repos/flowcy/output/plot-error-blocktype-", blocktype,"-datatype-", datatype)
pdf(file.path(destin, paste0("plot-error-blocktype-", blocktype, "-datatype-", datatype, ".pdf")), width=8, height=8)
matplot(loglik.mat, type='o', pch=1, ylab = "Test likelihoods", xlab = "Distance away from closest training time point")
obj = ksmooth(y = as.numeric(loglik.mat),
              x = as.numeric(matrix(rep(1:50, ncol(loglik.mat)), byrow=FALSE, ncol=ncol(loglik.mat))),
              x.points = 1:50, bandwidth = 3, kernel = "normal")
lines(obj, type='l', lwd=4)
title(main = "Five equally sized \n contiguous time blocks")
graphics.off()
}

############################################
## Load and plot results for type 2 block ##
############################################
blocktype = 2
## datatype = 24
for(datatype in c(21,24,25)){
nfold = 5
if(blocktype == 1) folds = blockcv_make_folds(ylist, nfold, verbose=TRUE)
if(blocktype == 2) folds = blockcv_hourlong_make_folds(ylist, nfold, verbose=TRUE)
plot(NA, xlim = c(0, length(ylist)), ylim = c(1,3))
loglik.list = lapply(c(1,2,3,4,5), function(ifold){
  if(ifold==1) return(NULL)
  destin = paste0("~/repos/flowcy/output/plot-error-blocktype-", blocktype,"-datatype-", datatype)
  load(file = file.path(destin, paste0("fitted-model-ifold", ifold, ".Rdata")))
  ## test.ii = folds[[ifold]]
  ## abline(v=(test.ii), lwd=3, col=ifold)
  ## points(y=test.loglik, x=test.ii, pch=16, cex=1)
  ## list(test.loglik[1:(nn/2)], test.loglik[nn:((nn/2) + 1)])
  return(test.loglik)
})

## Gather results into a matrix
inds = unlist(lapply(c(0,20,40,60,80), function(start) list(start + 1:10, start + 20:11)), recursive=FALSE)
loglik.mat.list = lapply(c(1,2,3,4,5), function(ifold){
  if(ifold==1) return(NULL) ## first fold failed.
  print(ifold)
  test.loglik = loglik.list[[ifold]]
  test.inds = folds[[ifold]]
  do.call(cbind, lapply(inds, function(ind){
    test.ind = test.inds[ind]
    if(1 %in% test.ind  | 500 %in% test.ind){ ## Making sure that the edges aren't included.
      return(NULL)
    } else {
      return(test.loglik[ind])
    }
  }))
})

loglik.mat = do.call(cbind, loglik.mat.list)
## loglik.mat = do.call(cbind, lapply(inds, function(ind) test.loglik[ind]))
rownames(loglik.mat) = paste("dist=", 1:10)


## datatype = 21
## blocktype = 2
destin = paste0("~/repos/flowcy/output/plot-error-blocktype-", blocktype,"-datatype-", datatype)
pdf(file.path(destin, paste0("plot-error-blocktype-", blocktype, "-datatype-", datatype, ".pdf")), width=8, height=8)
matplot(loglik.mat, type='o', pch=1, ylab = "Test likelihoods", xlab = "Distance away from closest training time point")
        ## col=col)
obj = ksmooth(y = as.numeric(loglik.mat),
              x = as.numeric(matrix(rep(1:10, ncol(loglik.mat)), byrow=FALSE, ncol=ncol(loglik.mat))),
              x.points = 1:10, bandwidth = 1, kernel = "normal")
lines(obj, type='l', lwd=4)
title(main = "Hour-long blocks")
graphics.off()
}
