##' Form the destination folder name. See blockcv-generate-data.R for more
##' details about argument types.
##' @param outputdir Output directory.
##' @param blocktype Block type.
##' @param datatype Data type.
##' @param numclust Number of clusters
##' @param append Anything else that needs to be appended.
##'
##' @return Full destination folder path.
##' @export
build_destin <- function(outputdir, blocktype, datatype, numclust, append=NULL){
  ## filename = paste0("blockcv-blocktype-", blocktype,
  ##                   "-datatype-", datatype, "-numclust-", numclust)
  filename = paste0("blockcv-", blocktype,
                    "-", datatype, "-", numclust)
  if(!is.null(append)) filename = paste0(filename, "-", append)
  destin = file.path(outputdir, filename)
}



## ##' Helper function to lag a vector.
## lagpad <- function(x, k) {
##   if (k>0) {
##     return (c(rep(NA, k), x)[1 : length(x)] );
##   }
##   else {
##     return (c(x[(-k+1) : length(x)], rep(NA, -k)));
##   }
## }


## ##' (to add to package) Add step function bases demarcated by transition region
## ##' crossings.
## add_transition <- function(X, lat){

##   stopifnot(length(lat) == nrow(X))

##   ## Define the crossings of the transition zone (latitude of 37).
##   TT = nrow(X)
##   ind = rep(0, TT)
##   north = which(lat > 37)
##   regions = sapply(1:3, function(ii){
##     endpt = c(0,range(north),TT)[ii:(ii+1)]
##     (endpt[1]+1):endpt[2]
##   })

##   ## Create the TF bases.
##   bases0 = lapply(regions, function(reg){
##     vec = rep(0, TT)
##     vec[reg] = 1
##     vec
##   })

##   ## Add them to X and return
##   bases = do.call(cbind, bases0[-1])
##   X = cbind(bases, X)
##   colnames(X)[1:ncol(bases)] = paste0("b", 1:ncol(bases))
##   return(X)
## }


## ##' (To add to package) add step function to the lags
## add_lagpar <- function(X, lags){
##   ## lags = c(0,3,6,9,12)
##   par = scale(X[,"par"])
##   par = par - min(par)
##   par = par / max(par)
##   parlist = lapply(lags, function(lag)lagpad(par, lag))
##   ## Make the additional columns
##   dat = do.call(cbind, parlist)
##   X = cbind(dat, X)
##   colnames(X)[1:ncol(dat)] = c(paste0("p", 1:ncol(dat)))
##   return(X)
## }
