## Generate artificial datasets once;
outputdir = "/scratch/sangwonh/output"
datadir = "/scratch/sangwonh/data"
scriptdir = "/scratch/sangwonh/scripts"
codedir = "/scratch/sangwonh/repos/flowcy/flowcy"
load_all(codedir)
source(file.path(scriptdir, "blockcv-helpers.R"))
datatype = 9
blocktype = 2
nsim = 20

## For each numclust and isim
numclusts = c(2:8)
## numclusts=8
for(numclust in numclusts){
  destin = build_destin(outputdir, blocktype, datatype, numclust)
  for(isim in 21:30){
  ## for(isim in 1:nsim){

    ## Generate new datasets
    set.seed(numclust * 100 + isim) ## This is the key for the dataset
    source(file.path(scriptdir, "blockcv-generate-data.R"))

    ## Save the responses in a fixed location
    if(!dir.exists(file.path(destin, paste0("sim-", isim)))){
      dir.create(file.path(destin, paste0("sim-", isim)))
    }
    save(X, file = file.path(destin, paste0("sim-", isim), "X.Rdata"))
    save(ylist, file = file.path(destin, paste0("sim-", isim), "ylist.Rdata"))
    save(countslist, file = file.path(destin, paste0("sim-", isim), "countslist.Rdata"))
  }
}
