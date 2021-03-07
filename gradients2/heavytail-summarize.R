## All directories
outputdir = "/scratch/sangwonh/output"
datadir = "/scratch/sangwonh/data"
scriptdir = "/scratch/sangwonh/scripts/gradients2" ## the current script's location
library(flowmix)

### Read in some command line arguments
parse_args(args = commandArgs(trailingOnly = TRUE), verbose=TRUE)


## Basic setup
maxdev = 1.5
numclust = 2
cv_gridsize = 10
nfold = 5
nrep = 10

## Basic setup
noisetype = c("gaussian", "heavytail", "skewed")[noisetype_ii]
folder1 = paste0("heavytail-", noisetype)
if(noisetype == "heavytail"){
  folder1 = paste0(folder, "-df-", DF)
}
if(noisetype == "skewed"){
  skew_alpha = skew_alpha_list[ialpha]
  folder1 = paste0(folder, "-ialpha-", ialpha)
}

## Form destination folder
folder2 = paste0("sim-", isim)
destin = file.path(outputdir, folder1, folder2)


## Summarize
cv_summary(destin = destin,
            nfold = 5,
            nrep = 10,
            save = TRUE,
            filename = "summary.RDS")

## Also save in "subsample-summaries" under destin
obj = readRDS(file = file.path(destin, "summary.RDS"))
create_destin(file.path(folder1, "subsample-summaries"))
saveRDS(obj, file.path(folder1, "subsample-summaries", paste0("summary-", isim, ".RDS")))
