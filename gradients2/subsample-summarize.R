## All directories
outputdir = "/scratch/sangwonh/output"
datadir = "/scratch/sangwonh/data"
scriptdir = "/scratch/sangwonh/scripts/gradients2" ## the current script's location
library(flowmix)

### Read in some command line arguments
parse_args(args = commandArgs(trailingOnly = TRUE), verbose=TRUE)

## Form destination folder
folder1 = paste0("subsample", "-b-", subsample_size)
folder2 = paste0("sim-", isim)
destin = file.path(outputdir, folder1, folder2)

## Run the summarization
cv_summary(destin = destin,
           nfold = 5,
           nrep = 10,
           save = TRUE,
           filename = "summary.RDS")

## Also save in "subsample-summaries" under destin
obj = readRDS(file = file.path(destin, "summary.RDS"))
create_destin(file.path(folder1, "subsample-summaries"))
saveRDS(obj, file = file.path(outputdir, folder1, "subsample-summaries", paste0("summary-", isim, ".RDS")))
