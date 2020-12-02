## Synopsis: form the data.
for(bootrep in 1:20){

  ## Load the original dataset
  load(file = file.path(datadir, paste0("MGL1704-binned-hourly.Rdata")), verbose = TRUE) ## Fix this.

  ## Do whatever processing was done before, from generate-data.R
  X = X[-(1:12),] %>% as.matrix()
  ylist = ylist[-(1:12)]
  ## countslist = qclist[-(1:12)]
  countslist = countslist[-(1:12)]

  ## Then, sample with replacement.


  ## Then, save the data
  save(ylist,
       X,
       countslist,
       file = file.path(datadir, paste0("MGL1704-binned-hourly-bootstrap1-", bootrep, ".Rdata")))

}
