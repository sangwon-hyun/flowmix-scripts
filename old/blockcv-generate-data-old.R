## Synopsis: after this script is run, three objects: |ylist|, |countslist| and
## |X| are ready to use.

####################
## Simulated data ##
####################

## if(datatype == 0){
if(datatype %in% c(0,100, 1000)){
  if(datatype==0){
    set.seed(0)
    obj = generate_data_generic(p=5, TT=500, fac=0.1, nt=3000, dimdat=3)
  }
  if(datatype==100){
    set.seed(0)
    obj = generate_data_generic(p=5, TT=500, fac=0.05, nt=3000, dimdat=3)
  }
  if(datatype==1000){
    set.seed(0)
    obj = generate_data_generic(p=5, TT=500, fac=0.01, nt=3000, dimdat=3)
  }
  ylist = obj$ylist
  X = obj$X
  countslist = NULL
  print("Simulated 3d particle data with T=500 is being used")



        ## arraynum=1;
        ## arraynum_max=10;
        ## blocktype=1;
        ## datatype=0;
        ## numclust=4;
        ## cv_gridsize=10;
        ## nfold=5;
        ## nrep=5;
        ## random_order=1;
        ## maxres_once=1;

  ######################
  #### Bin this data: ##
  ######################
  ## 1. Make grid:
  dat.gridsize = 40
  dat.grid = make_grid(ylist,
                       gridsize = dat.gridsize)

  ntlist = sapply(ylist, nrow)
  qclist = lapply(1:length(ylist), function(ii)rnorm(ntlist[[ii]], 1, 0.1))

  ## 2. Bin with QC to get biomass
  obj_qc = bin_many_cytograms(ylist, dat.grid, mc.cores = mc.cores,
                              verbose = TRUE, qclist = qclist)
  ybin_list = obj_qc$ybin_list
  biomass_list = obj_qc$counts_list

  ## 3. Assign to static names
  ylist = ybin_list
  countslist = biomass_list
}

#################
## Count data ###
#################
if(datatype %in% c(1,11, 12)){
  TT = 1000
  print(paste0("Count data with T=", TT, " is being used"))
  ## load(file.path(datadir, "MGL1704-from-CMAP-R.Rdata"))
  load(file.path(datadir, "MGL1704-binned.Rdata"))

  ## Clean X
  X = X[, -which(names(X) %in% c("time", "lat","lon"))]

  ## Get rid of the first 342 points
  ylist = ybin_list[-(1:342)]
  countslist = counts_list[-(1:342)]  ## COUNTS are being used
  X = X[-(1:342),]
  X = scale(X)
  assert_that(sum(is.na(X))==0)


  ## Trim to TT time points
  if(datatype %in% c(11,12)){
    ylist = ylist[1:TT]
    countslist = countslist[1:TT]
    X = X[1:TT,]

    ## Use only SSS, SST, PAR
    if(datatype == 12){
      print("Only sss, sst, par are being used.")
      X = X[,which(colnames(X) %in% c("sss", "sst", "par", "sss_cruise", "sst_cruise"))]
    }
  }
}


#############
## QC data ##
#############
if(datatype %in% c(2, 21, 22, 23, 24, 25)){
  print("loading QC data")
  ## load(file.path(datadir, "MGL1704-from-CMAP-R-only-qc.Rdata"))
  load(file.path(datadir, "MGL1704-binned.Rdata"))
  X = X[,-which(colnames(X) %in% c("time", "lat", "lon"))]
  X = as.matrix(X)
  ## qclist = lapply(ylist, function(y)y[,"Qc_mid"])
  ## ylist = lapply(ylist, function(y) y[,c("fsc_small", "chl_small", "pe")])

  print("Biomass data is being used.")


  ## Get rid of the first 342 points
  ylist = ybin_list[-(1:342)]
  countslist = biomass_list[-(1:342)] ## BIOMASS is being used
  X = X[-(1:342),]
  X = scale(X)
  assert_that(sum(is.na(X))==0)

  ## Trim to TT time points
  if(datatype %in% c(21, 22)){
    TT = 500
    X = X[1:TT,]
    ylist = ylist[1:TT]
    countslist = countslist[1:TT]

    ## Use only SSS, SST, PAR
    if(datatype == 22){
      print("Only sss, sst, par are being used.")
      X = X[,which(colnames(X) %in% c("sss", "sst", "par", "sss_cruise", "sst_cruise"))]
    }
  } else if (datatype == 23){

    TT = 200
    set.seed(0)
    times = sample(1:(nrow(X)), TT)
    set.seed(Sys.time())

    X = X[times,]
    ylist = ylist[times]
    countslist = countslist[times]

  } else if (datatype == 24){
      TT = 500
      X = X[500+(1:TT),]
      ylist = ylist[500+(1:TT)]
      countslist = countslist[500+(1:TT)]
  } else if (datatype == 25){
      TT = 500
      X = X[1000+(1:TT),]
      ylist = ylist[1000+(1:TT)]
      countslist = countslist[1000+(1:TT)]
  } else {
    stop("Datatype not recognized!")
  }

}



##########################
## Collapsed count data ##
##########################
if(datatype %in% 3){
  load(file.path(datadir, "MGL1704-from-CMAP-R-collapsed.Rdata"))

  ## Clean X
  X = X_collapsed
  X = X[,-which(colnames(X) %in% c("time", "lat", "lon"))]
  X = as.matrix(X)
  X = scale(X)

  print("COLLAPSED count data is being used.")

  ## Obtain ylist
  ## ybin_list = ybin_collapsed_list
  ylist = lapply(ybin_list_collapsed, function(a)a[,1:3])
  countslist = counts_list_collapsed

  ## Get rid of the first 342 points
  na.ind = which(is.na(X[,"par"]))
  ylist = ylist[-na.ind]
  countslist = countslist[-na.ind]
  X = X[-na.ind,]
  assert_that(sum(is.na(X))==0)
}


#########################
## Collapsed QC data ####
#########################
if(datatype %in% 4){
  print("loading QC data")
  load(file.path(datadir, "MGL1704-from-CMAP-R-collapsed.Rdata"))

  X = X_collapsed
  X = X[,-which(colnames(X) %in% c("time", "lat", "lon"))]
  X = as.matrix(X)
  X = scale(X)

  print("COLLAPSED Biomass (QC) data is being used.")

  ## Assign it to static names ylist and countslist
  ylist = lapply(ybin_list_collapsed, function(a)a[,1:3])
  countslist = biomass_list_collapsed

  ## Get rid of the first 342 points
  na.ind = which(is.na(X[,"par"]))
  ylist = ylist[-na.ind]
  countslist = countslist[-na.ind]
  X = X[-na.ind,]
  assert_that(sum(is.na(X))==0)
}


#########################
## Collapsed QC data ####
#########################
if(datatype %in% 5){

  ## Randomly picking time points
  TT = 500
  print("loading QC data")
  ## load(file.path(datadir, "MGL1704-from-CMAP-R-only-qc.Rdata"))
  load(file.path(datadir, "MGL1704-binned.Rdata"))
  X = X[,-which(colnames(X) %in% c("time", "lat", "lon"))]
  X = as.matrix(X)
  ## qclist = lapply(ylist, function(y)y[,"Qc_mid"])
  ## ylist = lapply(ylist, function(y) y[,c("fsc_small", "chl_small", "pe")])

  print("Biomass data is being used.")


  ## Get rid of the first 342 points
  ylist = ybin_list[-(1:342)]
  countslist = biomass_list[-(1:342)] ## BIOMASS is being used
  X = X[-(1:342),]
  X = scale(X)
  assert_that(sum(is.na(X))==0)

  ## Trim to TT time points
  if(datatype %in% c(21,22)){
    X = X[1:TT,]
    ylist = ylist[1:TT]
    countslist = countslist[1:TT]

    ## Use only SSS, SST, PAR
    if(datatype == 22){
      print("Only sss, sst, par are being used.")
      X = X[,which(colnames(X) %in% c("sss", "sst", "par", "sss_cruise", "sst_cruise"))]
    }
  }


}

######################
## Hourly QC data ####
######################
if(datatype %in% c(6, 61, 62, 63, 64, 65, 66)){

  ## Load and name two things and save
  load(file.path(datadir, "MGL1704-hourly-only-binned.Rdata"))
  ## load(file.path(datadir, "MGL1704-hourly.Rdata"))
  ## time = X[,"time"]
  ## names(counts_list) = time
  ## names(biomass_list) = time
  ## save(X, ylist, ybin_list, counts_list, biomass_list, file=file.path(datadir,
  #                                                                     "MGL1704-hourly.Rdata"))
  ## End of that

  ## Rescale diameters
  all.diam = unlist(lapply(ybin_list, function(y) y[,"diam_mid"]))
  range.diam = range(all.diam)
  width.diam = range.diam[2] - range.diam[1]
  all.chl = unlist(lapply(ybin_list, function(y) y[,"chl_small"]))
  range.chl = range(all.chl)
  all.pe = unlist(lapply(ybin_list, function(y) y[,"pe"]))
  range.pe = range(all.pe)


  ## Scale the diameter to be in roughly the same range that this is in.
  ybin_list = lapply(ybin_list, function(y){
    y_diam = y[,"diam_mid"]
    y_diam  = y_diam - min(range.diam)
    y_diam = y_diam / width.diam
    y_diam = y_diam * max(range.chl)
    y[,"diam_mid"] = y_diam
    return(y)
  })

  print("Biomass data is being used.")
  ## load(file.path(datadir, "MGL1704-hourly.Rdata"))
  time = X[,"time"]
  lat = X[,"lat"]
  X = X[,-which(colnames(X) %in% c("time", "lat", "lon"))]
  X = as.matrix(X)

  ## Remove unused large things
  rm(ylist)

  ## Get rid of the first 342 points
  exclude_ind = which(is.na(X[,"par"]))
  ylist = ybin_list[-exclude_ind]
  countslist = biomass_list[-exclude_ind]
  X = X[-exclude_ind,]
  X = scale(X)
  assert_that(sum(is.na(X))==0)
  time = time[-exclude_ind]
  lat = lat[-exclude_ind]

  ## ## Remove the date columns from the
  ## for(datecolumn in c("yyyy", "mm", "yy")){
  ##   if(datecolumn %in% colnames(X)){
  ##     X = X[,-which(colnames(X) == datecolumn)]
  ##   }
  ## }

  ## |ylist|, |X|, and |countslist| are ready.
  if(datatype == 61){
    print("Only sss, sst, par are being used.")
    X = X[,which(colnames(X) %in% c("sss", "sst", "par", "sss_cruise", "sst_cruise"))]
  }

  if(datatype == 62){
    ## Add the artificial covariates.
    print("Artificial covariates are being added.")
    X = add_growth(X)
    X = add_div(X)
    X = add_sine_to_X(X, time)
  }

  if(datatype == 63){
    print("Only sss, sst, par are being used.")
    X = X[,which(colnames(X) %in% c("sss", "sst", "par", "sss_cruise", "sst_cruise"))]
    print("Artificial covariates are being added.")
    X = add_growth(X)
    X = add_div(X)
    X = add_sine_to_X(X, time)
  }


  if(datatype == 64){

    ## Reduce to SSS, SST, PAR
    X = X[,which(colnames(X) %in% c("sss", "sst", "par"))]

    ## Add lagged PAR variable
    lags = c(3,6,9,12)
    X = add_lagpar(X, lags)
    X = scale(X)

    ## Add transition crossing
    X = add_transition(X, lat)

    ## Get rid of the first 12 points that have NAs due to lagging PAR
    na.rows = which(apply(X, 1, function(myrow) any(is.na(myrow))))
    stopifnot(all(na.rows==1:12))
    if(length(na.rows) > 0){
      X = X[-na.rows,]
      ylist = ylist[-na.rows]
      countslist = countslist[-na.rows]
    }

    ## ## Temporary: see the dimensions, for a sanity check.
    ## print(dim(X))
    ## print(sapply(ylist, dim))
    ## print(sapply(countslist, length))
    ## save(X, ylist, countslist, file=file.path(outputdir, "temp.Rdata"))
  }

  if(datatype == 65){

    ## Reduce to SSS, SST, PAR
    X = X[,which(colnames(X) %in% c("sss", "sst", "par"))]

    ## Add lagged PAR variable
    lags = c(3,6,9,12)
    X = add_lagpar(X, lags)
    X = scale(X)

    ## Add transition crossing
    X = add_transition(X, lat)

    ## Add interactions
    mat = model.matrix( ~ sss + sst + (b1 + b2)*(par + p1 + p2 + p3 + p4) - 1, data.frame(X))
    mat = mat[,-which(colnames(mat)%in% c("(Intercept)"))]

    ## Get rid of the first 12 points that have NAs due to lagging PAR
    na.rows = which(apply(X, 1, function(myrow) any(is.na(myrow))))
    stopifnot(all(na.rows==1:12))
    if(length(na.rows) > 0){
      X = X[-na.rows,]
      ylist = ylist[-na.rows]
      countslist = countslist[-na.rows]
    }
  }

  if(datatype == 66){

    ## Reduce to SSS, SST, PAR
    X = X[,which(colnames(X) %in% c("sss", "sst", "par"))]

    ## Add lagged PAR variable
    lags = c(3,6,9,12)
    X = add_lagpar(X, lags)
    X = scale(X)

    ## Isolate attention to the latter third-ish of the data
    ## X = add_transition(X, lat)
    north = range(which(lat > 37))
    third_half = (north[2]:nrow(X))
    X = X[third_half,]
    ylist = ylist[third_half]
    countslist = countslist[third_half]

    print(head(X))
    print(class(X))
  }
}


###################
## 1d real data ###
###################
if(datatype %in% c(7, 71, 72, 73, 74, 75, 76)){
  if(datatype==7){
    ## Load data from static names |X|, |ylist|, |countslist|
    load(file.path(datadir, "MGL1704-binned-1d.Rdata"))

    ## Add some artificial covariates
    print("Artificial covariates are being added.")
    X = add_growth(X)
    X = add_div(X)
    ## X = add_sine_to_X(X, time)

    ## plot(X[,"par"], type='l')
    ## lines(X[,"growth"]/max(X[,"growth"]), type='l', col='grey30', lwd=3, lty=2)
    ## lines(X[,"div"]/max(abs(X[,"div"])), type='l', col='grey30', lwd=3, lty=3)
    X = scale(X)
    ylist = lapply(ylist, cbind)
  }
  if(datatype==71){
    load(file.path(datadir, "MGL1704-hourly-only-binned-1d.Rdata"))

    ## Add some artificial covariates
    print("Hourly time resolution")
    print("Artificial covariates are being added.")
    X = add_growth(X)
    X = add_div(X)
    ## X = add_sine_to_X(X, time)

    ## plot(X[,"par"], type='l')
    ## lines(X[,"growth"]/max(X[,"growth"]), type='l', col='grey30', lwd=3, lty=2)
    ## lines(X[,"div"]/max(abs(X[,"div"])), type='l', col='grey30', lwd=3, lty=3)
    X = scale(X)
    ylist = lapply(ylist, cbind)
  }

  if(datatype == 72){
    load(file.path(datadir, "MGL1704-hourly-only-binned-1d.Rdata"))
    print("Hourly time resolution")
    print("Only sss, sst, par are being used.")
    X = X[,which(colnames(X) %in% c("sss", "sst", "par", "sss_cruise", "sst_cruise"))]
    print("Artificial covariates are being added.")
    X = add_growth(X)
    X = add_div(X)
    X = scale(X)
    ylist = lapply(ylist, cbind)
  }

  if(datatype == 73){

    load(file.path(datadir, "MGL1704-hourly-only-binned-1d.Rdata"))
    print("Hourly time resolution with bases and lagged PAR")
    print("Only sss, sst, par are being used.")
    X = X[,which(colnames(X) %in% c("sss", "sst", "par", "sss_cruise", "sst_cruise"))]
    print("Artificial covariates are being added.")

     ## Create the lagged sunlight variable
    lags = c(0,3,6,9,12)
    par = scale(X[,"par"])
    par = par - min(par)
    par = par / max(par)
    parlist = lapply(lags, function(lag)lagpad(par, lag))
    dat = do.call(cbind, parlist)
    colnames(dat) = paste0("p", 1:5)
    dat = as.data.frame(dat)

    ## Combine them with X
    X = X[,-which(colnames(X) == "par")]
    X = cbind(X, dat)
    print(colnames(X))

    ## Rid of NA rows.
    na.rows = which(apply(dat, 1, function(myrow)any(is.na(myrow))))
    if(length(na.rows) > 0){
      X = X[-na.rows,]
      ylist = ylist[-na.rows]
      countslist = countslist[-na.rows]
    }
    X = scale(X)
    ylist = lapply(ylist, cbind)
  }

  if(datatype == 74){

    load(file.path(datadir, "MGL1704-hourly-only-binned.Rdata"))
    lat = X[,"lat"]

    load(file.path(datadir, "MGL1704-hourly-only-binned-1d.Rdata"))
    print("Hourly time resolution with bases and lagged PAR")

    print("Only sss, sst, par are being used.")
    X = X[,which(colnames(X) %in% c("sss", "sst", "par", "sss_cruise", "sst_cruise"))]
    print("Artificial covariates are being added.")

     ## Create the lagged sunlight variable
    lags = c(0,3,6,9,12)
    par = scale(X[,"par"])
    par = par - min(par)
    par = par / max(par)
    parlist = lapply(lags, function(lag)lagpad(par, lag))


    ## Create the TF bases.
    TT = length(ylist)
    ind = rep(0, TT)
    north = which(lat > 37)
    regions = sapply(1:3, function(ii){
      endpt = c(0,range(north),TT)[ii:(ii+1)]
      (endpt[1]+1):endpt[2]
    })
    bases0 = lapply(regions, function(reg){
      vec = rep(0,TT)
      vec[reg] = 1
      vec
    })
    bases1 = lapply(regions, function(reg){
      vec = rep(0,TT)
      vec[reg] = (1:length(reg))/length(reg)
      vec
    })
    bases = c(bases0, bases1)

    ## Make the additional columns
    dat = cbind(do.call(cbind, bases),
                do.call(cbind, parlist))
    colnames(dat) = c(paste0("b", 1:6), paste0("p", 1:5))
    dat = as.data.frame(dat)


    ## Make model matrix
    mat = model.matrix( ~ b4 + b6 + b5 + (b1 + b2 + b3)*( p1 + p2 + p3 + p4 + p5),dat)
    mat = mat[,-which(colnames(mat)%in% c("(Intercept)"))]

    ## Rid of NA rows.
    na.rows = which(apply(dat, 1, function(myrow)any(is.na(myrow))))
    if(length(na.rows) > 0){
      X = X[-na.rows,]
      ylist = ylist[-na.rows]
      countslist = countslist[-na.rows]
    }

    ## Combine them with X
    X = X[,-which(colnames(X) == "par")]
    X = cbind(X, mat)
    print(colnames(X))
  }

  if(datatype == 75){

    ## Load the hourly level X
    load(file.path(datadir, "MGL1704-hourly-only-binned.Rdata"))
    exclude_ind = which(is.na(X[,"par"]))
    X = X[-exclude_ind,]
    lat = X[,"lat"]

    ## Create the TF bases.
    TT = nrow(X)
    ind = rep(0, TT)
    north = which(lat > 37)
    regions = sapply(1:3, function(ii){
      endpt = c(0,range(north),TT)[ii:(ii+1)]
      (endpt[1]+1):endpt[2]
    })

   load(file.path(datadir, "MGL1704-hourly-only-binned-1d-diam.Rdata"))
   print("Hourly time resolution with bases and lagged PAR")
   print("Only sss, sst, par are being used.")
   X = X[,which(colnames(X) %in% c("sss", "sst", "par", "sss_cruise", "sst_cruise"))]
   print("Artificial covariates are being added.")

    ## Create the lagged sunlight variable
   lags = c(0,3,6,9,12)
   par = scale(X[,"par"])
   par = par - min(par)
   par = par / max(par)
   parlist = lapply(lags, function(lag)lagpad(par, lag))

   ## Make the additional columns
   dat = do.call(cbind, parlist)
   colnames(dat) = paste0("p", 0:4)

   ## Combine them with X
   X = X[,-which(colnames(X) == "par")]
   X = cbind(X, dat)
   print(colnames(X))
   X = scale(X)

   ## Add bases
   bases0 = lapply(regions, function(reg){
     vec = rep(0, TT)
     vec[reg] = 1
     vec
   })
   bases = do.call(cbind, bases0)
   ## bases = do.call(cbind, bases0[2:3]) ## Omitting one basis.
   ## colnames(bases) = c("b2", "b3")
   colnames(bases) = c("b1", "b2", "b3")
   X = cbind(X, bases)

   ## Rid of NA rows.
   na.rows = which(apply(dat, 1, function(myrow)any(is.na(myrow))))
   if(length(na.rows) > 0){
     X = X[-na.rows,]
     ylist = ylist[-na.rows]
     countslist = countslist[-na.rows]
     ylist = lapply(ylist, cbind)
   }
  }
  if(datatype==76){



    ## Load the hourly level X, for the only intent of obtaining latitude
    load(file.path(datadir, "MGL1704-hourly-only-binned.Rdata"))
    exclude_ind = which(is.na(X[,"par"]))
    X = X[-exclude_ind,]
    lat = X[,"lat"]

    ## Identify third section of the cruise.
    TT = nrow(X)
    north = which(lat > 37)
    third = (max(north)):TT
    ## regions = sapply(1:3, function(ii){
    ##   endpt = c(0,range(north),TT)[ii:(ii+1)]
    ##   (endpt[1]+1):endpt[2]
    ## })

   load(file.path(datadir, "MGL1704-hourly-only-binned-1d-diam.Rdata"))
   print("Hourly time resolution with bases and lagged PAR")
   print("Only sss, sst, par are being used.")
   X = X[,which(colnames(X) %in% c("sss", "sst", "par", "sss_cruise", "sst_cruise"))]
   print("Artificial covariates are being added.")

    ## Isolate attention to third section of cruise
    X = X[third,]
    ylist = ylist[third]
    countslist = countslist[third]

    ## Create the lagged sunlight variable
   lags = c(0,3,6,9,12)
   par = scale(X[,"par"])
   par = par - min(par)
   par = par / max(par)
   parlist = lapply(lags, function(lag)lagpad(par, lag))

   ## Make the additional columns
   dat = do.call(cbind, parlist)
   colnames(dat) = paste0("p", 0:4)

   ## Combine them with X
   X = X[,-which(colnames(X) == "par")]
   X = cbind(X, dat)
   print(colnames(X))
   X = scale(X)

   ## Rid of NA rows.
   na.rows = which(apply(dat, 1, function(myrow)any(is.na(myrow))))
   if(length(na.rows) > 0){
     X = X[-na.rows,]
     ylist = ylist[-na.rows]
     countslist = countslist[-na.rows]
     ylist = lapply(ylist, cbind)
   }
  }
  maxdev = 0.1155245
}



## This is the simulation
if(datatype %in% c(8, 81, 82)){


  ## Sample settings
  arraynum_max=40;
  blocktype=1;
  datatype=8;
  numclust=2;
  cv_gridsize=7;
  nfold=5;
  nrep=5;
  random_order=0;
  maxres_once=1;
  fitonly=0;
  high_range_beta=1;
  sim=1;
  nsim=100;
  add_sigma_numer=0;
  add_sigma_denom=1;
  ## Sample settings

  ## Add data
  ## Generate fake, 1-mixture data.
  set.seed(0)
  TT = 100
  ntlist = c(rep(800, TT/2), rep(1000, TT/2))

  ## Generate covariate
  load(file.path(datadir, "MGL1704-hourly-only-binned.Rdata"))
  par = X[, "par"]
  par = par[!is.na(par)]
  par = ksmooth(x=1:length(par), y=par, bandwidth=5, x.points = 1:length(par))$y
  X = cbind( scale(par[1:TT]), rnorm(TT), c(rep(0, TT/2), rep(1, TT/2)))
  colnames(X) = c("par", "noise", "cp")



  ## Means over time
  p = 3
  numclust = 2
  beta = matrix(0, ncol=numclust, nrow=p+1)

  ## Beta coefficients
  beta[0+1,1] = 0
  beta[1+1,1] = 3
  beta[2+1,1] = 0
  beta[0+1,2] = 10
  beta[1+1,2] = -3
  beta[2+1,2] = 0
  ## beta[0+1,1] = 0
  ## beta[1+1,1] = 5
  ## beta[2+1,1] = 0
  ## beta[0+1,2] = 10
  ## beta[1+1,2] = -10
  ## beta[2+1,2] = 0

  ## (n x (p+1)) * ((p+1) x numclust)
  mnmat = cbind(1, X) %*% beta
  ## matplot(mnmat, type='l')

  ## Alpha coefficients
  ## Pies over time
  alpha = matrix(0, ncol=numclust, nrow=p+1)

  alpha[0+1, 1] = 0
  alpha[1+1, 1] = 0
  alpha[2+1, 1] = 0
  alpha[3+1, 1] = 0

  alpha[0+1, 2] = -10
  alpha[1+1, 2] = 0
  alpha[2+1, 2] = 0
  alpha[3+1, 2] = 10 + log(1/4)
  colnames(alpha) = paste0("clust", 1:numclust)
  rownames(alpha) = c("intercept", "par", "noise", "cp")
  pie = exp(cbind(1,X) %*% alpha)
  pie = pie/rowSums(pie)
  ## print(round(pie,3))
  ## matplot(pie, type='l', lwd=3)

  ## Samples nt memberships out of 1:numclust according to the probs in pie.
  ## Data is a probabilistic mixture from these two means, over time.
  ylist = lapply(1:TT,
                 function(tt){
                   ## print(round(pie[[tt]],3))
                   draws = sample(1:numclust,
                                  size = ntlist[tt], replace = TRUE,
                                  prob = c(pie[[tt]], 1-pie[[tt]]))
                   mns = mnmat[tt,]
                   means = mns[draws]
                   datapoints = means + rnorm(ntlist[tt], 0, 1)
                   cbind(datapoints)
                 })
  countslist=NULL
  maxdev = 6
}



if(datatype %in% c(9)){

  ## Generate new, realistic 3d data.

  ## Two clusters, one is constantly present, and one appears.

  ## Run cross-validation to get the problem

}










## A quick fix: there are two chlorophyll columns. One is bad. Get rid of that
## one.

## matplot(X[,c("CHL", "chl")], type='l', col=c('black', 'blue'), lwd=2)
if("chl" %in% colnames(X)){
  X = X[, -which(colnames(X) == "chl")]
  assertthat::assert_that(!("chl" %in% colnames(X) ))
}


## Another quick fix
## matplot(X[,c("sea_water_temp_WOA_clim", "sst")], col=c("black", "blue"), lwd=2, lty=1, type='l')
if("sea_water_temp_WOA_clim" %in% colnames(X)){
  X = X[, -which(colnames(X) == "sea_water_temp_WOA_clim")]
  assertthat::assert_that(!("sea_water_temp_WOA_clim" %in% colnames(X) ))
}
