## Synopsis: after this script is run, three objects: |ylist|, |countslist| and
## |X| are ready to use.

#########################################
## QC data at the original time scale  ##
#########################################
if(datatype %in% c(2, 21, 22, 23, 24, 25)){
  ## load(file.path(datadir, "MGL1704-from-CMAP-R-only-qc.Rdata"))
  load(file.path(datadir, "MGL1704-binned.Rdata"))
  X = X[,-which(colnames(X) %in% c("time", "lat", "lon"))]
  X = as.matrix(X)

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
  } else {
    stop("Datatype not recognized!")
  }
}

######################
## Hourly QC data ####
######################
if(datatype %in% c(6, 61, 62, 63, 64, 65, 66)){

  ## Load and name two things and save
  load(file.path(datadir, "MGL1704-hourly-only-binned.Rdata"))

  ## ## Rescale diameters
  ## all.diam = unlist(lapply(ybin_list, function(y) y[,"diam_mid"]))
  ## range.diam = range(all.diam)
  ## width.diam = range.diam[2] - range.diam[1]
  ## all.chl = unlist(lapply(ybin_list, function(y) y[,"chl_small"]))
  ## range.chl = range(all.chl)

  ## Scale the diameter to be in roughly the same range that other data are in.
  ybin_list = lapply(ybin_list, function(y){
    y_diam = y[,"diam_mid"]
    ## y_diam  = y_diam - min(range.diam)
    ## y_diam = y_diam / width.diam
    ## y_diam = y_diam * max(range.chl)
    ## y_diam = cruisedat::rescale_reshift_diam(y_diam, hardcode = TRUE)
    y_diam = rescale_reshift_diam(y_diam, hardcode = TRUE)
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
  if(length(exclude_ind)==0){
    ylist = ybin_list
    countslist = biomass_list
  }
  if(length(exclude_ind)>0){ ## WE MIGHT NOT NEED THIS ANYMORE
  ylist = ybin_list[-exclude_ind]
  countslist = biomass_list[-exclude_ind]
  X = X[-exclude_ind,]
  time = time[-exclude_ind]
  lat = lat[-exclude_ind]
  }
  X = scale(X)
  assert_that(sum(is.na(X))==0)

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

  if(datatype == 64){

    ## ## Reduce to SSS, SST, PAR
    ## X = X[,which(colnames(X) %in% c("sss", "sst", "par"))]

    ## Get rid of redundant cruise SST/SSS measurements
    X = X[,which(!(colnames(X) %in% c("sss_cruise", "sst_cruise")))]

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
    dat = readRDS(file = file.path(datadir, "MGL1704-hourly-paper.RDS"))
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

    oad(file.path(datadir, "MGL1704-hourly-only-binned-1d.Rdata"))
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

    ## Make the data generation a little bit cleaner
    load(file.path(datadir, "MGL1704-hourly-only-binned-1d-diam.Rdata"))
    print("Hourly time resolution with bases and lagged PAR")
    print("Artificial covariates are being added.")

    ## ## Extract time, lat, lon from X
    ## time = X[,"time"]
    ## lat = X[,"lat"]
    ## lon = X[,"lon"]
    ## X = X %>% dplyr::select(-lat, -lon, -time, sss_cruise, sst_cruise) %>% as.matrix

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

    ## Get rid of redundant cruise SST/SSS measurements
    X = X[,which(!(colnames(X) %in% c("sss_cruise", "sst_cruise")))]

    ## Cast ylist's elements as column vectors.
    ylist = lapply(ylist, cbind)

   ##  ## Create the lagged sunlight variable
   ##  lags = c(0,3,6,9,12)
   ##  par = scale(X[,"par"])
   ##  par = par - min(par)
   ##  par = par / max(par)
   ##  parlist = lapply(lags, function(lag)lagpad(par, lag))

   ## ## Make the additional columns
   ## dat = do.call(cbind, parlist)
   ## colnames(dat) = paste0("p", 0:4)

   ## ## Combine them with X
   ## X = X[,-which(colnames(X) == "par")]
   ## X = cbind(X, dat)
   ## print(colnames(X))
   ## X = scale(X)

   ## ## Rid of NA rows.
   ## na.rows = which(apply(dat, 1, function(myrow)any(is.na(myrow))))
   ## if(length(na.rows) > 0){
   ##   X = X[-na.rows,]
   ##   ylist = ylist[-na.rows]
   ##   countslist = countslist[-na.rows]
   ##   ylist = lapply(ylist, cbind)
   ## }
  }
  maxdev = 0.1155245
}

## This is the 1d simulation data for the paper
if(datatype %in% c(80, 81, 82, 83, 84, 85, 86, 87, 88, 89)){

  ## Generate fake 1d data
  obj = generate_data_1d_pseudoreal(bin = FALSE, datadir = datadir,
                                    nt1 = 200,
                                    beta_par = 0.3,
                                    p = 10)
  X = obj$X
  ylist = obj$ylist
  countslist = obj$countslist
  TT = length(ylist)
  print("range of ylist is")
  print(range(unlist(ylist)))

  ## Bin the data (again)
  dat.gridsize = 40
  dat.grid = flowcy::make_grid(ylist, gridsize = dat.gridsize)

  ## 2. Bin with just counts
  obj = flowcy::bin_many_cytograms(ylist, dat.grid, mc.cores = 1, verbose = TRUE)
  ylist = lapply(obj$ybin_list, cbind)
  print("range of binned ylist is")
  countslist = obj$counts_list

  ## maxdev = NULL ##maxdev = 1.5 ## Slightly larger than the truth, which is about (2.66-(-0.84)) * 0.3 = 1.0
}


if(datatype %in% c(9)){

  ## Generate new, realistic 3d data.

  ## Load data.
  ## datadir = "/staging/sh7/sangwonh/data"
  obj = generate_data_1d_pseudoreal_from_cv(datadir = datadir)
  ylist = obj$ylist
  X = obj$X
  countslist = obj$countslist

  ## Bin the data (again)
  dat.gridsize = 40
  dat.grid = flowcy::make_grid(ylist, gridsize = dat.gridsize)

  ## 2. Bin with just counts
  obj = flowcy::bin_many_cytograms(ylist, dat.grid, mc.cores = 1, verbose = TRUE)
  ylist = lapply(obj$ybin_list, cbind)
  countslist = obj$counts_list

  ## pdf("~/Desktop/pseudoreal.pdf", width=15, height=5)
  ## plot_1d(ylist=ylist, countslist=countslist)
  ## graphics.off()

  ## Other setup
  maxdev = 0.1155245
  final_range = lapply(ylist, range) %>% unlist %>% range
  sigma_fac = (final_range[2] - final_range[1])/numclust
  ## numclust = 5
  ## res = covarem(ylist = ylist, X = X, countslist = countslist,
  ##               numclust = 5,
  ##               nrep = 1,
  ##               pie_lambda = .01,
  ##               mean_lambda =.01,
  ##               verbose = TRUE,
  ##               maxdev = 0.1155245,
  ##               sigma_fac = sigma_fac)
  ## plot_1d(ylist=ylist, countslist=countslist, res=res)

  ## Fit data
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
