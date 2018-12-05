# ================================================
#
#  log sequence
# 

lseq <- function(from, to, length = 2, base = exp(1)) {
  .lseq <- base^(seq(log(from, base = base),log(to, base = base), length.out = length))
  return(.lseq)
}
l10seq <- function(from, to, length = 2, base = 10) {
  .lseq <- base^(seq(log(from, base = base),log(to, base = base), length.out = length))
  return(.lseq)
}


# ============================= for pairs plot
gatherpairs <- function(data, ..., 
                        xkey = '.xkey', xvalue = '.xvalue',
                        ykey = '.ykey', yvalue = '.yvalue',
                        na.rm = FALSE, convert = FALSE, factor_key = FALSE) {
  vars <- quos(...)
  xkey <- enquo(xkey)
  xvalue <- enquo(xvalue)
  ykey <- enquo(ykey)
  yvalue <- enquo(yvalue)
  
  data %>% {
    cbind(gather(., key = !!xkey, value = !!xvalue, !!!vars,
                 na.rm = na.rm, convert = convert, factor_key = factor_key),
          select(., !!!vars)) 
  } %>% gather(., key = !!ykey, value = !!yvalue, !!!vars,
               na.rm = na.rm, convert = convert, factor_key = factor_key)
}

#==============================================================================


psurv_predict <- function(MCMC_data,
                          newdata = NULL,
                          draws = NULL,
                          interpolate_length = 1e2,
                          interpolate_method = "linear",
                          mc.cores = 1){
  
  if(is.null(newdata)){
    newdata <- x$data
  }
  
  filter_newdata <- newdata %>%
    dplyr::filter(!is.na(conc))
  
  filter_newdata_list <- split(filter_newdata, filter_newdata$replicate)
  
  newdata_list <- split(newdata, newdata$replicate)
  
  if(length(newdata_list) != length(filter_newdata_list)){
    stop("A replicate has all 'conc' with 'NA'")
  }
  
  if(mc.cores == 1){
    cat('progress:', 0, '% \n')
  } else{
    cat('For information: No progress indicator is available when mc.cores > 1. \n')
  }
  
  n_replicate <- length(newdata_list)
  
  pp_list_solveTKTD <- parallel::mclapply(1:n_replicate,
                                          function(it) {
                                            ode_TKTD = solveTKTD(x,
                                                                 time_conc = filter_newdata_list[[it]]$time,
                                                                 conc = filter_newdata_list[[it]]$conc,
                                                                 draws = draws,
                                                                 interpolate_time = newdata_list[[it]]$time,
                                                                 interpolate_length = interpolate_length,
                                                                 interpolate_method = interpolate_method)
                                            ## progress
                                            percent_draws = it/n_replicate * 100
                                            cat('progress:', percent_draws, '% \n')
                                            ##
                                            return(ode_TKTD)                      
                                          },
                                          mc.cores = mc.cores)
  
  names(pp_list_solveTKTD) <- names(newdata_list)
  
  out_posterior_predict <- dplyr::bind_rows(pp_list_solveTKTD, .id = "replicate")
  
  return(out_posterior_predict)
}

# SD model solver
model_SD <- function(t, State, parms, input)  {
  with(as.list(c(parms, State)), {
    
    C = rep(input(t), length = draws)
    
    D = State[1:draws]
    H = State[(draws+1):(2*draws)] 
    
    dD <- kd*(C - D)     # internal concentration
    dH <- kk*pmax(0, D-z) + hb # risk function
    
    res <- c(dD, dH)
    
    list(res, signal = input(t))
    
  })
}

# IT model solver
model_IT <- function(t, State, parms, input) {
  with(as.list(c(parms, State)), {
    
    C = rep(input(t), length = draws)
    
    D = State[1:draws]
    
    dD <- kd*(C - D)    # internal concentration
    
    list(dD = dD, signal = input(t))
    
  })
}

# solver TKTD with deSolve
solveTKTD_svt <- function(MCMC_param,
                      model_type = NULL,
                      time_conc = NULL,
                      conc = NULL,
                      draws = NULL,
                      interpolate_time = NULL,
                      interpolate_length = 1e2,
                      interpolate_method = "linear"){
  
  MCMC_stanEstim <- MCMC_param
  
  if(is.null(draws)){
    draws <- nrow(MCMC_stanEstim)
  } else{
    if(draws > nrow(MCMC_stanEstim)) { stop('draws > size of MCMC.')}
    seq_MCMC <- round(seq(from = 1, to = nrow(MCMC_stanEstim), length.out = draws)) # with draws < size MCMC, elements of seq_MCMC are unique ! 
    MCMC_stanEstim <- MCMC_stanEstim[ seq_MCMC , ]
  }
  
  ## external signal with several rectangle impulses
  sigimp <- approxfun(time_conc, conc, method = interpolate_method, rule = 2)
  
  if(is.null(interpolate_time)){
    times <- seq(min(time_conc), max(time_conc), length.out = interpolate_length)
  } else{
    times <-  sort(unique(c(seq(min(interpolate_time), max(interpolate_time), length.out = interpolate_length), interpolate_time)))
  }
  
  ## model
  kd <- MCMC_stanEstim$kd
  hb <- MCMC_stanEstim$hb
  
  if(model_type == "IT"){
    
    alpha <- MCMC_stanEstim$alpha
    beta <- MCMC_stanEstim$beta
    
    ## The parameters
    parms  <- list( kd = kd,
                    alpha = alpha,
                    beta = beta,
                    draws = draws)
    
    ## Start values for steady state
    xstart <- c(D = rep(0,draws))
    
    ## Solve model
    out <- ode(y = xstart,
               times = times,
               func = model_IT,
               parms,
               input = sigimp)
    
  }
  if (model_type == "SD"){
    
    kk = MCMC_stanEstim$kk
    z = MCMC_stanEstim$z
    
    ## The parameters
    parms  <- list( kd = kd,
                    hb = hb,
                    kk = kk,
                    z = z,
                    draws = draws)
    
    ## Start values for steady state
    xstart <- c(D = rep(0,draws),
                H = rep(0,draws))
    
    ## Solve model
    out <- ode(y = xstart,
               times = times,
               func = model_SD,
               parms,
               input = sigimp)
  }
  
  if(model_type == "IT"){
    D_matrix = as.data.frame(out) %>%
      dplyr::select(-c(time, signal))%>%
      as.matrix() %>%
      matrixStats::colCummaxs()
    #S <- exp( - hb %*% t(times)) * ( 1-plogis(log(t(D_matrix)), location = log(parms$alpha), scale = 1/parms$beta))
    #S <- exp( - hb %*% t(times)) * ( 1-1 / (1+ (t(D_matrix)/parms$alpha)^(-parms$beta)))
    #S <- exp( - hb %*% t(times)) * ( 1-1 / (1+ exp((t(D_matrix) - parms$alpha)/(parms$beta))))
    #S <- exp( - hb %*% t(times)) * ( 1- t(D_matrix)^(parms$beta) / ( (parms$alpha)^(parms$beta) + t(D_matrix)^(parms$beta) ))
    S <- exp( - hb %*% t(times)) * ( 1 - 1 / (1+ exp( -parms$beta * (log(t(D_matrix)) - log(parms$alpha) )))) 
    dtheo <- t(S)
    
  }
  if(model_type=="SD"){
    
    H_matrix = as.data.frame(out) %>%
      dplyr::select(contains("H"),-c(time, signal))%>%
      as.matrix()
    S <- exp(-H_matrix)
    dtheo <- S
  }
  
  # OUTPUT --------------------------------------------------------------------
  
  pp_solveTKTD <- as_data_frame(dtheo) %>%
    mutate(time = out[, "time"], 
           conc = out[, "signal"])
  
  return(pp_solveTKTD)
}


# solver TKTD with deSolve
solveTKTD_MFx <- function(x,
                          model_type = NULL,
                          time_conc = NULL,
                          conc = NULL,
                          draws = NULL,
                          interpolate_time = NULL,
                          interpolate_length = 1e2,
                          interpolate_method = "linear"){
  
  # parameters
  mctot <- do.call("rbind", x$mcmc)

  MCMC_stanEstim <- data.frame(kd = 10^mctot[, "kd_log10"],
                               hb = 10^mctot[, "hb_log10"])
  
  if(model_type == "SD"){
    MCMC_stanEstim$z <- 10^mctot[, "z_log10"]
    MCMC_stanEstim$kk <- 10^mctot[, "kk_log10"]
  }
  if(model_type == "IT"){
    MCMC_stanEstim$alpha <- 10^mctot[, "alpha_log10"]
    MCMC_stanEstim$beta <- 10^mctot[, "beta_log10"]
  }
  
  if(is.null(draws)){
    draws <- nrow(MCMC_stanEstim)
  } else{
    if(draws > nrow(MCMC_stanEstim)) { stop('draws > size of MCMC.')}
    seq_MCMC <- round(seq(from = 1, to = nrow(MCMC_stanEstim), length.out = draws)) # with draws < size MCMC, elements of seq_MCMC are unique ! 
    MCMC_stanEstim <- MCMC_stanEstim[ seq_MCMC , ]
  }
  
  ## external signal with several rectangle impulses
  sigimp <- approxfun(time_conc, conc, method = interpolate_method, rule = 2)
  
  if(is.null(interpolate_time)){
    times <- seq(min(time_conc), max(time_conc), length.out = interpolate_length)
  } else{
    times <-  sort(unique(c(seq(min(interpolate_time), max(interpolate_time), length.out = interpolate_length), interpolate_time)))
  }
  
  ## model
  kd <- MCMC_stanEstim$kd
  hb <- MCMC_stanEstim$hb
  
  if(model_type == "IT"){
    
    alpha <- MCMC_stanEstim$alpha
    beta <- MCMC_stanEstim$beta
    
    ## The parameters
    parms  <- list( kd = kd,
                    alpha = alpha,
                    beta = beta,
                    draws = draws)
    
    ## Start values for steady state
    xstart <- c(D = rep(0,draws))
    
    ## Solve model
    out <- ode(y = xstart,
               times = times,
               func = model_IT,
               parms,
               input = sigimp)
    
  }
  if (model_type == "SD"){
    
    kk = MCMC_stanEstim$kk
    z = MCMC_stanEstim$z
    
    ## The parameters
    parms  <- list( kd = kd,
                    hb = hb,
                    kk = kk,
                    z = z,
                    draws = draws)
    
    ## Start values for steady state
    xstart <- c(D = rep(0,draws),
                H = rep(0,draws))
    
    ## Solve model
    out <- ode(y = xstart,
               times = times,
               func = model_SD,
               parms,
               input = sigimp)
  }
  
  if(model_type == "IT"){
    
    D_dataframe = as.data.frame(out) %>%
      dplyr::select(-c(time, signal))
    
    D_matrix = as.data.frame(out) %>%
      dplyr::select(-c(time, signal))%>%
      as.matrix() %>%
      matrixStats::colCummaxs()
    #S <- exp( - hb %*% t(times)) * ( 1-plogis(log(t(D_matrix)), location = log(parms$alpha), scale = 1/parms$beta))
    #S <- exp( - hb %*% t(times)) * ( 1-1 / (1+ (t(D_matrix)/parms$alpha)^(-parms$beta)))
    #S <- exp( - hb %*% t(times)) * ( 1-1 / (1+ exp((t(D_matrix) - parms$alpha)/(parms$beta))))
    #S <- exp( - hb %*% t(times)) * ( 1- t(D_matrix)^(parms$beta) / ( (parms$alpha)^(parms$beta) + t(D_matrix)^(parms$beta) ))
    S <- exp( - hb %*% t(times)) * ( 1 - 1 / (1+ exp( -parms$beta * (log(t(D_matrix)) - log(parms$alpha) )))) 
    dtheo <- t(S)
    
  }
  if(model_type=="SD"){
    
    D_dataframe = as.data.frame(out) %>%
      dplyr::select(contains("D"), -c(time, signal))

    H_matrix = as.data.frame(out) %>%
      dplyr::select(contains("H"),-c(time, signal))%>%
      as.matrix()
    S <- exp(-H_matrix)
    dtheo <- S
  }
  
  # OUTPUT --------------------------------------------------------------------
  
  df_internalConc <- D_dataframe %>%
    mutate(time = out[, "time"], 
           conc = out[, "signal"])
    
  
  df_survival <- as_data_frame(dtheo) %>%
    mutate(time = out[, "time"], 
           conc = out[, "signal"])
  
  return(list(df_internalConc = df_internalConc,
              df_survival = df_survival))
}

# -----------------------------------------------------------------------------
# LCx

LCx_along_X <- function(object, time_LCx = time_LCx, conc_range = NULL, npoints = 1e2, length = 1e2){
  
  if(is.null(conc_range)){
    conc_range = seq(0, max(object$jags.data$conc), length.out = npoints)
  } else{
    if(length(conc_range) != 2){
      stop('conc_range must a vector of length 2 with minimal and maximal value of the range of concentration')
    }
    conc_range = seq(conc_range[1], conc_range[2], length.out = npoints)
  }
  
  if(min(conc_range) != 0){
    stop("Minimal value of 'conc_range' must be 0.")
  }
  
  if(is.null(time_LCx)){
    time_LCx = max(object$jags.data$time)
  }
  
  df_dose <- doseResponse_survFitCstExp(x = object, time_LCx = time_LCx, conc_range, npoints)
  
  median_backgroundMortality_Conc0 = dplyr::filter(df_dose, concentration == 0)$q50
  
  X_prop_provided = seq(0,100, length.out = length)
  
  ls_LCx = list()
  
  for(i in 1:length(X_prop_provided)){
      X_prop <- (100-X_prop_provided[i])/100*median_backgroundMortality_Conc0
      ls_LCx[[i]] <- pointsLCx(df_dose, X_prop)
  }
  
  df_LCx <- dplyr::bind_rows(ls_LCx)
  
  return(df_LCx)
  
}


# dose response curve
# 
doseResponse_survFitCstExp <- function(x, time_LCx,
                                       conc_range, npoints){
  
  model_type = x$model_type
  
  # prameters
  mctot <- do.call("rbind", x$mcmc)
  kd <- 10^mctot[, "kd_log10"]
  hb <- 10^mctot[, "hb_log10"]
  
  # all theorical
  k <- 1:length(conc_range)
  j <- 1:npoints
  
  if(model_type == "SD"){
    
    if(is.null(time_LCx)){
      time_LCx = max(x$jags.data$time)
    }
    
    z <- 10^mctot[, "z_log10"]
    kk <- 10^mctot[, "kk_log10"]
    
    dtheo <- lapply(k, function(kit) { # conc
      Surv_SD_LCx(Cw = conc_range[kit],
                  time = time_LCx,
                  kk = kk,
                  kd = kd,
                  z = z,
                  hb = hb)
    })
  }
  if(model_type == "IT"){
    alpha <- 10^mctot[, "alpha_log10"]
    beta <- 10^mctot[, "beta_log10"]
    
    if(is.null(time_LCx)){
      time_LCx = max(x$jags.data$time)
    }
    
    dtheo <- lapply(k, function(kit) { # concentration pour chaque concentration
      Surv_IT_LCx(Cw = conc_range[kit],
                  time = time_LCx,
                  kd = kd,
                  hb = hb,
                  alpha = alpha,
                  beta = beta)
    })
  }
  
  # transpose dtheo
  dtheo <- do.call("rbind", lapply(dtheo, t))
  
  # quantile
  qinf95 <- apply(dtheo, 1, quantile, probs = 0.025, na.rm = TRUE)
  qsup95 <- apply(dtheo, 1, quantile, probs = 0.975, na.rm = TRUE)
  q50 <- apply(dtheo, 1, quantile, probs = 0.5, na.rm = TRUE)
  
  df_dose_Resp = data.frame(concentration = conc_range,
                            q50 = q50,
                            qinf95 = qinf95,
                            qsup95 = qsup95)
}

Surv_SD_LCx <- function(Cw, time, kk, kd, z, hb){
  S <- exp(-hb*time)
  x <- ifelse(Cw > z, 1 - z/Cw, NA)
  tz <- -(1/kd)*log(x)
  y <- ifelse(time > tz,
              exp( kk/kd*Cw*(exp(-kd*tz) -exp(-kd*time))
                   - kk*(Cw-z)*(time - tz)),
              NA)
  return(ifelse(!is.na(x) & !is.na(y), S * y, S))
}


Surv_IT_LCx <- function(Cw, time, kd, hb, alpha, beta){
  D <- Cw*(1-exp(-kd * time))
  S <- exp(-hb * time) * ( 1- 1/(1 + (D/alpha)^(- beta))) 
  return(S)
}

# points for LCx
# 

pointsLCx <- function(df_dose, X_prop){
  
  if(min(df_dose$q50) < X_prop & X_prop < max(df_dose$q50)){
    df.q50 = select(df_dose, c(concentration, q50)) %>%
      add_row(q50 = X_prop) %>%
      arrange(q50) %>%
      mutate(concentration = zoo::na.approx(concentration,q50, na.rm=FALSE))%>%
      filter(q50 == X_prop)
    
    LCX_q50 = df.q50$concentration
    
  } else {
    LCX_q50 = NA
    
    warning(paste("No median for LC", X_prop * 100,
                  " in the range of concentrations under consideration: [",
                  min(df_dose$concentration), ";", max(df_dose$concentration), "]"))
  }
  
  if(min(df_dose$qinf95) < X_prop & X_prop < max(df_dose$qinf95)){
    df.qinf95=select(df_dose, c(concentration,qinf95))%>%
      add_row(qinf95=X_prop)%>%
      arrange(qinf95)%>%
      mutate(concentration = na.approx(concentration,qinf95, na.rm=FALSE))%>%
      filter(qinf95==X_prop)
    
    LCX_qinf95 = df.qinf95$concentration
    
  } else{
    LCX_qinf95 = NA
    
    warning(paste("No 95%inf for LC", X_prop * 100,
                  " in the range of concentrations under consideration: [",
                  min(df_dose$concentration), ";", max(df_dose$concentration), "]"))
  }
  
  
  if(min(df_dose$qsup95) < X_prop & X_prop < max(df_dose$qsup95)){
    df.qsup95=select(df_dose, c(concentration,qsup95))%>%
      add_row(qsup95=X_prop)%>%
      arrange(qsup95)%>%
      mutate(concentration = na.approx(concentration,qsup95, na.rm=FALSE))%>%
      filter(qsup95==X_prop)
    
    LCX_qsup95 = df.qsup95$concentration
    
  } else{
    
    LCX_qsup95 = NA
    warning(paste("No 95%sup for LC", X_prop * 100,
                  " in the range of concentrations under consideration: [",
                  min(df_dose$concentration), ";", max(df_dose$concentration), "]"))
  }
  
  df_LCx <- data.frame(median = LCX_q50,
                       qinf95 = LCX_qinf95,
                       qsup95 = LCX_qsup95,
                       X_prop = X_prop)
  
  return(df_LCx)
  
}

