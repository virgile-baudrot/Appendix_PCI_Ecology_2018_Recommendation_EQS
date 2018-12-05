## Varying period between pulses


depurationTime_X_val <- function(x, X_dep = 0.95){
  df <- do.call("rbind", x$mcmc)
  depurationTime <- -log(1 - X_dep ) / 10^df[, "kd_log10"]
  return(depurationTime)
}


depuration_df <- function(x, pulse = 0, Maximal_conc = NULL, recovery_time = NULL, max_time = 15, model_type = NULL, draws = NULL){
  
  if( pulse == 2){
    profile_df <- data.frame(time = c(0, 0.1, 0.9, 1, recovery_time+1, recovery_time+1.1, recovery_time +1.9, recovery_time+2, max_time),
                             conc = Maximal_conc * c(0, 1,   1, 0, 0, 1,   1, 0, 0))
    
  }else{
    mult_time <- rep((recovery_time+1.1), 4)
    add_time <- c(0, 0.1, 0.9, 1)
    
    time_vec <- 0*mult_time + add_time
    i = 0
    max_tps = 1.1
    while(max_tps < max_time){
      i=i+1
      time_vec <-  c(time_vec, i*mult_time + add_time)
      max_tps <- max(time_vec)
    }
    profile_df <- data.frame(time = time_vec,
                             conc = Maximal_conc * rep(c(0, 1,   1, 0),i+1))%>%
      filter(time < max_time)
    
    profile_df <- bind_rows(profile_df, data.frame(time = max_time, conc = profile_df$conc[nrow(profile_df)]))
    
    cat("Indication: while loop done")
  }
  
  
  ## Simulation
  predict_df <- solveTKTD_MFx(x = x,
                              model_type = model_type,
                              time_conc = profile_df$time,
                              conc = profile_df$conc,
                              draws = draws,
                              interpolate_time = NULL,
                              interpolate_length = 1e2,
                              interpolate_method = "linear")
  
  return(list(profile_df = profile_df,
              predict_df = predict_df$df_survival,
              predict_df_internalConc = predict_df$df_internalConc))
  
}

depuration_plot <- function(depuration_df,  conc_max_CONC = NULL, conc_max_INTCONC = NULL, title){
  
  ##
  ## Plot profiles
  ##
  plt_profile = ggplot() + theme_bw()+
    labs(x = "Time", y = "Concentration", title = title) +
    lims( y = c(0,conc_max_CONC)) +
    geom_line(data = depuration_df$profile_df,
              aes(x = time, y = conc))
  
  ##
  ## Plot internal concentration
  ##
  
  df_quantile_IntConc <- depuration_df$predict_df_internalConc
  mat_IntConc <- select(df_quantile_IntConc, - c(time, conc))
  
  q50_IntConc <- apply(mat_IntConc, 1, quantile, prob = 0.5)
  qinf95_IntConc <- apply(mat_IntConc, 1, quantile, prob = 0.025)
  qsup95_IntConc <- apply(mat_IntConc, 1, quantile, prob = 0.975)
  
  plt_prediction_IntConc = ggplot() + theme_bw() +
    labs(x = "Time", y = "Internal Concentration") +
    lims( y = c(0,conc_max_INTCONC)) +
    geom_ribbon(aes(x = df_quantile_IntConc$time, ymin = qinf95_IntConc, ymax = qsup95_IntConc), fill = "lightgrey") +
    geom_line(aes(x = df_quantile_IntConc$time, y = q50_IntConc ), color = "orange")
  
  ##
  ## Plot internal concentration
  ##
  
  df_quantile <- depuration_df$predict_df
  mat <- select(df_quantile, - c(time, conc))
  
  q50 <- apply(mat, 1, quantile, prob = 0.5)
  qinf95 <- apply(mat, 1, quantile, prob = 0.025)
  qsup95 <- apply(mat, 1, quantile, prob = 0.975)
  
  plt_prediction = ggplot() + theme_bw() +
    labs(x = "Time", y = "Survival") +
    lims( y = c(0,1)) +
    geom_ribbon(aes(x = df_quantile$time, ymin = qinf95, ymax = qsup95), fill = "lightgrey") +
    geom_line(aes(x = df_quantile$time, y = q50 ), color = "orange")
  
  grid.arrange(plt_profile, plt_prediction_IntConc, plt_prediction, ncol = 1, heights = c(1,1, 1.5))
}

depuration_plot_SDandIT <- function(depuration_df_SD,depuration_df_IT,  conc_max_CONC = NULL, conc_max_INTCONC = NULL, title){
  
  ##
  ## Plot profiles
  ##
  plt_profile = ggplot() + theme_bw()+
    labs(x = "Time", y = "Concentration", title = title) +
    lims( y = c(0,conc_max_CONC)) +
    geom_line(data = depuration_df_SD$profile_df,
              aes(x = time, y = conc)) +
    geom_line(data = depuration_df_IT$profile_df,
              aes(x = time, y = conc))
  
  ##
  ## Plot internal concentration
  ##
  
  df_quantile_IntConc_SD <- depuration_df_SD$predict_df_internalConc
  mat_IntConc_SD <- select(df_quantile_IntConc_SD, - c(time, conc))
  df_quantile_IntConc_SD$q50_IntConc <- apply(mat_IntConc_SD, 1, quantile, prob = 0.5)
  df_quantile_IntConc_SD$qinf95_IntConc <- apply(mat_IntConc_SD, 1, quantile, prob = 0.025)
  df_quantile_IntConc_SD$qsup95_IntConc <- apply(mat_IntConc_SD, 1, quantile, prob = 0.975)
  
  df_quantile_IntConc_IT <- depuration_df_IT$predict_df_internalConc
  mat_IntConc_IT <- select(df_quantile_IntConc_IT, - c(time, conc))
  df_quantile_IntConc_IT$q50_IntConc <- apply(mat_IntConc_IT, 1, quantile, prob = 0.5)
  df_quantile_IntConc_IT$qinf95_IntConc <- apply(mat_IntConc_IT, 1, quantile, prob = 0.025)
  df_quantile_IntConc_IT$qsup95_IntConc <- apply(mat_IntConc_IT, 1, quantile, prob = 0.975)
  
  plt_prediction_IntConc = ggplot() + theme_bw() +
    labs(x = "Time", y = "Internal Concentration") +
    lims( y = c(0,conc_max_INTCONC)) +
    # SD
    geom_ribbon(data = df_quantile_IntConc_SD,
                aes(x = time, ymin = qinf95_IntConc, ymax = qsup95_IntConc), fill = "grey", alpha = 0.2) +
    geom_line(data = df_quantile_IntConc_SD,
              aes(x = time, y = qinf95_IntConc), color = "grey") +
    geom_line(data = df_quantile_IntConc_SD,
              aes(x = time, y = qsup95_IntConc), color = "grey") +
    geom_line(data = df_quantile_IntConc_SD,
              aes(x = time, y = q50_IntConc )) +
    # IT
    geom_ribbon(data = df_quantile_IntConc_IT,
                aes(x = time, ymin = qinf95_IntConc, ymax = qsup95_IntConc), fill = "grey", alpha = 0.2) +
    geom_line(data = df_quantile_IntConc_IT,
              aes(x = time, y = qinf95_IntConc), color = "grey", linetype = 2) +
    geom_line(data = df_quantile_IntConc_IT,
              aes(x = time, y = qsup95_IntConc), color = "grey", linetype = 2) +
    geom_line(data = df_quantile_IntConc_IT,
              aes(x = time, y = q50_IntConc ), linetype = 2)
  
  ##
  ## Plot internal concentration
  ##
  
  df_quantile_SD <- depuration_df_SD$predict_df
  mat_SD <- select(df_quantile_SD, - c(time, conc))
  df_quantile_SD$q50 <- apply(mat_SD, 1, quantile, prob = 0.5)
  df_quantile_SD$qinf95 <- apply(mat_SD, 1, quantile, prob = 0.025)
  df_quantile_SD$qsup95 <- apply(mat_SD, 1, quantile, prob = 0.975)
  
  df_quantile_IT <- depuration_df_IT$predict_df
  mat_IT <- select(df_quantile_IT, - c(time, conc))
  df_quantile_IT$q50 <- apply(mat_IT, 1, quantile, prob = 0.5)
  df_quantile_IT$qinf95 <- apply(mat_IT, 1, quantile, prob = 0.025)
  df_quantile_IT$qsup95 <- apply(mat_IT, 1, quantile, prob = 0.975)
  
  plt_prediction = ggplot() + theme_bw() +
    labs(x = "Time", y = "Survival") +
    lims( y = c(0,1)) +
    # SD
    geom_ribbon(data = df_quantile_SD,
                aes(x = time, ymin = qinf95, ymax = qsup95), fill = "grey", alpha = 0.2) +
    geom_line(data = df_quantile_SD,
              aes(x = time, y = qinf95 ), color = "grey") +
    geom_line(data = df_quantile_SD,
              aes(x = time, y = qsup95 ), color = "grey") +
    geom_line(data = df_quantile_SD,
              aes(x = time, y = q50 )) +
    # IT
    geom_ribbon(data = df_quantile_IT,
                aes(x = time, ymin = qinf95, ymax = qsup95), fill = "grey", alpha = 0.2) +
    geom_line(data = df_quantile_IT,
              aes(x = time, y = qinf95 ), color = "grey", linetype = 2) +
    geom_line(data = df_quantile_IT,
              aes(x = time, y = qsup95 ), color = "grey", linetype = 2) +
    geom_line(data = df_quantile_IT,
              aes(x = time, y = q50 ), linetype = 2) 
  
  grid.arrange(plt_profile, plt_prediction_IntConc, plt_prediction, ncol = 1, heights = c(1,1, 1.5))
}
