## ---------------------------
##
## Script name: Dose_Response_Functions.R
##
## Purpose of script: Functions used to fit dose-response curves to 
## host-to-vector transmission data. Different functional forms,
## with binomial and betabinomial likelihoods, selected based on AICc
##
## Author: Helene Cecilia
##
## Date Created: 2022-03-23

## -------------------------------
# The models : different functional form to describe dose-response -----
LogFunction = function(log_V, log_beta0, beta1){
  prob_inf = 1 / (1 + exp(-beta1*(log_V - log_beta0)) ) 
  
  return(prob_inf)
}

Ferguson = function(log_V, theta0, theta1){
  exponent = - (log_V / theta0)^theta1
  prob_inf = 1 - exp(exponent)

  return(prob_inf)
}

Hill = function(log_V, gamma0, gamma1){
  prob_inf = (log_V^gamma1) / (gamma0 + log_V^gamma1)
  
  return(prob_inf)
}

# The objective functions (-log(L)) ---------
LogNLL = function(params,k,N,log_V){

  log_beta0 <- params[1]
  beta1 <- params[2]
  
  prob_inf = LogFunction(log_V, log_beta0, beta1)
  
  NLL = 0
  log_like=numeric(length(N))

  for(i in 1:length(N)){
    log_like[i] = -dbinom(k[i], prob = prob_inf[i], size = N[i], log = TRUE)
    if(is.finite(log_like[i])){
      NLL = NLL + log_like[i]
    }else{
      NLL = NLL + 100000
    }
  }
  return(NLL)
}
parnames(LogNLL) <- c("log_beta0", "beta1")

LogBetaBinomNLL = function(params,k,N,log_V){

  log_beta0 <- params[1]
  beta1 <- params[2]
  overdispersion <- params[3]
  
  prob_inf = LogFunction(log_V, log_beta0, beta1)
  
  NLL = 0
  log_like=numeric(length(N))
  
  for(i in 1:length(N)){
    log_like[i] = -dbetabinom(k[i], prob = prob_inf[i], size = N[i], theta = overdispersion, log = TRUE)
    if(is.finite(log_like[i])){
      NLL = NLL + log_like[i]
    }else{
      NLL = NLL + 100000
    }
  }
  return(NLL)
  
}
parnames(LogBetaBinomNLL) <- c("log_beta0", "beta1", "overdispersion")

FergusonNLL = function(params,k,N,log_V){
  theta0 <- exp(params[1]) # theta0 must be positive
  theta1 <- params[2] 
  
  prob_inf = Ferguson(log_V, theta0, theta1)
  
  NLL = 0
  log_like=numeric(length(N))
  
  for(i in 1:length(N)){
    log_like[i] = -dbinom(k[i], prob = prob_inf[i], size = N[i], log = TRUE)
    if(is.finite(log_like[i])){
      NLL = NLL + log_like[i]
    }else{
      NLL = NLL + 100000
    }
  }
  return(NLL)
}

parnames(FergusonNLL) <- c("log.theta0", "theta1")

FergusonBetaBinomNLL = function(params,k,N,log_V){ 
  
  theta0 <- exp(params[1]) # theta0 must be positive
  theta1 <- params[2]
  overdispersion <- params[3]

  prob_inf = Ferguson(log_V, theta0, theta1)
  
  NLL = 0
  log_like=numeric(length(N))
  
  for(i in 1:length(N)){
    log_like[i] = -dbetabinom(k[i], prob = prob_inf[i], size = N[i], theta = overdispersion, log = TRUE)
    if(is.finite(log_like[i])){
      NLL = NLL + log_like[i]
    }else{
      NLL = NLL + 100000
    }
  }
  return(NLL)
}

parnames(FergusonBetaBinomNLL) <- c("log.theta0", "theta1","overdispersion")

HillNLL = function(params,k,N,log_V){ 

  gamma0 <- exp(params[1]) # gamma0 must be positive for proba to be inferior to 1
  gamma1 <- params[2]
  
  prob_inf = Hill(log_V, gamma0, gamma1)
  
  NLL = 0
  log_like=numeric(length(N))
  
  for(i in 1:length(N)){
    log_like[i] = -dbinom(k[i], prob = prob_inf[i], size = N[i], log = TRUE)
    if(is.finite(log_like[i])){
      NLL = NLL + log_like[i]
    }else{
      NLL = NLL + 100000
    }
  }
  return(NLL)
}

parnames(HillNLL) <- c("log.gamma0", "gamma1")

HillBetaBinomNLL = function(params,k,N,log_V){ 
  gamma0 <- exp(params[1]) # gamma0 must be positive
  gamma1 <- params[2] 
  overdispersion <- params[3]
  
  prob_inf = Hill(log_V, gamma0, gamma1)
  
  NLL = 0
  log_like=numeric(length(N))
  
  for(i in 1:length(N)){
    log_like[i] = -dbetabinom(k[i], prob = prob_inf[i], size = N[i], theta = overdispersion, log = TRUE)
    if(is.finite(log_like[i])){
      NLL = NLL + log_like[i]
    }else{
      NLL = NLL + 100000
    }
  }
  return(NLL)
  
}

parnames(HillBetaBinomNLL) <- c("log.gamma0", "gamma1", "overdispersion")


get.binom.model.characteristics <- function(model, type = 'binom'){
  p1.m = model@coef[1]
  p1.up = model@coef[1] + 1.96 * summary(model)@coef[1,2]
  p1.low = model@coef[1] - 1.96 * summary(model)@coef[1,2]
  p1.pv = summary(model)@coef[1,4]
  p2.m = model@coef[2]
  p2.up = model@coef[2] + 1.96 * summary(model)@coef[2,2]
  p2.low = model@coef[2] - 1.96 * summary(model)@coef[2,2]
  p2.pv = summary(model)@coef[2,4]
  LL = summary(model)@m2logL
  if (type == 'betabinom'){
    od.m  = model@coef[3]
    od.up = model@coef[3] + 1.96 * summary(model)@coef[3,2]
    od.low = model@coef[3] - 1.96 * summary(model)@coef[3,2]
    od.pv = summary(model)@coef[1,4]
    return(data.frame(p1.m, p1.low, p1.up, p1.pv,p2.m, p2.low, p2.up, p2.pv, od.m, od.low, od.up, od.pv, LL ))
  } else {
    return(data.frame(p1.m, p1.low, p1.up, p1.pv, p2.m, p2.low, p2.up, p2.pv, LL ))    
  }
}

# likelihood ratio test
lrt <- function (LL0, LL1, df0, df1) {
  L01 <- as.vector( 2 * (LL0 - LL1)) # assuming Negative log likelihoods
  df <- df1 - df0
  list(L01 = L01, df = df,
       "p-value" = pchisq(L01, df, lower.tail = FALSE))
}

optim_Log <- function(dataset, betabinom=FALSE, verbose = TRUE){
  N_df <- dataset$N
  k_df <- dataset$k
  log_V_df <- dataset$log_V
  if(!betabinom){
    start_ <- c(log_beta0 = 5, beta1 = 0.8)
    logl <- LogNLL
  }else{
    start_ <- c(log_beta0 = 5, beta1 = 0.8,
                overdispersion = 0.1)
    logl <- LogBetaBinomNLL
  }
  model_Log <- mle2(minuslogl = logl, 
                    start = start_,
                    data = list(N = N_df, k = k_df, log_V = log_V_df))
  print("first fit")
  if(verbose){
    print(summary(warnings()))
  }

  if(!betabinom){
    start_ <- c(log_beta0 = as.numeric(model_Log@coef[1]),
                beta1 = as.numeric(model_Log@coef[2]))
  }else{
    start_ <- c(log_beta0 = as.numeric(model_Log@coef[1]),
                beta1 = as.numeric(model_Log@coef[2]),
                overdispersion = as.numeric(model_Log@coef[3]))
  }
  
  dif = 1
  while(dif>tolerance) {
    old = model_Log@min
    model_Log <- mle2(minuslogl = logl, 
                      start = start_,
                      data = list(N = N_df, k = k_df, log_V = log_V_df))
    new = model_Log@min 
    dif = abs(old - new)
  }
  print("Optim")
  if(verbose){
    print(summary(warnings()))
  }
  return(model_Log)
}

optim_Ferg <- function(dataset, betabinom=FALSE, verbose = TRUE){
  N_df <- dataset$N
  k_df <- dataset$k
  log_V_df <- dataset$log_V
  if(!betabinom){
    start_ <- c(log.theta0 = log(5), theta1 = 3)
    logl <- FergusonNLL
  }else{
    start_ <- c(log.theta0 = log(5), theta1 = 3,
                overdispersion = 0.46)
    logl <- FergusonBetaBinomNLL
  }
  model_Ferg <- mle2(minuslogl = logl, 
                     start = start_,
                     data = list(N = N_df, k = k_df, log_V = log_V_df))
  print("first fit")
  if(verbose){
    print(summary(warnings()))
  }

  if(!betabinom){
    start_ <- c(log.theta0 = as.numeric(model_Ferg@coef[1]),
                theta1 = as.numeric(model_Ferg@coef[2]))
  }else{
    start_ <- c(log.theta0 = as.numeric(model_Ferg@coef[1]),
                theta1 = as.numeric(model_Ferg@coef[2]),
                overdispersion = as.numeric(model_Ferg@coef[3]))
  }
  
  dif = 1
  while(dif>tolerance) {
    old = model_Ferg@min
    model_Ferg <- mle2(minuslogl = logl, 
                       start = start_,
                       data = list(N = N_df, k = k_df, log_V = log_V_df))
    new = model_Ferg@min 
    dif = abs(old - new)
  }
  print("Optim")
  if(verbose){
    print(summary(warnings()))
  }
  return(model_Ferg)
}

optim_Hill <- function(dataset, betabinom=FALSE, verbose = TRUE){
  N_df <- dataset$N
  k_df <- dataset$k
  log_V_df <- dataset$log_V
  if(!betabinom){
    start_ <- c(log.gamma0 = log(137), gamma1 = 3.2)
    logl <- HillNLL
  }else{
    start_ <- c(log.gamma0 = log(137), gamma1 = 3.2,
                overdispersion = 4)
    logl <- HillBetaBinomNLL
  }
  model_Hill <- mle2(minuslogl = logl, 
                     start = start_,
                     data = list(N = N_df, k = k_df, log_V = log_V_df))
  print("first fit")
  if(verbose){
    print(summary(warnings()))
  }

  if(!betabinom){
    start_ <- c(log.gamma0 = as.numeric(model_Hill@coef[1]),
                gamma1 = as.numeric(model_Hill@coef[2]))
  }else{
    start_ <- c(log.gamma0 = as.numeric(model_Hill@coef[1]),
                gamma1 = as.numeric(model_Hill@coef[2]),
                overdispersion = as.numeric(model_Hill@coef[3]))
  }
  
  dif = 1
  while(dif>tolerance) {
    old = model_Hill@min
    model_Hill <- mle2(minuslogl = logl, 
                       start = start_,
                       data = list(N = N_df, k = k_df, log_V = log_V_df))
    new = model_Hill@min 
    dif = abs(old - new)
  }
  print("Optim")
  if(verbose){
    print(summary(warnings()))
  }
  return(model_Hill)
}

# Convert from RNA-emia (in the initial papers for DENV) to PFU, using Blaney et al. 2005 (hard-coded)
# returns the data (for plots)
conversion_and_small_values_correction <- function(data, serotype = NULL, solution){
  conversion_factor <- c(1.9,2.8,2.5,1.9) # from Blaney et al. 2005
  data_convert <- data
  if(!is.null(serotype)){
    data_convert <- data_convert[data_convert$serotype == serotype,]
    conv <- conversion_factor[serotype]
    data_convert$PFU_V <- 10**(data_convert$log_V - conv)
  }
  # Caution : currently if you use this function on a dataset that does not specify a serotype,
  # no unit conversion (from RNA-emia to PFU) is performed and the dataset must already have a PFU_V column
  if(solution == 1){
    data_convert <- filter(data_convert, !(PFU_V < 20 & k == 0))
    data_convert$PFU_V[data_convert$PFU_V < 20 & data_convert$k > 0] <- 10
    data_convert$log_V <- log10(data_convert$PFU_V + 1)
  }else if(solution == 2){
    # as the log_V was initially entered as 0 for absence, little trick to re-identify those after conversion
    data_convert$PFU_V[log10(data_convert$PFU_V) == -conv] <- 0
    data_convert$log_V <- log10(data_convert$PFU_V + 1)
  }
  return(data_convert)
}

# Convert from RNA-emia (in the initial papers for DENV) to PFU, using Blaney et al. 2005 (hard-coded)
# then performs functional form selection
dose_resp_conversion_and_selection <- function(data, serotype = NULL, solution, filename, verbose = TRUE){
  conversion_factor <- c(1.9,2.8,2.5,1.9) # from Blaney et al. 2005
  data_convert <- data
  if(!is.null(serotype)){
    data_convert <- data_convert[data_convert$serotype == serotype,]
    conv <- conversion_factor[serotype]
    data_convert$PFU_V <- 10**(data_convert$log_V - conv)
  }
  # Caution : currently if you use this function on a dataset that does not specify a serotype,
  # no unit conversion (from RNA-emia to PFU) is performed and the dataset must already have a PFU_V column
  if(solution == 1){
    data_convert <- filter(data_convert, !(PFU_V < 20 & k == 0))
    data_convert$PFU_V[data_convert$PFU_V < 20 & data_convert$k > 0] <- 10
    data_convert$log_V <- log10(data_convert$PFU_V + 1)
  }else if(solution == 2){
    # as the log_V was initially entered as 0 for absence, little trick to re-identify those after conversion
    data_convert$PFU_V[log10(data_convert$PFU_V) == -conv] <- 0
    data_convert$log_V <- log10(data_convert$PFU_V + 1)
  }
  dose_resp_function_selection(data_convert, filename, verbose)
  return(data_convert)
}

# functional form selection only (no assumption on units)
dose_resp_function_selection <- function(data, filename, verbose = TRUE){

  print("binom_Log")
  binom_Log <- optim_Log(data, verbose = verbose)
  print("binom_Log")
  betabinom_Log <- optim_Log(data, betabinom = T, verbose = verbose)

  print("binom_Ferg")
  binom_Ferg <- optim_Ferg(data, verbose = verbose)
  print("betabinom_Ferg")
  betabinom_Ferg <- optim_Ferg(data, betabinom = T, verbose = verbose)
  
  print("binom_Hill")
  binom_Hill <- optim_Hill(data, verbose = verbose)
  print("betabinom_Hill")
  betabinom_Hill <- optim_Hill(data, betabinom = T, verbose = verbose)

  saveRDS(betabinom_Hill, file = paste0(save_path,"betabinom_Hill_",filename,".rds"))
  saveRDS(betabinom_Ferg, file = paste0(save_path,"betabinom_Ferg_",filename,".rds"))
  saveRDS(betabinom_Log, file = paste0(save_path,"betabinom_Log_",filename,".rds"))
  saveRDS(binom_Hill, file = paste0(save_path,"binom_Hill_",filename,".rds"))
  saveRDS(binom_Ferg, file = paste0(save_path,"binom_Ferg_",filename,".rds"))
  saveRDS(binom_Log, file = paste0(save_path,"binom_Log_",filename,".rds"))

  AICs = AICc(binom_Log,
             betabinom_Log,
             binom_Ferg,
             betabinom_Ferg,
             binom_Hill,
             betabinom_Hill)
  row.names(AICs) = c("binom_Log",
                       "betabinom_Log",
                       "binom_Ferg",
                       "betabinom_Ferg",
                       "binom_Hill",
                       "betabinom_Hill")
  print(AICs)
}

# computes uncertainty around the curve by sampling many trajectories and extracting quantiles
# and estimates dose 50
compute_uncertainty_dose_resp <- function(best_model, samples, LogV_vector, fct, fct_name = "",
                                          type = "betabinom"){
  Prob <- matrix(NA,ncol = length(LogV_vector), nrow = samples)
  dose50 <- NULL
  if(fct_name == "Log"){
    MV.Coefs = rmvn(n=samples, mu = best_model@coef, V = vcov(best_model))
  }else{
    if(type == "betabinom"){
      MV.Coefs = rtmvnorm(n = samples, mean = best_model@coef, sigma = vcov(best_model),
                          lower = c(0,-Inf,-Inf))
    }else{
      MV.Coefs = rtmvnorm(n = samples, mean = best_model@coef, sigma = vcov(best_model),
                          lower = c(0,-Inf))
    }
  }
  correct <- 0
  for (ii in 1:dim(MV.Coefs)[1]){
    if(fct_name != "Log"){
      replicate <- fct(LogV_vector, exp(MV.Coefs[ii,1]), MV.Coefs[ii,2])
    }else{
      replicate <- fct(LogV_vector, MV.Coefs[ii,1], MV.Coefs[ii,2])
    }
    # browser()
    while(sum(replicate>1) > 0 || sum(replicate < 0) > 0 || sum(is.na(replicate)) > 0){ 
      # hard fix : resample as long as we get trajectories with proba > 1 or < 0 and with NAs
      # condition should not be met now that we use rtmvnorm (recorded in correct)
      print(MV.Coefs[ii,])
      new_coef = rmvn(n=1, mu = best_model@coef, V = vcov(best_model))
      if(fct_name != "Log"){
        replicate <- fct(LogV_vector, exp(new_coef[1]), new_coef[2])
      }else{
        replicate <- fct(LogV_vector, new_coef[1], new_coef[2])
      }
      correct <- correct + 1
    }
    Prob[ii,] = replicate
    if(fct_name == "Log"){
      d <- MV.Coefs[ii,1] 
      dose50 <- c(dose50,d)
    }else{
      d <- LogV_vector[which.min(abs(replicate-0.5))]
      dose50 <- c(dose50,d)
    }
  }
  if(correct > 0){
    print(paste("nb of artificial resampling",correct, sep = ": "))
  }
  traj <- matrix_to_stack(Prob,
                          col_to_col = "log_V", row_to_col = "id",
                          value_col = "prob_inf")
  traj$log_V <- rep(LogV_vector, each = samples)
  
  if(fct_name == "Log"){
    fit <- fct(LogV_vector,
               best_model@coef[1],
               best_model@coef[2])
    dose50_estim <- best_model@coef[1] 
  }else{
    fit <- fct(LogV_vector,
               exp(best_model@coef[1]),
               best_model@coef[2])
    dose50_estim <- LogV_vector[which.min(abs(fit-0.5))]
  }
  envelope <- NULL
  for(i in seq(length(LogV_vector))){
    time_point <- Prob[,i]
    
    max95_time <- quantile(time_point, probs = 0.975, na.rm = T) 
    min95_time <- quantile(time_point, probs = 0.025, na.rm = T)
    
    envelope <- rbind(envelope, data.frame(dose = LogV_vector[i],
                                           min_95 = min95_time,
                                           max_95 = max95_time,
                                           fit = fit[i],
                                           type = fct_name))
    
  }
  return(list("dose50_CI" = dose50, "traj" = traj, "envelope" = envelope,
              "fit" = fit, "dose50_estim" = dose50_estim))
}


