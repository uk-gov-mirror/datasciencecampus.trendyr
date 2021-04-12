##' Function is used for trend extraction of equally spaced time series data
##'
##' @param inpupt_df dataframe with date in the first column and the variable of interest in the second column
##' @param nowcast Whether the data has weekend effects (default is FALSE)
##'
##' this version has an option for nowcasts: if nowcast = k then the last k observations will be dropped from the data,
##' k forecasts will be produced using net growth and then added to the data and the trend will be finally extracted from this data
##'
##' @author Sonia Mazzi
##' @export
trend_extractor_nowcast <- function(input_df, nowcast = 0){
  #prepare the data
  variable_of_interest <- names(input_df)[2]
  names(input_df) <- c("date", "var_int")
  original_length <- nrow(input_df)
  shorter_length <- original_length - nowcast
  YY <<- input_df$var_int[1:shorter_length]
  #
  #Hyperparameter estimation
  #
  source("likelihoodLLMvard.R")
  min_out <- nlminb(c(0, 0, 1), lik.llm.vard, lower = c(-0.2, -0.5, 0.75), upper = c(0.2, 0.5, 1))$par
  #
  #Diffuse Kalman Filter
  #
  source("dkfLLMvard.R")
  dkf_out <- dkf.llm.vard(min_out, YY)
  #
  #The smoothing stage
  #
  source("smfilter.R")
  sm_out <- smfilt(dkf_out)
  #
  #alpha is the predicted state vector
  alpha <- sm_out$alpha
  ll <- dkf_out$ll
  # trend
  mu <- alpha[1,]
  trend <- ifelse(mu < 0, 0, mu)
  # first derivative
  first_derivative <- alpha[2,]
  first_derivative <- ifelse(mu < 0, 0, first_derivative)
  # net growth
  net_growth <- ifelse(trend >= 1, first_derivative/trend, 0)
  #
  if (nowcast > 0){
    yhat_ng <- rep(0, nowcast)
    NG <- net_growth[shorter_length]
    TREND <- trend[shorter_length]
    for (i in 1:nowcast){
      yhat_ng[i] <- TREND * (1 + NG)^i
    }
    YY <<- c(YY, yhat_ng)
    #
    #Hyperparameter estimation
    #
    source("likelihoodLLMvard.R")
    min_out <- nlminb(c(0, 0, 1), lik.llm.vard, lower = c(-0.2, -0.5, 0.75), upper = c(0.2, 0.5, 1))$par
    #
    #Diffuse Kalman Filter
    #
    source("dkfLLMvard.R")
    dkf_out <- dkf.llm.vard(min_out, YY)
    #
    #The smoothing stage
    #
    source("smfilter.R")
    sm_out <- smfilt(dkf_out)
    #
    #alpha is the predicted state vector
    alpha <- sm_out$alpha
    ll <- dkf_out$ll
    # trend
    mu <- alpha[1,]
    trend <- ifelse(mu < 0, 0, mu)
    # first derivative
    first_derivative <- alpha[2,]
    first_derivative <- ifelse(mu < 0, 0, first_derivative)
    # net growth
    net_growth <- ifelse(trend >= 2, first_derivative/trend, 0)
  }
  #put data together in a data frame
  df_out <- data.frame(Date = as.Date(input_df$date),
                       observed = input_df$var_int, trend, first_derivative, net_growth)
  #write out the results to a csv file
  fn <- paste0("df_out_", variable_of_interest, ".csv")
  write_csv(df_out, fn)
}




