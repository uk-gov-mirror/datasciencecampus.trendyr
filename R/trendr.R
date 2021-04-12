##' Extracts the trend and first derivative using a local linear model in state-space form and the Diffuse Kalman. Gives argument for weekend effect.
##' Filter.
##'
##' @param df The data frame with a time series
##' @param weekend_effects Whether the data has weekend effects (default is FALSE)
##' @param value.colname The name of the column with the values of the series
##' @param time.colname The name of the column with the time values
##' @param output.save Whether the results should save locally or now (default is FALSE)
##' @param output.dir The directory for the output file, by default your working directory
##' @param output.file The file name for the output file, by default out.csv
##' @param output.plot Specify whether you want the output to be plotted using R's basic plotting
##'
##' @author Sonia Mazzi and Michael Hodge
##' @examples
##'   #trendr(
##'     #df = randomData,
##'     #weekend_effects = TRUE,
##'     #value.colname = 'count',
##'     #time.colname = 'time',
##'     #output.dir = 'C:\Users\User\Documents\',
##'     #output.file = 'out.csv',
##'     #output.plot = T)
##' @export
trendr <- function(df = randomData,
                   weekend_effects = FALSE,
                   value.colname = "value",
                   time.colname = "date",
                   output.save = FALSE,
                   output.dir = "./",
                   output.file = "out.csv",
                   output.plot = FALSE) {
  `%>%` <- magrittr::`%>%`

  # Suppress readr read_csv message
  options(readr.num_columns = 0)

  y <- get('df')[value.colname][[1]]
  time <- get('df')[time.colname][[1]]

  DOW <- weekdays(as.Date(time))

  #Hyperparameter estimation
  if (weekend_effects == TRUE) {
    x0 <- c(0, 0, 1, 0)

    #Likelihood
    min_out_values <- nlminb(start = x0, lik.llm.we.vard,
                             lower = c(-0.2, -0.5, 0.85, -0.5),
                             upper = c(0.2, 0.5, 1, 0.5), y = y, DOW = DOW)$par
    # Diffuse Kalman Filter
    dkf_out_values <- dkf.llm.we.vard(min_out_values, y, DOW)


  } else {
    x0 <- c(0, 0, 1)

    #Likelihood
    min_out_values <- nlminb(start = x0, lik.llm.vard, y = y,
                             lower = c(-0.2, -0.5, 0.85),
                             upper = c(0.2, 0.5, 1))$par

    # Diffuse Kalman Filter
    dkf_out_values <- dfkLLMvard(min_out_values, y)
  }

  # The smoothing stage
  sm_out_values <- smfilt(dkf_out_values)

  # Alpha is the predicted state vector
  alpha <- sm_out_values$alpha

  ll <- dkf_out_values$ll

  # Trend
  mu <- alpha[1,]

  # First derivative
  beta <- round(alpha[2,], 5)

  # Put data together in a data frame
  # write to output directory

  df_out <- dplyr::tibble(date = time, observed = y, trend = mu, first_derivative = beta)

  # Write to output directory
  if (output.save == TRUE) {
    write.csv(df_out, paste0(output.dir,output.file), row.names=FALSE)
  }

  # Plot output
  if (output.plot ==  T) {

    timeseries <- as.Date(df_out$date)
    values <- df_out$observed
    trend <- df_out$trend
    first_der <- df_out$first_derivative

    par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
    plot(timeseries, values, type = "l", xlab = time.colname,
         ylab = value.colname, ylim = c(min(first_der), max(values)))

    lines(timeseries, trend, col='red')

    par(new = T)
    plot(timeseries, first_der, col = 'blue', type = "l", xaxt = "n", yaxt = "n",
         ylab = "", xlab = "")
    axis(side = 4)
    mtext("first derivative", side = 4, line = 3)
    legend("bottomright", c("raw", "trend", "first der"),
           col = c("black", "red", "blue"), lty = c(1,1,1), inset = c(0,1), xpd = TRUE, horiz = TRUE, bty = "n")

  }

  return(df_out)
}


