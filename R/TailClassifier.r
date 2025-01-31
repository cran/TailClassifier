# The function is developed under R/4.4.1 Platform: x86_64-apple-darwin20 (64-bit).
# Author/Maintainer: Jialin Zhang (JZ); jzhang (at) math.msstate.edu

#  The function TailClassifier() suggests one of the following types of tail for your discrete data: 1) Power decaying tail; 2) Sub-exponential decaying tail; and 3) Near-exponential decaying tail. The function also provides an estimate of the parameter for the classified-distribution as a reference.

## Guidance:

#' @title
#' Tail Classifier
#' @description
#' The function TailClassifier() suggests one of the following types of tail for your discrete data: 1) Power decaying tail; 2) Sub-exponential decaying tail; and 3) Near-exponential decaying tail. The function also provides an estimate of the parameter for the classified-distribution as a reference.
#'
#' @param sample_frequencies The frequency counts for your discrete sample data.
#' @param v_left The starting point of tail profile. 20 is recommended. A smaller v_left may lead to unreliable results. A larger v_left might be adopted if the sample size is extremely large.
#' @param v_right The ending point of tail profile. Recommendation is no more than 100 regardless of sample size.
#' @param plot_lower The lower range of v-axis.
#' @param plot_upper The upper range of v-axis.
#' @param Plot0_title The title for Plot0. The default is ``Plot 0 of Heavy Tail Detection''.
#' @param Plot1_title The title for Plot1. The default is ``Plot 1 of Heavy Tail Detection''.
#' @param Plot2_title The title for Plot2. The default is ``Plot 2 of Heavy Tail Detection''.
#' @param Plot3_title The title for Plot3. The default is ``Plot 3 of Heavy Tail Detection''.
#' @param C_Level The confidence level of confidence intervals in results. The default is 0.95.
#' @param ConfidenceBand TRUE if a confidence band is requested. FALSE otherwise.
#' @param Plot_0_y_limit_lower_extend Modify the y limit in Plot 0 to allow the confidence band to correctly display in different scenarios.
#' @param Plot_0_y_limit_upper_extend Modify the y limit in Plot 1 to allow the confidence band to correctly display in different scenarios.
#' @param Plot_1_y_limit_lower_extend Modify the y limit in Plot 2 to allow the confidence band to correctly display in different scenarios.
#' @param Plot_1_y_limit_upper_extend Modify the y limit in Plot 3 to allow the confidence band to correctly display in different scenarios.
#' @param Plot_2_y_limit_lower_extend Modify the y limit in Plot 0 to allow the confidence band to correctly display in different scenarios.
#' @param Plot_2_y_limit_upper_extend Modify the y limit in Plot 1 to allow the confidence band to correctly display in different scenarios.
#' @param Plot_3_y_limit_lower_extend Modify the y limit in Plot 2 to allow the confidence band to correctly display in different scenarios.
#' @param Plot_3_y_limit_upper_extend Modify the y limit in Plot 3 to allow the confidence band to correctly display in different scenarios.
#' @param subtitle_size Controls the subtitle font size.
#' @param axis_label_size Controls the size of axis labels.
#' @param axis_ticks_size Controls the size of axis tick numbers.

#' @return A statement on the type of tail.
#' @examples

#' ## Power Example
#' # Generate data from power decaying distribution with parameter 1.5
#' rpar <- function(n, a, xm = 1) {
#'   v <- runif(n)
#'   xm / v^(1.0/a)
#' }
#' dpar <- function(x, a, xm = 1){
#' return(a*xm^a/(x^(a+1)))
#' }
#' set.seed(2023)
#' data <- floor(rpar(1000, 0.5)) # lambda = 1.5
#' Result <- TailClassifier(table(data), plot_lower = 5, plot_upper = 400, v_left = 20, v_right = 54,
#'  Plot_0_y_limit_upper_extend = 8)
#' ## display the results
#' Result
#' ## call the classification decision
#' Result$Type
#' ## call the confidence intervals for the parameters
#' data.frame(Result$Results[3])[,c(1,3:4)]
#' ## call a specific plot
#' Result$Results[[1]][1]
#' Result$Results[[1]][2]
#' Result$Results[[1]][3]
#' Result$Results[[1]][4]
#' ## check the rank of possible type of tails
#' Result$Rank
#'

#'
#' @import ggplot2 cowplot utils
#' @importFrom scales extended_breaks
#' @importFrom stats qnorm
#' @export

TailClassifier <- function(sample_frequencies, v_left= 20, v_right= min(floor(sum(sample_frequencies)/20), sum(sample_frequencies[sample_frequencies > 1]) - 1),
                           plot_lower = v_left, plot_upper = v_right,
                           Plot0_title = "Plot 0 of Heavy Tail Detection \n \n",
                           Plot1_title = "Plot 1 of Heavy Tail Detection",
                           Plot2_title = "Plot 2 of Heavy Tail Detection",
                           Plot3_title = "Plot 3 of Heavy Tail Detection",
                           C_Level = 0.95, ConfidenceBand = TRUE,
                           Plot_0_y_limit_lower_extend = 1.5, Plot_0_y_limit_upper_extend = 1.5,
                           Plot_1_y_limit_lower_extend = 0.25, Plot_1_y_limit_upper_extend = 0.25,
                           Plot_2_y_limit_lower_extend = 0.25, Plot_2_y_limit_upper_extend = 0.25,
                           Plot_3_y_limit_lower_extend = 0.25, Plot_3_y_limit_upper_extend = 0.25,
                           subtitle_size = 14, axis_label_size = 12, axis_ticks_size = 10) {

  if(ConfidenceBand == T){alpha_geom = 0.25} else{alpha_geom = 0}

  ## sample size
  n=sum(sample_frequencies)
  ## stop possible mis-use
  if (is.integer(sample_frequencies) == FALSE) {warning("Your input must be sample frequencies (counts)!")}
  if (floor(v_left) != v_left) {stop("v_left must be an integer!")}
  if (floor(v_right) != v_right) {stop("v_right must be an integer!")}
  if (v_left < 16) {warning("v_left must be great than exp(exp(1))!")}
  if (v_right <= v_left) {stop("Argument 'v_right' should be much larger than 'v_left'!")}
  if (v_right > n-1) {stop("Argument 'v_right' should be less than n (sample size)!")}
  if (v_right > 5000) {warning("Argument 'v_right' may not surpass 3000. Higher v_right leads to inreliable results due to computation accuracy.")}
  if( v_right > floor(sum(sample_frequencies)/2)){warning("No interval inference when v_right is greater than half of the sample size!")}
  ## zif() for the method
  z1f <- function(sample_freq)
  {
    khat<-length(sample_freq)
    n<-sum(sample_freq)
    prod<-rep(1,khat)
    zf<-rep(0,n)
    for (v in 1:n){
      for (k in 1:khat){
        if (sample_freq[k]>=1){
          zf[v] = zf[v]+sample_freq[k]/n*prod[k];
          prod[k] = prod[k]*(1-(sample_freq[k]-1)/(n-v));
        }
      }
    }
    return(zf)
  }
  ## obtain the Tv and ln (Tv)
  z1v = z1f(sample_frequencies)
  zeros = sum(z1v==0)
  z1v = z1v[z1v!=0]
  tv = z1v*(0:(n-1-zeros))
  if((v_right < zeros) == T){stop("The data suggest an exponential or thinner tail.")}
  vf = ((v_left+1):(v_right-zeros+1))
  ts = tv[vf]
  lts = log(ts)

  ## function for z_u,v variance
  sigma_z_v <- function(z1v, v){ # v >= 3, 2v <= sum(z1v>0)
    v <- v[v>=3 & (2*v<=length(z1v))]
    sigmas <- data.frame(V = v, Sigmas = numeric(length(v)))
    for (i in 1:length(v)) {
      tmp <- z1v[2*v[i]]/(v[i]+1)^2 - 2*v[i]/(v[i]+1)^2*(z1v[2*v[i]-2]-z1v[2*v[i]-1])+v[i]^2/(v[i]+1)^2*(z1v[2*v[i]-4]-2*z1v[2*v[i]-3]+z1v[2*v[i]-2])-v[i]^2/(v[i]+1)^2*(z1v[i]/v[i]-(z1v[v[i]-2]-z1v[v[i]-1]))^2
      if(tmp<0) {sigmas[i,2] <- 0} else{
        sigmas[i,2] <- sqrt(tmp)
      }
    }
    return(sigmas)
  }

  v_max <- floor(length(z1v[-1])/2)
  sigmas <- sigma_z_v(z1v, 1:v_max)
  margin_errors <- qnorm(1-(1-C_Level)/2)*(sigmas$V+1)*sigmas$Sigmas/sqrt(n)
  ConfidenceIntervals <- data.frame(v = sigmas$V, Lower = sigmas$V*(z1v[sigmas$V+1]-margin_errors), Upper = sigmas$V*(z1v[sigmas$V+1]+margin_errors), PointEstimate = sigmas$V*(z1v[sigmas$V+1]))

  plot_theme <- theme(
    plot.title = element_text(hjust = 0.5, size = subtitle_size),
    axis.title = element_text(size = axis_label_size),
    axis.text = element_text(size = axis_ticks_size)
  )

  ## The four plots
  ## plot 0 is Tv ~ v; heavier than exponential tails converge to infinite, thinner than exponential tails converge to zero
  df_0 <- merge(data.frame(v = (0:(length(tv)-1))[plot_lower:plot_upper], tv = tv[plot_lower:plot_upper]), ConfidenceIntervals, by = "v", all = TRUE)
  df_0$Lower <- ifelse(is.na(df_0$Lower), df_0$tv, df_0$Lower)
  df_0$Upper <- ifelse(is.na(df_0$Upper), df_0$tv, df_0$Upper)
  df_0 <- df_0[-1,]

  plot0 <- ggplot(df_0, aes(x=v, y=tv)) + theme_bw() + xlab("v") + ylab(bquote(T[v])) +
    ylim((min(tv[plot_lower:plot_upper])) - Plot_0_y_limit_lower_extend,(max(tv[plot_lower:plot_upper])) + Plot_0_y_limit_upper_extend) +
    scale_x_continuous(limits = c(plot_lower, plot_upper)) + ggtitle(Plot0_title) + plot_theme +
    geom_line(na.rm=TRUE) + geom_ribbon(aes(ymin=Lower,ymax=Upper), alpha=alpha_geom)

  ## plot 1 is ln Tv ~ ln v; power decaying tails become linear; concave lines indicate tails are thinner than power decaying.
  log_v <- log(1:length(tv))[plot_lower:plot_upper]
  log_tv <- log(tv)[(plot_lower+1):(plot_upper+1)]
  log_v <- log_v[!is.infinite(log_tv) & !is.na(log_tv)]
  log_tv <- log_tv[!is.infinite(log_tv) & !is.na(log_tv)]

  ConfidenceIntervals1 <- suppressWarnings(log(ConfidenceIntervals))
  names(ConfidenceIntervals1)[1] <- "log_v"
  df_1 <- merge(data.frame(log_v = log_v, log_tv = log_tv), ConfidenceIntervals1, by = "log_v", all = TRUE)
  df_1$Lower <- ifelse(is.na(df_1$Lower), df_1$log_tv, df_1$Lower)
  df_1$Upper <- ifelse(is.na(df_1$Upper), df_1$log_tv, df_1$Upper)
  df_1$PointEstimate <- ifelse(is.na(df_1$PointEstimate), df_1$log_tv, df_1$PointEstimate)

  plot1 <- ggplot(df_1, aes(x=log_v, y=log_tv)) + theme_bw() + xlab("v") + ylab(bquote(ln~T[v])) +
    ylim(floor(4*min(log_tv[plot_lower:plot_upper], na.rm = TRUE))/4 - Plot_1_y_limit_lower_extend,
         floor(4*max(log_tv[plot_lower:plot_upper], na.rm = TRUE))/4+Plot_1_y_limit_upper_extend) +
    scale_x_continuous(limits = c((log(plot_lower)), (log(plot_upper))), labels = ~floor(exp(.)), breaks = ~ c(scales::extended_breaks()(.x)), sec.axis = sec_axis(~., name = "ln v")) + ggtitle(Plot1_title) + plot_theme +
    geom_line(na.rm=TRUE) + geom_ribbon(aes(ymin=Lower,ymax=Upper), alpha=alpha_geom)

  ## plot 2 is ln Tv ~ ln ln v; sub-exp tails become linear; thinner than sub-exp tails are concave; heavier than sub-exp tails are convex.
  plot_lower <- max(2, plot_lower) # this is to avoid an undefined ln ln (1)
  loglog_v <- log(log(1:length(tv)))[plot_lower:plot_upper]
  log_tv <- log(tv)[(plot_lower+1):(plot_upper+1)]
  loglog_v <- loglog_v[!is.infinite(log_tv) & !is.na(log_tv)]
  log_tv <- log_tv[!is.infinite(log_tv) & !is.na(log_tv)]
  log_tv <- log_tv[!is.infinite(loglog_v) & !is.na(loglog_v)]
  loglog_v <- loglog_v[!is.infinite(loglog_v) & !is.na(loglog_v)]

  ConfidenceIntervals2 <- ConfidenceIntervals1
  names(ConfidenceIntervals2)[1] <- "loglog_v"
  ConfidenceIntervals2$loglog_v <- suppressWarnings(log(ConfidenceIntervals1$log_v))
  df_2 <- merge(data.frame(loglog_v = loglog_v, log_tv = log_tv), ConfidenceIntervals2, by = "loglog_v", all = TRUE)
  df_2$Lower <- ifelse(is.na(df_2$Lower), df_2$log_tv, df_2$Lower)
  df_2$Upper <- ifelse(is.na(df_2$Upper), df_2$log_tv, df_2$Upper)
  df_2$PointEstimate <- ifelse(is.na(df_2$PointEstimate), df_2$log_tv, df_2$PointEstimate)


  plot2 <- ggplot(df_2, aes(x=loglog_v, y=log_tv)) + theme_bw() + xlab("v") + ylab(bquote(ln~T[v])) +
    ylim(floor(4*min(log_tv[plot_lower:plot_upper], na.rm = TRUE))/4 - Plot_2_y_limit_lower_extend,
         floor(4*max(log_tv[plot_lower:plot_upper], na.rm = TRUE))/4+Plot_2_y_limit_upper_extend) +
    scale_x_continuous(limits = c((log(log((plot_lower)))), (log(log(plot_upper)))), labels = ~floor(exp(exp(.))), breaks = ~ c(scales::extended_breaks()(.x)), sec.axis = sec_axis(~., name = "ln ln v")) + ggtitle(Plot2_title) + plot_theme +
    geom_line(na.rm=TRUE) + geom_ribbon(aes(ymin=Lower,ymax=Upper), alpha=alpha_geom)

  ## plot 3 is ln Tv ~ ln ln ln v; near-exp tails become linear, thinner than near-exp tails are concave; heavier than near-exp tails are convex.
  plot_lower <- max(3, plot_lower)
  suppressWarnings({logloglog_v <- log(log(log(1:length(tv))))[plot_lower:plot_upper]})
  log_tv <- log(tv)[(plot_lower+1):(plot_upper+1)]
  logloglog_v <- logloglog_v[!is.infinite(log_tv) & !is.na(log_tv)]
  log_tv <- log_tv[!is.infinite(log_tv) & !is.na(log_tv)]
  log_tv <- log_tv[!is.infinite(logloglog_v) & !is.na(logloglog_v)]
  logloglog_v <- logloglog_v[!is.infinite(logloglog_v) & !is.na(logloglog_v)]

  ConfidenceIntervals3 <- ConfidenceIntervals2
  names(ConfidenceIntervals3)[1] <- "logloglog_v"
  ConfidenceIntervals3$logloglog_v <- suppressWarnings(log(ConfidenceIntervals2$loglog_v))
  df_3 <- merge(data.frame(logloglog_v = logloglog_v, log_tv = log_tv), ConfidenceIntervals3, by = "logloglog_v", all = TRUE)
  df_3$Lower <- ifelse(is.na(df_3$Lower), df_3$log_tv, df_3$Lower)
  df_3$Upper <- ifelse(is.na(df_3$Upper), df_3$log_tv, df_3$Upper)
  df_3$PointEstimate <- ifelse(is.na(df_3$PointEstimate), df_3$log_tv, df_3$PointEstimate)

  plot3 <- ggplot(df_3, aes(x=logloglog_v, y=log_tv)) + theme_bw() + xlab("v") + ylab(bquote(ln~T[v])) +
    ylim(floor(4*min(log_tv[plot_lower:plot_upper], na.rm = TRUE))/4 - Plot_3_y_limit_lower_extend,
         floor(4*max(log_tv[plot_lower:plot_upper], na.rm = TRUE))/4+Plot_3_y_limit_upper_extend) +
    scale_x_continuous(limits = c((log(log(log((plot_lower))))), log(log(log(plot_upper)))), labels = ~floor(exp(exp(exp(.)))), breaks = ~ c(scales::extended_breaks()(.x)), sec.axis = sec_axis(~., name = "ln ln ln v")) + ggtitle(Plot3_title) + plot_theme +
    geom_line(na.rm=TRUE) + geom_ribbon(aes(ymin=Lower,ymax=Upper), alpha=alpha_geom)

  final_plot <- structure(list(plot0, plot1, plot2, plot3), class = "multiplot")

  lwr_1 <- max(0, 1/((ConfidenceIntervals1$Upper[which(ConfidenceIntervals1$log_v==log(v_right))] - ConfidenceIntervals1$Upper[which(ConfidenceIntervals1$log_v==log(v_left))])/(log(v_right) - log(v_left))))

  if(is.na(ConfidenceIntervals1$Lower[which(ConfidenceIntervals1$log_v==log(v_right))])){warning("v_right is too large!"); upr_1 <- Inf} else{
    upr_1_tmp <- 1/((ConfidenceIntervals1$Lower[which(ConfidenceIntervals1$log_v==log(v_right))] - ConfidenceIntervals1$Upper[which(ConfidenceIntervals1$log_v==log(v_left))])/(log(v_right) - log(v_left)))
    if((upr_1_tmp < 0) == T){upr_1 <- Inf} else{upr_1 <- upr_1_tmp}
  }

  point_estimate_1 <- (log(v_right) - log(v_left))/(ConfidenceIntervals1$PointEstimate[which(ConfidenceIntervals1$log_v==log(v_right))] - ConfidenceIntervals1$PointEstimate[which(ConfidenceIntervals1$log_v==log(v_left))])

#  lwr_2 <- 1/(1/((max(ConfidenceIntervals2$Lower[which(ConfidenceIntervals2$loglog_v==log(log(v_left))):which(ConfidenceIntervals2$loglog_v==log(log(v_right)))]) - min(ConfidenceIntervals2$Upper[which(ConfidenceIntervals2$loglog_v==log(log(v_left))):which(ConfidenceIntervals2$loglog_v==log(log(v_right)))]))/(log(log(v_right)) - log(log(v_left)))) + 1)

  lwr_2 <- max(0, 1/((ConfidenceIntervals2$Upper[which(ConfidenceIntervals2$loglog_v==log(log(v_right)))] - ConfidenceIntervals2$Lower[which(ConfidenceIntervals2$loglog_v==log(log(v_left)))])/(log(log(v_right)) - log(log(v_left))) + 1))

#  upr_2 <- 1/(1/((max(ConfidenceIntervals2$Upper[which(ConfidenceIntervals2$loglog_v==log(log(v_left))):which(ConfidenceIntervals2$loglog_v==log(log(v_right)))]) - min(ConfidenceIntervals2$Lower[which(ConfidenceIntervals2$loglog_v==log(log(v_left))):which(ConfidenceIntervals2$loglog_v==log(log(v_right)))]))/(log(log(v_right)) - log(log(v_left)))) + 1)

  if(is.na(ConfidenceIntervals2$Lower[which(ConfidenceIntervals2$loglog_v==log(log(v_right)))])){warning("v_right is too large!"); upr_2 <- 1} else{
    upr_2 <- min(1/((ConfidenceIntervals2$Lower[which(ConfidenceIntervals2$loglog_v==log(log(v_right)))] - ConfidenceIntervals2$Upper[which(ConfidenceIntervals2$loglog_v==log(log(v_left)))])/(log(log(v_right)) - log(log(v_left))) + 1), 1)
    if((upr_2 < 0) == T){upr_2 <- 1}
  }
  if((lwr_2 > 1) == T){lwr_2 <- 0}

#  point_estimate_2 <- 1/(1/((ConfidenceIntervals2$PointEstimate[which(ConfidenceIntervals2$loglog_v==log(log(v_right)))] - ConfidenceIntervals2$PointEstimate[which(ConfidenceIntervals2$loglog_v==log(log(v_left)))])/(log(log(v_right)) - log(log(v_left)))) + 1)

  point_estimate_2 <- 1/((ConfidenceIntervals2$PointEstimate[which(ConfidenceIntervals2$loglog_v==log(log(v_right)))] - ConfidenceIntervals2$PointEstimate[which(ConfidenceIntervals2$loglog_v==log(log(v_left)))])/(log(log(v_right)) - log(log(v_left))) + 1)

  if(is.na(ConfidenceIntervals3$Lower[which(ConfidenceIntervals3$logloglog_v==log(log(log(v_right))))])){warning("v_right is too large!"); lwr_3 <- 0} else{
    lwr_3 <- max(0, (ConfidenceIntervals3$Lower[which(ConfidenceIntervals3$logloglog_v==log(log(log(v_right))))] - ConfidenceIntervals3$Upper[which(ConfidenceIntervals3$logloglog_v==log(log(log(v_left))))])/(log(log(log(v_right))) - log(log(log(v_left)))))
  }


  upr_3 <- max(0, (ConfidenceIntervals3$Upper[which(ConfidenceIntervals3$logloglog_v==log(log(log(v_right))))] - ConfidenceIntervals3$Lower[which(ConfidenceIntervals3$logloglog_v==log(log(log(v_left))))])/(log(log(log(v_right))) - log(log(log(v_left)))))


  point_estimate_3 <- (ConfidenceIntervals3$PointEstimate[which(ConfidenceIntervals3$logloglog_v==log(log(log(v_right))))] - ConfidenceIntervals3$PointEstimate[which(ConfidenceIntervals3$logloglog_v==log(log(log(v_left))))])/(log(log(log(v_right))) - log(log(log(v_left))))

  df_CI <- data.frame(Tail_Type = c("Power", "Sub-Exponential", "Near-Exponential"), CI_level = rep(C_Level, 3), lwr = c(lwr_1, lwr_2, lwr_3), upr = c(upr_1, upr_2, upr_3), PointEstimate = c(point_estimate_1, point_estimate_2, point_estimate_3))
  attr(df_CI, "title") <- paste0(C_Level*100, "% CI based on information from v_left = ",v_left," to v_right = ", v_right)
  class(df_CI) <- c("titled_df", class(df_CI))

  if(any(df_CI$upr < df_CI$lwr) == T){warning("The confidence intervals are not reliable. You may try with different v_left or v_right values.")}

  rs <- c(mean((lts-mean(lts))*(log(vf)-mean(log(vf))))/sqrt((mean(lts^2)-mean(lts)^2)*(mean((log(vf))^2)-mean(log(vf))^2)),mean((lts-mean(lts))*(log(log(vf))-mean(log(log(vf)))))/sqrt((mean(lts^2)-mean(lts)^2)*(mean((log(log(vf)))^2)-mean(log(log(vf)))^2)),mean((lts-mean(lts))*(log(log(log(vf)))-mean(log(log(log(vf))))))/sqrt((mean(lts^2)-mean(lts)^2)*(mean((log(log(log(vf))))^2)-mean(log(log(log(vf))))^2)))

  if(mean((lts-mean(lts))*(log(vf)-mean(log(vf))))/sqrt((mean(lts^2)-mean(lts)^2)*(mean((log(vf))^2)-mean(log(vf))^2))<=0 | diff(range(ConfidenceIntervals[,4][v_left:v_right])) < 0){ # this is trying to capture exponential
    list_print <- list(Plot = final_plot, Conclusion = "The data suggest an exponential or thinner tail.", CI = df_CI)
    out <- list(Results = list_print, Type = 0, CI_0 = ConfidenceIntervals, CI_1 = ConfidenceIntervals1, CI_2 = ConfidenceIntervals2, CI_3 = ConfidenceIntervals3, z1v = z1v, Rank = c(0, rank(-rs)), Pearson_r = rs)
    class(out) <- "function_output"
    return(out)
  } else{
    index <- which.max(rs)
    # if ((diff(range(z1v[v_left:v_right])) < 0.05) == T) {
    #   if(v_left < exp(point_estimate_1)){warning(paste0("v_left must be greater than exp(exp(theta_hat)) = ", exp(exp(point_estimate_1)), " if power tail is to be selected!"))}
    #   list_print <- list(Plot = final_plot, Conclusion = "The data suggest an exponential or thinner tail.", CI = df_CI)
    #   out <- list(Results = list_print, Type = 1, CI_0 = ConfidenceIntervals, CI_1 = ConfidenceIntervals1, CI_2 = ConfidenceIntervals2, CI_3 = ConfidenceIntervals3, z1v = z1v, Rank = c(0, rank(-rs)), Pearson_r = rs)
    #   class(out) <- "function_output"
    #   return(out)}
    if (index == 1) {
      if(v_left < exp(point_estimate_1)){warning(paste0("v_left must be greater than exp(exp(theta_hat)) = ", exp(exp(point_estimate_1)), " if power tail is to be selected!"))}
      list_print <- list(Plot = final_plot, Conclusion = "The data suggest a power decaying tail.", CI = df_CI)
      out <- list(Results = list_print, Type = 1, CI_0 = ConfidenceIntervals, CI_1 = ConfidenceIntervals1, CI_2 = ConfidenceIntervals2, CI_3 = ConfidenceIntervals3, z1v = z1v, Rank = rank(-rs), Pearson_r = rs)
      class(out) <- "function_output"
      return(out)}
    if (index == 2) {
      list_print <- list(Plot = final_plot, Conclusion = "The data suggest a sub-exponential decaying tail.", CI = df_CI)
      out <- list(Results = list_print, Type = 2, CI_0 = ConfidenceIntervals, CI_1 = ConfidenceIntervals1, CI_2 = ConfidenceIntervals2, CI_3 = ConfidenceIntervals3, z1v = z1v, Rank = rank(-rs), Pearson_r = rs)
      class(out) <- "function_output"
      return(out)}
    if (index == 3) {
      list_print <- list(Plot = final_plot, Conclusion = "The data suggest a near-exponential decaying tail.", CI = df_CI)
      out <- list(Results = list_print, Type = 3, CI_0 = ConfidenceIntervals, CI_1 = ConfidenceIntervals1, CI_2 = ConfidenceIntervals2, CI_3 = ConfidenceIntervals3, z1v = z1v, Rank = rank(-rs), Pearson_r = rs)
      class(out) <- "function_output"
      return(out)}
  }
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("v", "Lower", "Upper")) # Silent false note

#' @exportS3Method
print.multiplot <- function(x, ...) {
  print(cowplot::plot_grid(x[[1]], x[[2]], x[[3]], x[[4]]))
}
#' @exportS3Method
print.titled_df <- function(x, ...){
  # Print the title and the data frame
  print(attr(x, "title"))
  print.data.frame(x)
}
#' @exportS3Method
print.function_output <- function(x, ...) {
  print(x$Result)  # Only print Result by default
}

