# The function is developed under R/4.1.0 Platform: x86_64-apple-darwin17.0 (64-bit).
# Author/Maintainer: Jialin Zhang (JZ); jzhang (at) math.msstate.edu

# TailClassifier is the Tail-Classifier function. The function suggests one of the following types of tail for your discrete data: 1) Power decaying tail; 2) Sub-exponential decaying tail; and 3) Near-exponential decaying tail.

## Guidance:

#' @title
#' Tail Classifier for Thick-Tailed Discrete Data
#' @description
#' Function TailClassifier() in this package is a Tail-Classifier function. The function suggests one of the following types of tail for your discrete data: 1) Power decaying tail; 2) Sub-exponential decaying tail; and 3) Near-exponential decaying tail.
#'
#' @param sample_frequencies The frequency counts for your discrete sample data.
#' @param v.left The starting point of tail profile. 5 is recommended. A smaller v.left may lead to unreliable results. A larger v.left might be adopted if the sample size is extremely large.
#' @param v.right The ending point of tail profile. Recommendation is 5\% of the sample size but no more than 500. For example, a sample with size 1000 could choose v.right to be 50; and a sample with size 20000 could choose v.right to be 500.
#' @return A statement on the type of tail.
#' @examples
#' ## read built-in random sample that was generated under a sub-exponential distribution
#' csv <- system.file("extdata", "sample_data.csv", package = "TailClassifier")
#' sample_data <- readr::read_csv(csv)
#' ## generate the frequency table of the sample
#' sample_freq=table(sample_data)
#' ## make a classification
#' TailClassifier(sample_freq)

#' @export

TailClassifier <- function(sample_frequencies, v.left=5, v.right=min(floor(sum(sample_frequencies)/20), 500)) {
  n=sum(sample_frequencies)
  if (is.integer(sample_frequencies) == FALSE) {stop("Your input must be sample frequencies (counts)!")}
  if (floor(v.left) != v.left) {stop("v.left must be an integer!")}
  if (floor(v.right) != v.right) {stop("v.right must be an integer!")}
  if (v.left < 4) {warning("Argument 'v.left' should be at least 4!")}
  if (v.right <= v.left) {stop("Argument 'v.right' should be much larger than 'v.left'!")}
  if (v.right > n-1) {stop("Argument 'v.right' should be less than n (sample size)!")}
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
  z1v = z1f(sample_frequencies)
  zeros = sum(z1v==0)
  z1v = z1v[z1v!=0]
  tv = z1v*(0:(n-1-zeros))
  vf = ((v.left+1):(v.right-zeros+1))
  ts = tv[vf]
  lts = log(ts)
  if(mean((lts-mean(lts))*(log(vf)-mean(log(vf))))/sqrt((mean(lts^2)-mean(lts)^2)*(mean((log(vf))^2)-mean(log(vf))^2))<0){
    return("The data suggest an exponentially decaying tail.")
  } else{
    rs=c(mean((lts-mean(lts))*(log(vf)-mean(log(vf))))/sqrt((mean(lts^2)-mean(lts)^2)*(mean((log(vf))^2)-mean(log(vf))^2)),mean((lts-mean(lts))*(log(log(vf))-mean(log(log(vf)))))/sqrt((mean(lts^2)-mean(lts)^2)*(mean((log(log(vf)))^2)-mean(log(log(vf)))^2)),mean((lts-mean(lts))*(log(log(log(vf)))-mean(log(log(log(vf))))))/sqrt((mean(lts^2)-mean(lts)^2)*(mean((log(log(log(vf))))^2)-mean(log(log(log(vf))))^2)))
    index <- which.max(rs)
    if (index == 1) {return("The data suggest a power decaying tail.")}
    if (index == 2) {return("The data suggest a sub-exponential decaying tail.")}
    if (index == 3) {return("The data suggest a near-exponential decaying tail.")}
  }
}

