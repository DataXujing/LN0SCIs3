\name{LN0SCIs3-package}
\alias{LN0SCIs3-package}
\alias{LN0SCIs3}
\docType{package}
\title{
\packageTitle{LN0SCIs3}
}
\description{
\packageDescription{LN0SCIs3}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{LN0SCIs3}
\packageIndices{LN0SCIs3}
 At present, these three function perform better than other methods that can be used to calculate the simultaneous confidence interval of log-normal populations with excess zeros.
}
\author{
\packageAuthor{LN0SCIs3}

Maintainer: \packageMaintainer{LN0SCIs3}
}
\references{
  [1] Besag I, Green P, Higdon D, Mengersen K, 1995. Bayesian computation and Stochastic-
      systems.

  [2] Hannig J, Abdel-Karim A, Iyer H, 2006. Simultaneous fiducial generalized confidence
      intervals for ratios of means of lognormal distribution.

  [3] Hannig J, Lee T C M, 2009.Generalized fiducial inference for wavelet regression.

  [4] Li X, Zhou X, Tian L, 2013. Interval estimation for the mean of lognormal data with
      excess zeros.

  [5] Schaarschmidt F, 2013. Simultaneous confidence intervals for multiple comparisons a-
      mong expected values of log-normal variables.

  [6] Jing Xu, Xinmin Li, Hua Liang. Simultaneous Confidence Intervals for Ratios of Means of       Log-normal Populations with Zeros.


}
\examples{
 #-----GPQW Example---------------------------------
alpha<-0.05
N<-5000;
M<-5000;

p<-c(0.1,0.15,0.1)
n<-c(30,15,10)
mu<-c(1,1.3,2)
sigma<-c(1,1,1)

GPQW(n[1],n[2],n[3],p[1],p[2],p[3],mu[1],mu[2],mu[3],sigma[1],sigma[2],sigma[3],alpha,N);


#-----GPQH Example-------------------------------

alpha<-0.05
N<-5000;
M<-5000;

p<-c(0.1,0.1,0.1)
n<-c(10,10,10)
mu<-c(0,0,0)
sigma<-c(1,1,1)

GPQH(n[1],n[2],n[3],p[1],p[2],p[3],mu[1],mu[2],mu[3],sigma[1],sigma[2],sigma[3],alpha,N);


#-------------GPQB_H Example-----------------------
alpha<-0.05
N<-5000;
M<-5000;

p<-c(0.1,0.1,0.1)
n<-c(10,10,10)
mu<-c(0.5,0.5,0.5)
sigma<-c(1,1,1)

GPQB_H(n[1],n[2],n[3],p[1],p[2],p[3],mu[1],mu[2],mu[3],sigma[1],sigma[2],sigma[3],alpha,N);


}
