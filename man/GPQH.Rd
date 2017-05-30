\name{GPQH}
\alias{GPQH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  GPQH
}
\description{
   A function to compute SCIs for lon-normal populations with zeros.
}
\usage{
GPQH(n1, n2, n3, p1, p2, p3, mu1, mu2, mu3, sigma1, sigma2, sigma3, alpha, N, C2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ni}{
     the ith sample size,must be an integer
}

  \item{pi}{
    Sample i's zero probability
}

  \item{mui}{
     Non-zeros sample i's mean after log transformation
}

  \item{sigmai}{
     Non-zeros sample i's variance after log transformation
}

  \item{N}{
    The number of independent generated data sets
}

  \item{C2}{
     Matrix C, You can refer to reference [6] for specific forms
}

}
\details{
%%  Non
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% Non
}
\author{
  Jing Xu, Xinmin Li, HuaLiang
}
\note{
%%  Non
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

#-----GPQH Example-------------------------------

alpha<-0.05
N<-5000;
M<-5000;

p<-c(0.1,0.1,0.1)
n<-c(10,10,10)
mu<-c(0,0,0)
sigma<-c(1,1,1)

GPQH(n[1],n[2],n[3],p[1],p[2],p[3],mu[1],mu[2],mu[3],sigma[1],sigma[2],sigma[3],alpha,N);
}
