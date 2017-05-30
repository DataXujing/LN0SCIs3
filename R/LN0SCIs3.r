
#------------LN0SCIs3------------------------------

###############################

##  By  Jing Xu & Xinmin Li & HuaLiang##

##  ----UTF-8--------------##

##############################

#-------------------------------------------------------GPQW--------------------------------------
GPQW0<-function(n1,n2,n3,p1,p2,p3,mu1,mu2,mu3,sigma1,sigma2,sigma3,alpha,N,C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))  )                
{
  t1<-Sys.time()
  result<-list();
  
  n10<-rbinom(1,n1,p1);
  n20<-rbinom(1,n2,p2);
  n30<-rbinom(1,n3,p3);
  
  n11<-n1-n10;
  n21<-n2-n20;
  n31<-n3-n30;
  
  y1bar<-rnorm(1,mu1,(sigma1/sqrt(n11)));
  y2bar<-rnorm(1,mu2,(sigma2/sqrt(n21)));
  y3bar<-rnorm(1,mu3,(sigma3/sqrt(n31)));
  
  s1sq<-sigma1^2*rchisq(1,df=(n11-1),ncp=0);
  s2sq<-sigma2^2*rchisq(1,df=(n21-1),ncp=0);
  s3sq<-sigma3^2*rchisq(1,df=(n31-1),ncp=0);
  
  Z1<-rnorm(N,0,1);
  Z2<-rnorm(N,0,1);
  Z3<-rnorm(N,0,1);
  U1<-rchisq(N,df=(n11-1),ncp=0);
  U2<-rchisq(N,df=(n21-1),ncp=0);
  U3<-rchisq(N,df=(n31-1),ncp=0);
  
  Z1W<-rnorm(N,0,1);
  Z2W<-rnorm(N,0,1);
  Z3W<-rnorm(N,0,1);
  
  #C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))   
  
  TPW1=TPW2=TPW3=TW1=TW2=TW3=numeric(N);
  r=r1=r2=matrix(0,choose(3,2),N);
  
  for(k in 1:N)
  {
    TPW1[k]<-(n10+0.5*(Z1W[k])^2)/(n1+(Z1W[k])^2)-(Z1W[k]*sqrt(n10*(1-n10/n1)+(Z1W[k])^2/4))/(n1+(Z1W[k])^2);
    TPW2[k]<-(n20+0.5*(Z2W[k])^2)/(n2+(Z2W[k])^2)-(Z2W[k]*sqrt(n20*(1-n20/n2)+(Z2W[k])^2/4))/(n2+(Z2W[k])^2);
    TPW3[k]<-(n30+0.5*(Z3W[k])^2)/(n3+(Z3W[k])^2)-(Z3W[k]*sqrt(n30*(1-n30/n3)+(Z3W[k])^2/4))/(n3+(Z3W[k])^2);
    
    TW1[k]<-log(1-TPW1[k])+(y1bar-Z1[k]*sqrt(s1sq)/(sqrt(n11*U1[k]))+s1sq/(2*U1[k]));
    TW2[k]<-log(1-TPW2[k])+(y2bar-Z2[k]*sqrt(s2sq)/(sqrt(n21*U2[k]))+s2sq/(2*U2[k]));
    TW3[k]<-log(1-TPW3[k])+(y3bar-Z3[k]*sqrt(s3sq)/(sqrt(n31*U3[k]))+s3sq/(2*U3[k]));
    
    r[,k]<-C2%*%c(TW1[k],TW2[k],TW3[k])
    
  }
  
  r1[1,]<-sort(r[1,]);               
  r1[2,]<-sort(r[2,]);   
  r1[3,]<-sort(r[3,]);   
  
  r2[1,]<-rank(r[1,]);
  r2[2,]<-rank(r[2,]);
  r2[3,]<-rank(r[3,]);
  
  mik1=mak1=numeric(N);
  for(k in 1:N)
  {
    mik1[k]<-min(c(r2[1,k],r2[2,k],r2[3,k]));
    mak1[k]<-max(c(r2[1,k],r2[2,k],r2[3,k]));
  }
  
  mik2<-sort(mik1);
  mak2<-sort(mak1);
  
  kl<-mik2[N*alpha/2]
  ku<-mak2[N*(1-alpha/2)]
  
  
  the.first.lower.limit<-round(r1[1,kl],6);
  the.first.upper.limit<-round(r1[1,ku],6);

  the.second.lower.limit<-round(r1[2,kl],6);
  the.second.upper.limit<-round(r1[2,ku],6);

  the.third.lower.limit<-round(r1[3,kl],6);
  the.third.upper.limit<-round(r1[3,ku],6);
  
  result$title1<-"====================Method: GPQW===================="
  result$title2<-"The Simultaneous Confidence Intervals are:          "

  interval0<--data.frame(matrix(c(round(r1[1,kl],6),round(r1[1,ku],6),
                                  round(r1[2,kl],6),round(r1[2,ku],6),
                                  round(r1[3,kl],6),round(r1[3,ku],6)),3,2,byrow=T))
  names(interval0)<-c("LCL","UCL")

  
  result$interval<-interval0
  
  
  t2<-Sys.time()
  result$star<-'**********************Time**************************'
  result$t=t2-t1;
  result;
  
}
GPQW<-function(n1,n2,n3,p1,p2,p3,mu1,mu2,mu3,sigma1,sigma2,sigma3,alpha,N,C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))  ) 
{
  GPQW1<-GPQW0(n[1],n[2],n[3],p[1],p[2],p[3],mu[1],mu[2],mu[3],sigma[1],sigma[2],sigma[3],alpha,N);
  for(i in c('title1','title2','interval', 'star',"t"))
  {
    print(GPQW1[[i]])
    
  }
  
}
  






#----------------GPQH -----------------------------------------------------

rndMixture<-function(N0,N1,p,N)
{
  u<-numeric(N);
  u<-runif(N,0,1);
  r<-rep(0,N);
  for(i in 1:N)
  {
    if(u[i]<p)
    {r[i]<-rbeta(1,N0+1,N1);}
    else{r[i]<-rbeta(1,N0,N1+1);}
  }
  r          
}      
#---------
GPQH0<-function(n1,n2,n3,p1,p2,p3,mu1,mu2,mu3,sigma1,sigma2,sigma3,alpha,N,C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))   )               
{
  
  t1<-Sys.time()
  result<-list();
  
  n10<-rbinom(1,n1,p1);
  n20<-rbinom(1,n2,p2);
  n30<-rbinom(1,n3,p3);
  
  n11<-n1-n10;
  n21<-n2-n20;
  n31<-n3-n30;
  
  y1bar<-rnorm(1,mu1,(sigma1/sqrt(n11)));
  y2bar<-rnorm(1,mu2,(sigma2/sqrt(n21)));
  y3bar<-rnorm(1,mu3,(sigma3/sqrt(n31)));
  
  s1sq<-sigma1^2*rchisq(1,df=(n11-1),ncp=0);
  s2sq<-sigma2^2*rchisq(1,df=(n21-1),ncp=0);
  s3sq<-sigma3^2*rchisq(1,df=(n31-1),ncp=0);
  
  Z1<-rnorm(N,0,1);
  Z2<-rnorm(N,0,1);
  Z3<-rnorm(N,0,1);
  U1<-rchisq(N,df=(n11-1),ncp=0);
  U2<-rchisq(N,df=(n21-1),ncp=0);
  U3<-rchisq(N,df=(n31-1),ncp=0);
  
  
  TPH1<-rndMixture(n10,n11,0.5,N);
  TPH2<-rndMixture(n20,n21,0.5,N);
  TPH3<-rndMixture(n30,n31,0.5,N);
  
  #C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))   
  
  TPW1=TPW2=TPW3=TW1=TW2=TW3=numeric(N);
  r=r1=r2=matrix(0,choose(3,2),N);
  
  for(k in 1:N)
  {
    
    TW1[k]<-log(1-TPH1[k])+(y1bar-Z1[k]*sqrt(s1sq)/(sqrt(n11*U1[k]))+s1sq/(2*U1[k]));
    TW2[k]<-log(1-TPH2[k])+(y2bar-Z2[k]*sqrt(s2sq)/(sqrt(n21*U2[k]))+s2sq/(2*U2[k]));
    TW3[k]<-log(1-TPH3[k])+(y3bar-Z3[k]*sqrt(s3sq)/(sqrt(n31*U3[k]))+s3sq/(2*U3[k]));
    
    r[,k]<-C2%*%c(TW1[k],TW2[k],TW3[k])
    
  }
  
  r1[1,]<-sort(r[1,]);               
  r1[2,]<-sort(r[2,]);   
  r1[3,]<-sort(r[3,]);   
  
  r2[1,]<-rank(r[1,]);
  r2[2,]<-rank(r[2,]);
  r2[3,]<-rank(r[3,]);
  
  mik1=mak1=numeric(N);
  for(k in 1:N)
  {
    mik1[k]<-min(c(r2[1,k],r2[2,k],r2[3,k]));
    mak1[k]<-max(c(r2[1,k],r2[2,k],r2[3,k]));
  }
  
  mik2<-sort(mik1);
  mak2<-sort(mak1);
  
  kl<-mik2[N*alpha/2]
  ku<-mak2[N*(1-alpha/2)]
  
  # result$the.first.lower.limit<-r1[1,kl];
  # result$the.first.upper.limit<-r1[1,ku];
  # 
  # result$the.second.lower.limit<-r1[2,kl];
  # result$the.second.upper.limit<-r1[2,ku];
  # 
  # result$the.third.lower.limit<-r1[3,kl];
  # result$the.third.upper.limit<-r1[3,ku];
  
  result$title1<-"====================Method: GPQH===================="
  result$title2<-"The Simultaneous Confidence Intervals are:          "
  
  interval0<--data.frame(matrix(c(round(r1[1,kl],6),round(r1[1,ku],6),
                                  round(r1[2,kl],6),round(r1[2,ku],6),
                                  round(r1[3,kl],6),round(r1[3,ku],6)),3,2,byrow=T))
  names(interval0)<-c("LCL","UCL")
  
  
  result$interval<-interval0
  
  
  t2<-Sys.time()
  result$star<-'**********************Time**************************'
  result$t=t2-t1;
  result;

  
  
}


GPQH<-function(n1,n2,n3,p1,p2,p3,mu1,mu2,mu3,sigma1,sigma2,sigma3,alpha,N,C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))  ) 
{
  GPQH1<-GPQH0(n[1],n[2],n[3],p[1],p[2],p[3],mu[1],mu[2],mu[3],sigma[1],sigma[2],sigma[3],alpha,N);
  for(i in c('title1','title2','interval', 'star',"t"))
  {
    print(GPQH1[[i]])
    
  }
  
}




#---------------------GPQB_H----------------------
rndMixture<-function(N0,N1,p,N)
{
  u<-numeric(N);
  u<-runif(N,0,1);
  r<-rep(0,N);
  for(i in 1:N)
  {
    if(u[i]<p)
    {r[i]<-rbeta(1,N0+1,N1);}
    else{r[i]<-rbeta(1,N0,N1+1);}
  }
  r          
}     


GPQB_H0<-function(n1,n2,n3,p1,p2,p3,mu1,mu2,mu3,sigma1,sigma2,sigma3,alpha,N,C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))    )                
{
  t1<-Sys.time()
  result<-list();
  
  n10<-rbinom(1,n1,p1);
  n20<-rbinom(1,n2,p2);
  n30<-rbinom(1,n3,p3);
  
  n11<-n1-n10;
  n21<-n2-n20;
  n31<-n3-n30;
  
  y1bar<-rnorm(1,mu1,(sigma1/sqrt(n11)));
  y2bar<-rnorm(1,mu2,(sigma2/sqrt(n21)));
  y3bar<-rnorm(1,mu3,(sigma3/sqrt(n31)));
  
  s1sq<-sigma1^2*rchisq(1,df=(n11-1),ncp=0);
  s2sq<-sigma2^2*rchisq(1,df=(n21-1),ncp=0);
  s3sq<-sigma3^2*rchisq(1,df=(n31-1),ncp=0);
  
  Z1<-rnorm(N,0,1);
  Z2<-rnorm(N,0,1);
  Z3<-rnorm(N,0,1);
  U1<-rchisq(N,df=(n11-1),ncp=0);
  U2<-rchisq(N,df=(n21-1),ncp=0);
  U3<-rchisq(N,df=(n31-1),ncp=0);
  
  
  TPH1<-rndMixture(n10,n11,0.5,N);
  TPH2<-rndMixture(n20,n21,0.5,N);
  TPH3<-rndMixture(n30,n31,0.5,N);
  
  #C2<-rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))   
  
  TW1=TW2=TW3=numeric(N);
  r=r1=matrix(0,choose(3,2),N);
  
  for(k in 1:N)
  {
    
    TW1[k]<-log(1-TPH1[k])+(y1bar-Z1[k]*sqrt(s1sq)/(sqrt(n11*U1[k]))+s1sq/(2*U1[k]));
    TW2[k]<-log(1-TPH2[k])+(y2bar-Z2[k]*sqrt(s2sq)/(sqrt(n21*U2[k]))+s2sq/(2*U2[k]));
    TW3[k]<-log(1-TPH3[k])+(y3bar-Z3[k]*sqrt(s3sq)/(sqrt(n31*U3[k]))+s3sq/(2*U3[k]));
    
    r[,k]<-C2%*%c(TW1[k],TW2[k],TW3[k])
    
  }
  
  r1[1,]<-sort(r[1,]);               
  r1[2,]<-sort(r[2,]);   
  r1[3,]<-sort(r[3,]);   
  
  # result$the.first.lower.limit<-r1[1,N*alpha/(2*2)];
  # result$the.first.upper.limit<-r1[1,N*(1-alpha/(2*2))];
  # 
  # result$the.second.lower.limit<-r1[2,N*alpha/(2*3)];
  # result$the.second.upper.limit<-r1[2,N*(1-alpha/(2*3))];
  # 
  # result$the.third.lower.limit<-r1[3,N*alpha/(2*2)];
  # result$the.third.upper.limit<-r1[3,N*(1-alpha/(2*2))];
  
  result$title1<-"===================Method: GPQB_H==================="
  result$title2<-"The Simultaneous Confidence Intervals are:          "
  
  interval0<--data.frame(matrix(c(round(r1[1,N*alpha/(2*2)],6),round(r1[1,N*(1-alpha/(2*2))],6),
                                  round(r1[2,N*alpha/(2*3)],6),round(r1[2,N*(1-alpha/(2*3))],6),
                                  round(r1[3,N*alpha/(2*2)],6),round(r1[3,N*(1-alpha/(2*2))],6)),3,2,byrow=T))
  names(interval0)<-c("LCL","UCL")
  
  
  result$interval<-interval0
  
  
  t2<-Sys.time()
  result$star<-'**********************Time**************************'
  result$t=t2-t1;
  result;
  
  
  
  
}


GPQB_H<-function(n1,n2,n3,p1,p2,p3,mu1,mu2,mu3,sigma1,sigma2,sigma3,alpha,N,C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))  ) 
{
  GPQB_H1<-GPQB_H0(n[1],n[2],n[3],p[1],p[2],p[3],mu[1],mu[2],mu[3],sigma[1],sigma[2],sigma[3],alpha,N);
  for(i in c('title1','title2','interval', 'star',"t"))
  {
    print(GPQB_H1[[i]])
    
  }
  
}



