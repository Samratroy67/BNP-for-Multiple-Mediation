#==== Librarys needed======

library(MASS)
library(mvnfast)
library(invgamma)
#setwd("C:/Study_2023/BNP/Finding_data")
source("Cluster_Complex_Setup.R")
data <- read.csv("data_for_sec8.csv")

#=====Define mathematical functions=====

expit <- function(x) {1/(1+exp(-x)) }
logit <- function(p) {log(p)-log(1-p)}

inv<-function(x){chol2inv(chol(x))}

dinvchisq <- function (x, df , scale = 1/df,log = FALSE)
{
  x <- as.vector(x)
  df <- as.vector(df)
  scale <- as.vector(scale)
  if (any(x <= 0)) 
    stop("x must be positive.")
  if (any(df <= 0)) 
    stop("The df parameter must be positive.")
  if (any(scale <= 0)) 
    stop("The scale parameter must be positive.")
  NN <- max(length(x), length(df), length(scale))
  x <- rep(x, len = NN)
  df <- rep(df, len = NN)
  scale <- rep(scale, len = NN)
  nu <- df/2
  dens <- nu * log(nu) - log(gamma(nu)) + nu * log(scale) - 
    (nu + 1) * log(x) - (nu * scale/x)
  if (log == FALSE) 
    dens <- exp(dens)
  dens
}

attr(dinvchisq,"ex") <- function(){
  x <- dinvchisq(1,1,1)
  x <- rinvchisq(10,1)
  
  #Plot Probability Functions
  x <- seq(from=0.1, to=5, by=0.01)
  plot(x, dinvchisq(x,0.5,1), ylim=c(0,1), type="l", main="Probability Function",
       ylab="density", col="red")
  lines(x, dinvchisq(x,1,1), type="l", col="green")
  lines(x, dinvchisq(x,5,1), type="l", col="blue")
  legend(3, 0.9, expression(paste(nu==0.5, ", ", lambda==1),
                            paste(nu==1, ", ", lambda==1), paste(nu==5, ", ", lambda==1)),
         lty=c(1,1,1), col=c("red","green","blue"))	
}


rinvchisq <- function (n , df , scale=1/df)
{
  df <- rep(df, len = n)
  scale <- rep(scale, len = n)
  if (any(df <= 0)) 
    stop("The df parameter must be positive.")
  if (any(scale <= 0)) 
    stop("The scale parameter must be positive.")
  x <- (df * scale)/rchisq(n, df = df)
  return(x)
}

#================== define function to Generate the data===========

#====================================== Functions to update parameters==============================

#(1)==regression variance update=========

newregvar<-function(tempx,tempy, prec0, beta0){
  an<-betaa0+length(tempy)/2
  ifelse((ncol(as.matrix(tempx))>1),(xtx<-t(tempx)%*%tempx),(xtx<-as.matrix(tempx)%*%t(as.matrix(tempx))))
  newprec<-xtx+prec0
  if(ncol(as.matrix(tempx))>1){betan<-(ginv(xtx+prec0))%*%(prec0%*%beta0+t(tempx)%*%tempy)}
  if(ncol(as.matrix(tempx))<=1){betan<-(ginv(xtx+prec0))%*%(prec0%*%beta0+(tempx*tempy))}
  bn<-betab0+.5*(sum(tempy*tempy)+t(beta0)%*%prec0%*%beta0-t(betan)%*%newprec%*%betan)
  return(rinvgamma(1,an,bn))
}

#(2)==regression betas update===========

newbet<-function(tempx,tempy,sig, prec0, beta0){
  ifelse((ncol(as.matrix(tempx))>1),(xtx<-t(tempx)%*%tempx),(xtx<-as.matrix(tempx)%*%t(as.matrix(tempx))))
  newprec<-xtx+prec0
  if(ncol(as.matrix(tempx))>1){betan<-(ginv(xtx+prec0))%*%(prec0%*%beta0+t(tempx)%*%tempy)}
  if(ncol(as.matrix(tempx))<=1){betan<-(ginv(xtx+prec0))%*%(prec0%*%beta0+(tempx*tempy))}
  return(rmvn(1,betan,sig*inv(newprec)))
}


#(3)===================================
updatevar<-function(tempx){
  newdf<-nu0+length(tempx)
  if(length(tempx)==1){varx<-0}
  if(length(tempx)>1){varx<-var(tempx)}
  numer<-nu0*tau0+(length(tempx)-1)*varx+(c0*length(tempx)/(c0+length(tempx)))*(mean(tempx)-mu0)^2
  newval<-rinvchisq(1,newdf,numer/newdf)
  return(newval)
}

#==(4)==================================
updatemean<-function(tempx,curtau){
  newvar<-1/(c0/curtau+length(tempx)/curtau)
  newmean<-(mu0*c0/curtau+mean(tempx)*length(tempx)/curtau)*newvar
  newval<-rnorm(1,newmean,sqrt(newvar))
  return(newval)
}

#==(5)==================================

newalpthet<-function(numclus,prioralp){
  eta<-rbeta(1,prioralp+1,n)
  pieta<-(numclus/(n*(1-log(eta))))/(1+numclus/(n*(1-log(eta))))
  whichmix<-rbinom(1,1,pieta)
  newalp<-whichmix*rgamma(1,(alpa0+numclus),(alpb0-log(eta)))+
    (1-whichmix)*rgamma(1,(alpa0+numclus-1),(alpb0-log(eta)))
  return(newalp)
}

#==(6)==================================

#update alpha_psi using MH

newalppsi<-function(alpha){
  sortuniquey<-sort(unique(s[,1]))
  ss<-numeric(length(sortuniquey))
  for(j in 1:length(sortuniquey)){
    ss[j]<-sum(s[,1]==sortuniquey[j])
  }
  k1 <- length(ss)
  likecur<-dgamma(alpha,alpa0,alpb0)*(alpha^(nrow(unique(s[,c(1:2)]))-k1))*prod((alpha+ss)*beta(alpha+1,ss))
  proposed<-rgamma(1,2,1)
  likeprop<-dgamma(proposed,alpa0,alpb0)*(proposed^(nrow(unique(s[,c(1:2)]))-k1))*prod((proposed+ss)*beta(proposed+1,ss))
  rat<-(likeprop)/(likecur)
  if(runif(1,0,1)<rat){alpha<-proposed}
  return(alpha)
}

#==(7)==================================

#update alpha_omega using MH

newalpomega <-function(alpha){
  sortuniqueym<-unique(s[,c(1:2)])
  ss<-numeric(nrow(sortuniqueym))
  for(j in 1:nrow(sortuniqueym)){
    ss[j]<-sum(s[,1]==sortuniqueym[j,1] & s[,2]==sortuniqueym[j,2])
  }
  likecur<-dgamma(alpha,alpa0,alpb0)*(alpha^(nrow(unique(s))-nrow(sortuniqueym)))*prod((alpha+ss)*beta(alpha+1,ss))
  proposed<-rgamma(1,2,1)
  likeprop<-dgamma(proposed,alpa0,alpb0)*(proposed^(nrow(unique(s))-nrow(sortuniqueym)))*prod((proposed+ss)*beta(proposed+1,ss))
  rat<-(likeprop)/(likecur)
  if(runif(1,0,1)<rat){alpha<-proposed}
  return(alpha)
}

#===============================Process the data=============================================

n<- 1701
p1<-1
p2<-12
q <- 2
ptx<-1
p<-1+ptx+p1+p2+q
pp <- 1+ptx+p1+p2
#====================== Here we specify seed====================


Y <- data$CESD
M <- data[,c(18,20)]
A <- data$trt
L <- data[,3:15]

#standardize continuous L and M
L[,(p1+1):(p1+p2)]<-scale(L[,(p1+1):(p1+p2)])
M[,1:q]<-scale(M[,1:q])

X<-cbind(A,L)
matX<-cbind(1,X)
full_cov <- cbind(M,A,L)
full_mat <- cbind(1,full_cov)
full_cov <- data.matrix(full_cov, rownames.force = NA)
full_mat <- data.matrix(full_mat, rownames.force = NA)
X <- data.matrix(X, rownames.force = NA)
L <- data.matrix(L, rownames.force = NA)
matX <- data.matrix(matX, rownames.force = NA)
M <- data.matrix(M, rownames.force = NA)

reg1 <- lm(Y~full_cov)
bs <- reg1$coefficients
as <- vcov(reg1)
ps<-inv(n/5*(as)) 

bs1 <- list()
ps1 <- list()
for (r in 1:q)
{
  reg2 <- lm(M[,r]~X)
  bs1[[r]] <- reg2$coefficients
  as1 <- vcov(reg2)
  ps1[[r]]<-inv(n/5*(as1))
}

betaa0<-1 #hyperprior regression variance parameter
betab0<-1 #other hyperprior regression variance parameter

#======================================= Initial Clusters==============================

#initially put everyone in one of the 3*2*2 =12 clusters

s<-matrix(nrow=n,ncol=3,1)

distymx<-dist(cbind(Y,M,A,L))
hi<-hclust(distymx,method="ward.D")
s[,1]<-cutree(hi,k=3)

for(j in 1:3){
  distmx<-dist(cbind(M[s[,1]==j,],A[s[,1]==j],L[s[,1]==j,]))
  hi<-hclust(distmx,method="ward.D")
  groups<-cutree(hi,k=2)
  s[s[,1]==j,2]<-groups
}
for(j in 1:3){
  for ( k in 1:2)
  {
    distx<-dist(cbind(A[s[,1]==j & s[,2]==k],L[s[,1]==j & s[,2]==k,]))
    hi<-hclust(distx,method="ward.D")
    groups<-cutree(hi,k=2)
    if(length(A[s[,1]==j & s[,2]==k])==1)
    {
      groups <- 1
    }
    s[s[,1]==j & s[,2]==k,3]<-groups
  }
}
s[1,] <- c(3,2,2)

#======================update values of beta_y's for outcome regressions==================

k<-length(unique(s[,1])) #number of y clusters
sig2_y<-numeric(k)
beta_y<-matrix(nrow=k,ncol=(p))
for(i in 1:k)
{
  sig2_y[i]<-newregvar(full_mat[s[,1]==unique(s[,1])[i],],Y[s[,1]==unique(s[,1])[i]], ps, bs)
  beta_y[i,]<-newbet(full_mat[s[,1]==unique(s[,1])[i],],Y[s[,1]==unique(s[,1])[i]],sig2_y[i], ps, bs)
}
sig2_y <- t(rbind(seq(1,k,1), sig2_y))
temp <- matrix(seq(1,k,1))
beta_y <- cbind(temp, beta_y)


#============update values of beta_m 's for mediator regressions============

l <- length(unique(s[,2]))

final <- numeric()
final1 <- numeric()
for ( r in 1:q)
{
  for ( i in 1: k)
  {
    for ( j in 1:l)
    {
      temp1 <- newregvar(matX[s[,1]==unique(s[,1])[i] & s[,2]==unique(s[,2])[j],],M[s[,1]==unique(s[,1])[i] & s[,2]==unique(s[,2])[j],r], ps1[[r]], bs1[[r]])
      temp2 <- newbet(matX[s[,1]==unique(s[,1])[i] & s[,2]==unique(s[,2])[j],],M[s[,1]==unique(s[,1])[i] & s[,2]==unique(s[,2])[j],r],temp1, ps1[[r]], bs1[[r]])
      final <- rbind(final, c(r,i,j,temp1))
      final1 <- rbind(final1, c(r,i,j,temp2))
    }
  }
}
beta_m <- final1
sig_m <- final

#=========update covariate parameters within each cluster============
#tauinit<-1
a0 <- 1
b0 <- 1
nu0<-2
tau0<-1
c0<-0.5
mu0<-0
alpa0<-1
alpb0<-1

nk<-nrow(unique(s))
m <- length(unique(s[,3]))
xpipars<-matrix(nrow=nk,ncol=(p1+ptx)) #p1 L's +treatment 
xmupars<-matrix(nrow=nk,ncol=p2) 
xsigpars<-matrix(nrow=nk,ncol=p2)  


# ===update binary predictor parameters====

for(i in 1:(ptx+p1)){
  #beta(1,1) prior
  count<-1
  for(j in 1:k){
    
    for(t in 1:l){
      
      for ( u in 1: m){
        
        tempx <- X[s[,1]==j & s[,2]==t & s[,3]==u,i]
        xpipars[count,i] <- rbeta(1,(sum(tempx)+a0),(length(tempx)-sum(tempx)+b0))
        count <- count+1
      }
    }
  }
}

#=====update parameters for continuous predictors=======

for(i in 1:(p2)){
  #beta(1,1) prior
  count<-1
  for(j in 1:k){
    
    for(t in 1:l){
      
      for ( u in 1: m){
        
        tempx<-X[s[,1]==j & s[,2]==t & s[,3]==u,(p1+ptx+i)]
        xsigpars[count,i]<-updatevar(tempx)
        xmupars[count,i]<-updatemean(tempx,xsigpars[count,i])
        count <- count+1
      }
    }
  }
}

uniques <- unique(s)
uniqueS <- uniques[order(uniques[,1],uniques[,2], uniques[,3]), ,drop=FALSE]

#======== Add uniqueS to the covariate parameters=======
xpipars <- cbind(uniqueS,xpipars)
xmupars <- cbind(uniqueS,xmupars)
xsigpars <- cbind(uniqueS,xsigpars)

#========= h0x calculation =============================

h0x<-matrix(nrow=n,ncol=ncol(X),0)

denom <- (sqrt(2*pi)/sqrt(c0))*gamma(nu0/2)*(2/(nu0*tau0))^(nu0/2)
cn <- c0 + 1
nun <- nu0 + 1
for(i in 1:n){
  for(j in 1:p2){
    x<-X[i,(p1+ptx+j)]
    taun <- (1/nun)*(nu0*tau0 + (c0/(c0+1))*(mu0-x)^2)
    num <- (sqrt(2*pi)/sqrt(cn))*gamma(nun/2)*(2/(nun*taun))^(nun/2)
    h0x[i,(p1+ptx+j)]<-num/(denom*sqrt(2*pi))  #sqrt(2*pi) is from the data part of likelihood 
  }
}

for(i in 1:n){
  for(j in 1:(p1+ptx)){
    h0x[i,j]<-beta(a0+X[i,j],b0-X[i,j]+1)
  }
}

h0i<-apply(h0x,1,prod)


#========= h0m calculation ========


h0m <- matrix(0, nrow = n, ncol = q)
for ( r in 1:q)
{
  for(i in 1:100){
    betaprior<-as.numeric(rmvn(1,bs1[[r]],inv(ps1[[r]])))
    regvar<-rinvgamma(1,betaa0,betab0)
    h0m[,r]<-h0m[,r]+dnorm(M[,r],(matX%*%betaprior),sqrt(regvar))
  }
}

h0m<-h0m/100
h0m<-apply(h0m,1,prod)

#========= h0y calculation ========

h0y<-numeric(n)
for(i in 1:100){
  betaprior<-as.numeric(rmvn(1,bs,inv(ps)))
  regvar<-rinvgamma(1,betaa0,betab0)
  h0y<-h0y+dnorm(Y,(full_mat%*%betaprior),sqrt(regvar))
}
h0y<-h0y/100



#==========for now just set to some value======== 
alphatheta<-2
alphapsi<-2
alphaomega <- 2
ones <- rep(1,n)
full_dat <- cbind(ones,A,L)

##=================================== UPDATE CLUSTER MEMBERSHIP==================================
##############################################

gibbstot<- 2000
burnin<-50
y <- list()
NDE <- numeric()
NIE <- numeric()
TE <- numeric()
for(gibbsreps in 1:gibbstot){
  
  uniques <- unique(s)
  sortedunique <- uniques[order(uniques[,1],uniques[,2], uniques[,3]), ,drop=FALSE]
  cluster_res <- cluster(X,Y,M,matX, full_cov, s[,1],s[,2],s[,3], sortedunique,beta_y, sig2_y,beta_m, sig_m,xpipars,xmupars,xsigpars,alphatheta, alphapsi,alphaomega,h0y,h0i,h0m,q,nu0,tau0,c0,mu0,a0,b0,betaa0, betab0, bs, bs1)
  s <- cbind(cluster_res$sy,cluster_res$sm,cluster_res$sx)
  beta_y <- cluster_res$beta_y
  sig2_y <- cluster_res$sig2_y
  beta_m <- cluster_res$beta_m
  sig_m <- cluster_res$sig_m
  xpipars <- cluster_res$xpipars
  xmupars <- cluster_res$xmupars
  xsigpars <- cluster_res$xsigpars
  k <- length(unique(s[,1]))
  for(cl in 1:k){
    sig2_y[cl,2]<-newregvar(full_mat[s[,1]==cl,],Y[s[,1]==cl],ps,bs)
    sig <- sig2_y[cl,2]
    beta_y[cl,-1]<-newbet(full_mat[s[,1]==cl,],Y[s[,1]==cl], sig,ps,bs)
  }
  for(row in 1:nrow(sig_m))
  {
    index <- which(s[,1]== sig_m[row,2] & s[,2]==sig_m[row,3])
    indep <- matX[index,]
    dep <- M[index,sig_m[row,1]]
    w <- sig_m[row,1]
    sig_m[row,-c(1:3)] <- newregvar(indep,dep,ps1[[w]],bs1[[w]])
    sig <- sig_m[row,-c(1:3)]
    beta_m[row,-c(1:3)] <- newbet(indep,dep,sig,ps1[[w]],bs1[[w]])
  }
  
  for(row in 1:nrow(xmupars))
  {
    for(j in 1: p2)
    {
      index <- which(s[,1]== xmupars[row,1] & s[,2]== xmupars[row,2] & s[,3]== xmupars[row,3])
      indep <- X[index,p1+ptx+j]
      xsigpars[row,3+j] <- updatevar(indep)
      xmupars[row,3+j] <- updatemean(indep,xsigpars[row,3+j])
    }
  }
  for(row in 1:nrow(xpipars))
  {
    for(j in 1: p1+ptx)
    {
      index <- which(s[,1]== xpipars[row,1] & s[,2]== xpipars[row,2] & s[,3]== xpipars[row,3])
      indep <- X[index,j]
      xpipars[row,3+j] <- rbeta(1,sum(indep)+a0,length(indep)-sum(indep)+b0)
    }
  }
  
  alphapsi <- newalppsi(alphapsi)
  alphaomega <- newalpomega(alphaomega)
  alphatheta <- newalpthet(length(unique(s[,1])),alphatheta)
  
  #================= Computation of Causal Effects===================
  #=========== The first step is to sample L==================
  outcome <- numeric()
  outcome_NIE <- numeric()
  outcome_TE <- numeric()
  for(rep in 1: 1000)
  {
    prob1 <- c(table(s[,1])/(n+ alphatheta), alphatheta/(alphatheta+n))
    f1 <- rmultinomf(prob1)
    if(f1 == length(prob1))
      #========The above means completely new cluster=====  
    {
      l <- numeric()
      for(j in 1:p1)
      {
        temp <- rbeta(1,a0,b0)
        l[j] <- rbinom(1,1,temp)
      }
      for(j in (p1+1):(p1+p2))
      {
        
        temp1 <- rinvchisq(1,nu0, tau0)
        temp2 <- rnorm(1,mu0, sqrt(temp1/c0))
        l[j] <-  rnorm(1,temp2, sqrt(temp1)) 
      }
    }
    #========= Second if condition===========
    if(f1 < length(prob1))
    {
      prob2 <- c(table(s[which(s[,1]==f1),2])/(sum(table(s[which(s[,1]==f1),2]))+ alphapsi), alphapsi/(sum(table(s[which(s[,1]==f1),2]))+ alphapsi))
      f2 <- rmultinomf(prob2)
      #==========Third if within second if=====
      if(f2 ==length(prob2))
      {
        l <- numeric()
        for(j in 1:p1)
        {
          temp <- rbeta(1,a0,b0)
          l[j] <- rbinom(1,1,temp)
        }
        for(j in (p1+1):(p1+p2))
        {
          
          temp1 <- rinvchisq(1,nu0, tau0)
          temp2 <- rnorm(1,mu0, sqrt(temp1/c0))
          l[j] <-  rnorm(1,temp2, sqrt(temp1)) 
        }
      }
      #=========Fourth if within second if====
      if(f2 < length(prob2))
      {
        prob3 <- c(table(s[which(s[,1]==f1 & s[,2]==f2),3])/(sum(table(s[which(s[,1]==f1 & s[,2]==f2),3]))+ alphaomega), alphaomega/(sum(table(s[which(s[,1]==f1 & s[,2]==f2),3]))+ alphaomega))
        f3 <- rmultinomf(prob3)
        #===========Fifth within 4th within 2nd if====
        if(f3==length(prob3))
        {
          l <- numeric()
          for(j in 1:p1)
          {
            temp <- rbeta(1,a0,b0)
            l[j] <- rbinom(1,1,temp)
          }
          for(j in (p1+1):(p1+p2))
          {
            
            temp1 <- rinvchisq(1,nu0, tau0)
            temp2 <- rnorm(1,mu0, sqrt(temp1/c0))
            l[j] <-  rnorm(1,temp2, sqrt(temp1)) 
          }
        }
        if(f3 < length(prob3))
        {
          temp_index <- which(xmupars[,1]==f1 & xmupars[,2]==f2 & xmupars[,3]==f3)
          pi_par <- xpipars[temp_index,-c(1:4)]
          mu_par <- xmupars[temp_index,-c(1:3)]
          sig_par <- xsigpars[temp_index,-c(1:3)]
          x1 <- numeric()
          x2 <- numeric()
          for(f in 1: length(pi_par))
          {
            x1[f] <- rbinom(1,1,pi_par[f]) 
          }
          for(g in 1: length(mu_par))
          {
            x2[g] <- rnorm(1,mu_par[g], sqrt(sig_par[g])) 
          }
          l <- c(x1,x2)
        }
      }
    }
    #print("First part is OK")
    #===========The second step is to sample M given A=0 and A=1 and L====
    a1l <- c(1,l)
    a0l <- c(0,l)
    #=== First we do for A=0=========
    #============= Find the cluster number======
    probs <- numeric()
    #====== Fill in probs for existing clusters============
    g0x<- numeric(p1+ptx+p2)
    
    denom <- (sqrt(2*pi)/sqrt(c0))*gamma(nu0/2)*(2/(nu0*tau0))^(nu0/2)
    cn <- c0 + 1
    nun <- nu0 + 1
    for(j1 in 1:p2){
      x<- a0l[p1+ptx+j1]
      taun <- (1/nun)*(nu0*tau0 + (c0/(c0+1))*(mu0-x)^2)
      num <- (sqrt(2*pi)/sqrt(cn))*gamma(nun/2)*(2/(nun*taun))^(nun/2)
      g0x[(p1+ptx+j1)]<-num/(denom*sqrt(2*pi))  #sqrt(2*pi) is from the data part of likelihood 
    }
    
    for(j2 in 1:(p1+ptx)){
      g0x[j2]<-beta(a0+a0l[j2],b0-a0l[j2]+1)
    }
    
    g0xa1 <- g0x 
    g0xa0 <- g0x
    
    for(j in 1:(ptx)){
      piprior<-rbeta(100,a0,b0)
      g0xa1[j]<-mean(dbinom(1,1,piprior))
      g0xa0[j]<-mean(dbinom(0,1,piprior))
    }
    k_0 <- prod(g0x)
    k_00 <- prod(g0xa0)
    k_01 <- prod(g0xa1)
    
    temp_cl <- xmupars[,1:3]
    temp_cl_2 <- unique(temp_cl[,1:2])
    count <- 1
    for(j in 1: nrow(temp_cl_2))
    {
      n_j_star <- length(which(s[,1]==temp_cl_2[j,1] & s[,2]== temp_cl_2[j,2]))
      x_cl_num <- length(which(temp_cl[,1]==temp_cl_2[j,1] & temp_cl[,2]==temp_cl_2[j,2]))
      for(l in 1: x_cl_num)
      {
        n_l_j_star <- length(which(s[,1]==temp_cl_2[j,1] & s[,2]== temp_cl_2[j,2] & s[,3]==l))
        prodx1 <- 1
        prodx2 <- 1
        index <- which(xpipars[,1]==temp_cl_2[j,1] & xpipars[,2]== temp_cl_2[j,2] & xpipars[,3]==l)
        for(t in 1:(p1+ptx))
        {
          prodx1 <- prodx1* dbinom(a0l[t],1,xpipars[index,3+t])
        }
        for(t in 1:p2)
        {
          prodx2 <- prodx2* dnorm(a0l[p1+ptx+t],xmupars[index,3+t], sqrt(xsigpars[3+t]))
        }
        probs[count] <- (n_j_star/(alphapsi+n)) *(n_l_j_star/(alphaomega+n_j_star))*prodx1*prodx2
        count <- count +1
      }
    }
    
    #====== Fill in probs for existing m but new x clusters============
    count1 <- 1
    for(j in 1: nrow(temp_cl_2))
    {
      n_j_star <- length(which(s[,1]==temp_cl_2[j,1] & s[,2]== temp_cl_2[j,2])) 
      
      probs[nrow(xmupars)+count1] <- (n_j_star/(alphapsi+n))*(alphaomega/(alphaomega+n_j_star))* k_0
      count1 <- count1 + 1  
    }
    
    #====== Fill in probs for completely new cluster============
    probs[nrow(xmupars)+nrow(temp_cl_2)+1] <- (alphapsi/(alphapsi+n))*k_0
    
    pro <- exp(log(probs, base=exp(1)) - max(log(probs,base=exp(1))))
    new_m_cluster <- rmultinomf(pro)
    
    if(new_m_cluster==length(pro))
    {
      temp_mat <- c(1, a0l)
      m <- numeric()
      for(r in 1: q)
      {
        temp_mean <- mvrnorm(1,bs1[[r]],inv(ps1[[r]]))
        temp_sig_2 <- rinvchisq(1,1,1)
        m[r] <- rnorm(1,sum(temp_mean * temp_mat), sqrt(temp_sig_2))
      }
    }
    if(new_m_cluster < length(pro))
    {
      if(new_m_cluster < (nrow(xmupars)+1))
      {
        ff <- xmupars[new_m_cluster,1:2]
        temp_mat <- c(1,a0l)
        m <- numeric()
        m_coef <- beta_m[which(beta_m[,2]==ff[1] & beta_m[,3]==ff[2]),-c(1:3)]
        m_sig <- sig_m[which(sig_m[,2]==ff[1] & sig_m[,3]==ff[2]),-c(1:3)]
        for(r in 1: q)
        {
          m[r] <- rnorm(1,sum(m_coef[r,] * temp_mat), sqrt(m_sig[r]))
        }
      }else{
        f <- new_m_cluster - nrow(xmupars)
        ff <- unique(xmupars[,1:2])[f,]
        temp_mat <- c(1,a0l)
        m <- numeric()
        m_coef <- beta_m[which(beta_m[,2]==ff[1] & beta_m[,3]==ff[2]),-c(1:3)]
        m_sig <- sig_m[which(sig_m[,2]==ff[1] & sig_m[,3]==ff[2]),-c(1:3)]
        for(r in 1: q)
        {
          m[r] <- rnorm(1,sum(m_coef[r,] * temp_mat), sqrt(m_sig[r]))
        }
      }
    }
    
    #=== Now we do for A=1=========
    #============= Find the cluster number======
    probs <- numeric()
    #====== Fill in probs for existing clusters============
    g0x<- numeric(p1+ptx+p2)
    
    denom <- (sqrt(2*pi)/sqrt(c0))*gamma(nu0/2)*(2/(nu0*tau0))^(nu0/2)
    cn <- c0 + 1
    nun <- nu0 + 1
    for(j1 in 1:p2){
      x<- a1l[p1+ptx+j1]
      taun <- (1/nun)*(nu0*tau0 + (c0/(c0+1))*(mu0-x)^2)
      num <- (sqrt(2*pi)/sqrt(cn))*gamma(nun/2)*(2/(nun*taun))^(nun/2)
      g0x[(p1+ptx+j1)]<-num/(denom*sqrt(2*pi))  #sqrt(2*pi) is from the data part of likelihood 
    }
    
    for(j2 in 1:(p1+ptx)){
      g0x[j2]<-beta(a0+a1l[j2],b0-a1l[j2]+1)
    }
    
    g0xa1 <- g0x 
    g0xa0 <- g0x
    
    for(j in 1:(ptx)){
      piprior<-rbeta(100,a0,b0)
      g0xa1[j]<-mean(dbinom(1,1,piprior))
      g0xa0[j]<-mean(dbinom(0,1,piprior))
    }
    k_1 <- prod(g0x)
    k_10 <- prod(g0xa0)
    k_11 <- prod(g0xa1)
    
    temp_cl <- xmupars[,1:3]
    temp_cl_2 <- unique(temp_cl[,1:2])
    count <- 1
    for(j in 1: nrow(temp_cl_2))
    {
      n_j_star <- length(which(s[,1]==temp_cl_2[j,1] & s[,2]== temp_cl_2[j,2]))
      x_cl_num <- length(which(temp_cl[,1]==temp_cl_2[j,1] & temp_cl[,2]==temp_cl_2[j,2]))
      for(l in 1: x_cl_num)
      {
        n_l_j_star <- length(which(s[,1]==temp_cl_2[j,1] & s[,2]== temp_cl_2[j,2] & s[,3]==l))
        prodx1 <- 1
        prodx2 <- 1
        index <- which(xpipars[,1]==temp_cl_2[j,1] & xpipars[,2]== temp_cl_2[j,2] & xpipars[,3]==l)
        for(t in 1:(p1+ptx))
        {
          prodx1 <- prodx1* dbinom(a1l[t],1,xpipars[index,3+t])
        }
        for(t in 1:p2)
        {
          prodx2 <- prodx2* dnorm(a1l[p1+ptx+t],xmupars[index,3+t], sqrt(xsigpars[3+t]))
        }
        probs[count] <- (n_j_star/(alphapsi+n)) *(n_l_j_star/(alphaomega+n_j_star))*prodx1*prodx2
        count <- count +1
      }
    }
    
    #====== Fill in probs for existing m but new x clusters============
    count1 <- 1
    for(j in 1: nrow(temp_cl_2))
    {
      n_j_star <- length(which(s[,1]==temp_cl_2[j,1] & s[,2]== temp_cl_2[j,2])) 
      
      probs[nrow(xmupars)+count1] <- (n_j_star/(alphapsi+n))*(alphaomega/(alphaomega+n_j_star))* k_1
      count1 <- count1 + 1  
    }
    
    #====== Fill in probs for completely new cluster============
    probs[nrow(xmupars)+nrow(temp_cl_2)+1] <- (alphapsi/(alphapsi+n))*k_1
    
    pro <- exp(log(probs, base=exp(1)) - max(log(probs,base=exp(1))))
    new_m_cluster <- rmultinomf(pro)
    
    if(new_m_cluster==length(pro))
    {
      temp_mat <- c(1, a1l)
      m1 <- numeric()
      for(r in 1: q)
      {
        temp_mean <- mvrnorm(1,bs1[[r]],inv(ps1[[r]]))
        temp_sig_2 <- rinvchisq(1,1,1)
        m1[r] <- rnorm(1,sum(temp_mean * temp_mat), sqrt(temp_sig_2))
      }
    }
    if(new_m_cluster < length(pro))
    {
      if(new_m_cluster < (nrow(xmupars)+1))
      {
        ff <- xmupars[new_m_cluster,1:2]
        temp_mat <- c(1,a1l)
        m1 <- numeric()
        m_coef <- beta_m[which(beta_m[,2]==ff[1] & beta_m[,3]==ff[2]),-c(1:3)]
        m_sig <- sig_m[which(sig_m[,2]==ff[1] & sig_m[,3]==ff[2]),-c(1:3)]
        for(r in 1: q)
        {
          m1[r] <- rnorm(1,sum(m_coef[r,] * temp_mat), sqrt(m_sig[r]))
        }
      }else{
        f <- new_m_cluster - nrow(xmupars)
        ff <- unique(xmupars[,1:2])[f,]
        temp_mat <- c(1,a1l)
        m1 <- numeric()
        m_coef <- beta_m[which(beta_m[,2]==ff[1] & beta_m[,3]==ff[2]),-c(1:3)]
        m_sig <- sig_m[which(sig_m[,2]==ff[1] & sig_m[,3]==ff[2]),-c(1:3)]
        for(r in 1: q)
        {
          m1[r] <- rnorm(1,sum(m_coef[r,] * temp_mat), sqrt(m_sig[r]))
        }
      }
    }
    #=========================
    
    #print("second part is OK")
    #============= The last step is to compute the conditional expectation of Y=========
    k0m <- rep(0,q)
    k1m <- rep(0,q)
    k0m1 <- rep(0,q)
    k1m1 <- rep(0,q)
    
    for ( r in 1:q)
    {
      for(i in 1:100){
        betaprior <- as.numeric(rmvn(1,bs1[[r]],inv(ps1[[r]])))
        regvar <- rinvgamma(1,betaa0,betab0)
        k0m[r] <- k0m[r]+dnorm(m[r],sum(c(1,a0l)*betaprior),sqrt(regvar))
        k1m[r] <- k1m[r]+dnorm(m[r],sum(c(1,a1l)*betaprior),sqrt(regvar))
        k0m1[r] <- k0m1[r]+dnorm(m1[r],sum(c(1,a0l)*betaprior),sqrt(regvar))
        k1m1[r] <- k1m1[r]+dnorm(m1[r],sum(c(1,a1l)*betaprior),sqrt(regvar))
      }
    }
    
    k0m<-k0m/100
    k1m<-k1m/100
    k0m<- prod(k0m)
    k1m<- prod(k1m)
    
    k0m1<-k0m1/100
    k1m1<-k1m1/100
    k0m1<- prod(k0m1)
    k1m1<- prod(k1m1)
    #Part_1============================
    sum1 <- 0
    sum1_new <- 0
    sum1_second <- 0
    sum11 <- 0
    sum11_new <- 0
    sum11_second <- 0
    sum20 <- 0
    sum21 <- 0
    sum21_new <- 0
    sum30 <- 0
    sum31 <- 0
    sord <- s[order(s[,1],s[,2], s[,3]), ,drop=FALSE]
    for(j in 1: length(unique(sord[,1])))
    {
      
      sub_s <- sord[which(sord[,1]==j),]
      if(is.null(nrow(sub_s))==T)
      {
        k_j <- 1
      }else{
        k_j <- length(unique(sub_s[,2]))
      }
      for(l in 1: k_j)
      {
        sub_sub_s <- sord[which(sord[,1]==j & sord[,2]==l),]
        if(is.null(nrow(sub_sub_s))==T)
        {
          k_jl <- 1
        }else{
          k_jl <- length(unique(sub_sub_s[,3]))
        }
        for(u in 1: k_jl)
        {
          part31 <- length(which(sord[,1]==j & sord[,2]==l & sord[,3]==u))/(alphaomega +length(which(sord[,1]==j & sord[,2]==l)))
          prodx10 <-1
          prodx11 <-1
          prodx2 <- 1
          for(t1 in 1: p1+ptx)
          {
            index <- which(xpipars[,1]==j & xpipars[,2]==l & xpipars[,3]==u)
            prodx10 <- prodx10* dbinom(a0l[t1],1,xpipars[index,3+t1])
            prodx11 <- prodx11* dbinom(a1l[t1],1,xpipars[index,3+t1])
          }
          for(t2 in 1:p2)
          {
            prodx2 <- prodx2* dnorm(a0l[p1+ptx+t],xmupars[index,3+t2], sqrt(xsigpars[3+t2])) 
          }
          
          sum30 <- sum30 + part31 * prodx10 * prodx2
          sum31 <- sum31 + part31 * prodx11 * prodx2
        }
        part21 <- length(which(sord[,1]==j & sord[,2]==l))/(alphapsi+length(which(sord[,1]==j)))
        m_coef <- beta_m[which(beta_m[,2]==j & beta_m[,3]==l),-c(1:3)]
        m_sig <- sig_m[which(sig_m[,2]==j & sig_m[,3]==l),-c(1:3)]
        likem0 <- 1
        likem1 <- 1
        likem1_new <- 1
        for(r in 1: q)
        {
          likem0 <- likem0 * dnorm(m[r], sum(c(1,a0l)*m_coef[r,]),sqrt(m_sig[r]))
          likem1 <- likem1 * dnorm(m[r], sum(c(1,a1l)*m_coef[r,]),sqrt(m_sig[r]))
          likem1_new <- likem1_new * dnorm(m1[r], sum(c(1,a1l)*m_coef[r,]),sqrt(m_sig[r]))
        }
        sum20 <- sum20 + part21 * likem0 * sum30
        sum21 <- sum21 + part21 * likem1 * sum30
        sum21_new <- sum21_new + part21 * likem1_new * sum30
      }
      E <- sum(c(1,m,a1l) * beta_y[j,-1])
      E_new <- sum(c(1,m1,a1l) * beta_y[j,-1])
      E_second <- sum(c(1,m,a0l) * beta_y[j,-1])
      sum1 <- sum1 + length(which(s[,1]==j)) /(alphatheta+n) * sum21 * E
      sum1_new <- sum1_new + length(which(s[,1]==j)) /(alphatheta+n) * sum21_new * E_new
      sum1_second <- sum1_second + length(which(s[,1]==j)) /(alphatheta+n) * sum20 * E_second
      sum11 <- sum11 + length(which(s[,1]==j)) /(alphatheta+n) * sum21
      sum11_new <- sum11_new + length(which(s[,1]==j)) /(alphatheta+n) * sum21_new
      sum11_second <- sum11_second + length(which(s[,1]==j)) /(alphatheta+n) * sum20 
    }
    part1_num <- as.numeric(sum1)
    part1_num_new <- as.numeric(sum1_new)
    part1_num_second <- as.numeric(sum1_second)
    part1_denom <- sum11
    part1_denom_new <- sum11_new
    part1_denom_second <- sum11_second
    #print("part1 is ok")
    #Part_2=============================
    sum1 <- 0
    sum1_new <- 0
    sum1_second <- 0
    sum11 <- 0
    sum11_new <- 0
    sum11_second <- 0
    sum20 <- 0
    sum21 <- 0
    sum21_new <- 0
    
    
    for(j in 1: length(unique(sord[,1])))
    {
      
      sub_s <- sord[which(sord[,1]==j),]
      if(is.null(nrow(sub_s))==T)
      {
        k_j <- 1
      }else{
        k_j <- length(unique(sub_s[,2]))
      }
      for(l in 1: k_j)
      {
        
        part21 <- length(which(sord[,1]==j & sord[,2]==l))/(alphapsi+length(which(sord[,1]==j)))
        part22 <- alphaomega/(alphaomega+length(which(sord[,1]==j & sord[,2]==l)))
        m_coef <- beta_m[which(beta_m[,2]==j & beta_m[,3]==l),-c(1:3)]
        m_sig <- sig_m[which(sig_m[,2]==j & sig_m[,3]==l),-c(1:3)]
        likem0 <- 1
        likem1 <- 1
        likem1_new <- 1
        for(r in 1: q)
        {
          likem0 <- likem0 * dnorm(m[r], sum(c(1,a0l)*m_coef[r,]),sqrt(m_sig[r]))
          likem1 <- likem1 * dnorm(m[r], sum(c(1,a1l)*m_coef[r,]),sqrt(m_sig[r]))
          likem1_new <- likem1_new * dnorm(m1[r], sum(c(1,a1l)*m_coef[r,]),sqrt(m_sig[r]))
        }
        sum20 <- sum20 + part21 * part22 * likem0
        sum21 <- sum21 + part21 * part22 * likem1
        sum21_new <- sum21_new + part21 * part22 * likem1_new
        
      }
      E <- sum(c(1,m,a1l) * beta_y[j,-1])
      E_new <- sum(c(1,m1,a1l) * beta_y[j,-1])
      E_second <- sum(c(1,m,a0l) * beta_y[j,-1])
      sum1 <- sum1 + length(which(sord[,1]==j)) /(alphatheta+n) * sum21 * E * k_01
      sum1_new <- sum1_new + length(which(sord[,1]==j)) /(alphatheta+n) * sum21_new * E_new * k_11
      sum1_second <- sum1_second + length(which(sord[,1]==j)) /(alphatheta+n) * sum20 * E_second * k_00
      sum11 <- sum11 + length(which(sord[,1]==j)) /(alphatheta+n) * sum21 * k_01
      sum11_new <- sum11_new + length(which(sord[,1]==j)) /(alphatheta+n) * sum21_new * k_11
      sum11_second <- sum11_second + length(which(sord[,1]==j)) /(alphatheta+n) * sum20 * k_00
    }
    part2_num <- as.numeric(sum1)
    part2_num_new <- as.numeric(sum1_new)
    part2_num_second <- as.numeric(sum1_second)
    part2_denom <- sum11
    part2_denom_new <- sum11_new
    part2_denom_second <- sum11_second
    #print("part2 is ok")
    #part3===============================================
    sum1 <- 0
    sum1_new <- 0
    sum1_second <- 0
    sum11 <- 0
    sum11_new <- 0
    sum11_second <- 0
    for(j in 1: length(unique(sord[,1])))
    {
      E <- sum(c(1,m,a1l) * beta_y[j,-1])
      E_new <- sum(c(1,m1,a1l) * beta_y[j,-1])
      E_second <- sum(c(1,m,a0l) * beta_y[j,-1])
      sum1 <- sum1 + (length(which(sord[,1]==j)) /(alphatheta+n)) *(alphapsi/(alphapsi+length(which(sord[,1]==j))))  * E * k_01 * k1m
      sum1_new <- sum1_new + (length(which(sord[,1]==j)) /(alphatheta+n)) *(alphapsi/(alphapsi+length(which(sord[,1]==j))))  * E_new * k_11 * k1m1
      sum1_second <- sum1_second + (length(which(sord[,1]==j)) /(alphatheta+n)) *(alphapsi/(alphapsi+length(which(sord[,1]==j))))  * E_second * k_00 * k0m
      sum11 <- sum11 +(length(which(sord[,1]==j)) /(alphatheta+n)) * (alphapsi/(alphapsi+length(which(sord[,1]==j))))  * k_01 * k1m
      sum11_new <- sum11_new +(length(which(sord[,1]==j)) /(alphatheta+n)) * (alphapsi/(alphapsi+length(which(sord[,1]==j))))  * k_11 * k1m1
      sum11_second <- sum11_second +(length(which(sord[,1]==j)) /(alphatheta+n)) * (alphapsi/(alphapsi+length(which(sord[,1]==j))))  * k_00 * k0m
    }
    part3_num <- as.numeric(sum1)
    part3_num_new <- as.numeric(sum1_new)
    part3_num_second <- as.numeric(sum1_second)
    part3_denom <- sum11
    part3_denom_new <- sum11_new
    part3_denom_second <- sum11_second
    #print("part3 is ok")
    #part4===============================================
    E0y<- 0
    E1y <- 0
    E1y_new <- 0
    for(i in 1:100){
      betaprior<-as.numeric(rmvn(1,bs,inv(ps)))
      regvar<-rinvgamma(1,betaa0,betab0)
      E0y<-E0y+ sum(c(1,m,a0l) * betaprior)
      E1y <- E1y+ sum(c(1,m,a1l) * betaprior)
      E1y_new <- E1y_new + sum(c(1,m1,a1l) * betaprior)
    }
    E0y<-E0y/100
    E1y<-E1y/100
    E1y_new<-E1y_new/100
    part4_num <- (alphatheta/(alphatheta+n))* k1m * k_01 * E1y
    part4_num_new <- (alphatheta/(alphatheta+n))* k1m1 * k_11 * E1y_new
    part4_num_second <- (alphatheta/(alphatheta+n))* k0m * k_00 * E0y
    part4_denom <- (alphatheta/(alphatheta+n))* k1m * k_01
    part4_denom_new <- (alphatheta/(alphatheta+n))* k1m1 * k_11
    part4_denom_second <- (alphatheta/(alphatheta+n))* k0m * k_00
    #print("part 4 is ok")
    
    Pot_outcome_1 <- (part1_num+part2_num+part3_num+part4_num)/(part1_denom+part2_denom+part3_denom+part4_denom)
    Pot_outcome_2 <- (part1_num_second+part2_num_second+part3_num_second+part4_num_second)/(part1_denom_second+part2_denom_second+part3_denom_second+part4_denom_second)
    Pot_outcome_3 <- (part1_num_new+part2_num_new+part3_num_new+part4_num_new)/(part1_denom_new+part2_denom_new+part3_denom_new+part4_denom_new)
    
    diff <- Pot_outcome_1 - Pot_outcome_2
    diff_NIE <- Pot_outcome_3 - Pot_outcome_1
    diff_TE <- Pot_outcome_3 - Pot_outcome_2
    outcome <- c(outcome, diff)
    outcome_NIE <- c(outcome_NIE,diff_NIE)
    outcome_TE <- c(outcome_TE,diff_TE)
  }
  NDE_iter <- mean(outcome, na.rm = T)
  NIE_iter <- mean(outcome_NIE, na.rm = T)
  TE_iter <- mean(outcome_TE, na.rm = T)
  NDE <- c(NDE,NDE_iter)
  NIE <- c(NIE,NIE_iter)
  TE <- c(TE, TE_iter)
  print(c(gibbsreps))
}
#=======Here we gather the final output====

NIE <- NIE[-c(1:1000)]

print(mean(NIE, na.rm = T))
print(quantile(NIE, probs = c(0.025, 0.975), na.rm = T))
