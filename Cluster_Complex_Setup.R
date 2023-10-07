library(geometry)

#=======insertvec function=========
insertvec <- function(v,p,e)
{
  if(p==1)
  {
    result <- c(e,v)
  }
  if(1<p & p<(length(v)+1))
  {
    result <- c(v[c(1:(p-1))],e,v[c(p:length(v))])
  }
  if(p > length(v))
  {
    result <- c(v,e)
  }
  result
}
#=======insertmat function =========
insertmat <- function(m,p,r)
{
  if (p==1)
  {
    result <- rbind(r,m)
  }
  if(1 < p & p < (length(m[,1])+1))
  {
    result <- rbind(m[c(1:(p-1)),],r,m[c(p:length(m[,1])),])
  }
  if(p>length(m[,1]))
  {
    result <- rbind(m,r)  
  }
  
  result
}

#======= rmultinomf function===================
rmultinomf <- function(p)
{
  csp <- cumsum(p)/sum(p)
  rnd <- runif(1)
  res <- 0
  for ( i in 1: length(p))
  {
    if(rnd>csp[i])
    {
      res <- res+1
    }
  }
  return(res+1)
}

#(1)============updtvar function======================

updtvar<-function(tempx,nu0,tau0,c0,mu0){
  newdf<-nu0+length(tempx)
  if(length(tempx)==1){varx<-0}
  if(length(tempx)>1){varx<-var(tempx)}
  numer<-nu0*tau0+(length(tempx)-1)*varx+(c0*length(tempx)/(c0+length(tempx)))*(mean(tempx)-mu0)^2
  newval<- rinvchisq(1,newdf,numer/newdf)
  return(newval)
}

#==(2)===========updtmean function=======================

updtmean<-function(tempx,curtau,c0,mu0){
  newvar<-1/(c0/curtau+length(tempx)/curtau)
  newmean<-(mu0*c0/curtau+mean(tempx)*length(tempx)/curtau)*newvar
  newval<-rnorm(1,newmean,sqrt(newvar))
  return(newval)
}

#===============nwrgvar function===============
nwrgvar<-function(tempx,tempy, betaa0, betab0,beta0){
  an<-betaa0+length(tempy)/2
  ifelse((ncol(as.matrix(tempx))>1),(xtx<-t(tempx)%*%tempx),(xtx<-as.matrix(tempx)%*%t(as.matrix(tempx))))
  prec0 <- 1* diag(ncol(xtx))
  #beta0 <- rep(0,ncol(xtx))   
  newprec<-xtx+prec0
  if(ncol(as.matrix(tempx))>1){betan<-(ginv(xtx+prec0))%*%(prec0%*%beta0+t(tempx)%*%tempy)}
  if(ncol(as.matrix(tempx))<=1){betan<-(ginv(xtx+prec0))%*%(prec0%*%beta0+(tempx*tempy))}
  bn<-betab0+.5*(sum(tempy*tempy)+t(beta0)%*%prec0%*%beta0-t(betan)%*%newprec%*%betan)
  return(rinvgamma(1,an,bn))
}

#==========================nwbet function=====================================
nwbet<-function(tempx,tempy,sig, beta0){
  ifelse((ncol(as.matrix(tempx))>1),(xtx<-t(tempx)%*%tempx),(xtx<-as.matrix(tempx)%*%t(as.matrix(tempx))))
  prec0 <- 1* diag(ncol(xtx))
  #beta0 <- rep(0,ncol(xtx))   
  newprec<-xtx+prec0
  if(ncol(as.matrix(tempx))>1){betan<-(ginv(xtx+prec0))%*%(prec0%*%beta0+t(tempx)%*%tempy)}
  if(ncol(as.matrix(tempx))<=1){betan<-(ginv(xtx+prec0))%*%(prec0%*%beta0+(tempx*tempy))}
  #print(newprec)
  return(rmvn(1,betan,sig*inv(newprec)))
}


#===================== Cluster function===============
cluster <- function(X,y,M,matx, full_cov, sy,sm,sx,uniqueS, beta_y, sig2_y,beta_m, sig_m,xPiPars, xMuPars, xSigPars,alpha, mu,psai, h0y, h0i, h0m,q,nu0,tau0,c0,mu0,a0,b0,betaa0, betab0, bs, bs1)
{
  n <- nrow(X)
  p <- ncol(X)
  count_mat <- numeric()
  check <- list()
  s_check <- list()
  
  #=====loop through every person=======    
  #===== Check whether he is alone in his cluster and act accordingly=====
  multr <- numeric()
  #for (i in 1:n)
  for(i in 1:n)
  {
    
    # First get the "drop_index" , which will be used later
    
    get_index <- as.numeric(which(uniqueS[,1]==sy[i] & uniqueS[,2]==sm[i] & uniqueS[,3]==sx[i]))
    get_index1 <- as.numeric(which(uniqueS[,1]==sy[i] & uniqueS[,2]==sm[i]))
    get_index2 <- as.numeric(which(uniqueS[,1]==sy[i]))
    
    drop_index <- NULL
    
    if(length(which(sy==sy[i] & sm == sm[i] & sx== sx[i])) == 1)
    {
      drop_index <- get_index
      
      if(length(which(sy==sy[i] & sm == sm[i]))==1)
      {
        drop_index <- get_index1
        
        if(length(which(sy==sy[i]))==1)
          
        {
          drop_index <- get_index2
        }
      }
    }
    # You have the drop_index now. You can now go ahead
    #Now we will use three "if" conditions in nested form.
    
    # ===== The first "if" condition================
    if(length(which(sy==sy[i] & sm == sm[i] & sx== sx[i])) == 1)
    {
      
      #====This is the second "if" condition============
      
      if(length(which(sy==sy[i] & sm == sm[i]))==1)
      {
        #--------------Delete beta, sigma, mu, pi---------
        #---------Note that, we need to treat both vector and matrix case------
        
        if(is.matrix(beta_m)==T)
        {
          beta_m <- beta_m[-which(beta_m[,2]==sy[i] & beta_m[,3]==sm[i]),]
        }else{
          beta_m <- beta_m[-which(beta_m[,2]==sy[i] & beta_m[,3]==sm[i])]
          beta_m <- t(as.matrix(beta_m))
        }
        if(is.matrix(sig_m)==T)
        {
          sig_m <- sig_m[-which(sig_m[,2]==sy[i]& sig_m[,3]==sm[i]),]
        }else{
          sig_m <- sig_m[-which(sig_m[,2]==sy[i]& sig_m[,3]==sm[i])]
          sig_m <- t(as.matrix(sig_m))
        }
        if(is.matrix(xPiPars)==T)
        {
          xPiPars <- xPiPars[-which(xPiPars[,1]==sy[i] & xPiPars[,2]==sm[i]),]
        }else{
          xPiPars <- xPiPars[-which(xPiPars[,1]==sy[i] & xPiPars[,2]==sm[i])] 
          xPiPars <- t(as.matrix(xPiPars))
        }
        if(is.matrix(xMuPars)==T)
        {
          xMuPars <- xMuPars[-which(xMuPars[,1]==sy[i] & xMuPars[,2]==sm[i]),]
        }else{
          xMuPars <- xMuPars[-which(xMuPars[,1]==sy[i] & xMuPars[,2]==sm[i])]  
          xMuPars <- t(as.matrix(xMuPars))
        }
        if(is.matrix(xSigPars)==T)
        {
          xSigPars <- xSigPars[-which(xSigPars[,1]==sy[i] & xSigPars[,2]==sm[i]),]
        }else{
          xSigPars <- xSigPars[-which(xSigPars[,1]==sy[i] & xSigPars[,2]==sm[i])] 
          xSigPars <- t(as.matrix(xSigPars))
        }
        
        #-----This makes sure that everything is in matrix form before going into next "if" condition---
        if(is.matrix(beta_m)==F)
        {
          beta_m <- t(as.matrix(beta_m))
        }
        if(is.matrix(sig_m)==F)
        {
          sig_m <- t(as.matrix(sig_m))
        }
        if(is.matrix(xPiPars)==F)
        {
          xPiPars <- t(as.matrix(xPiPars))
        }
        if(is.matrix(xMuPars)==F)
        {
          xMuPars <- t(as.matrix(xMuPars))
        }
        if(is.matrix(xSigPars)==F)
        {
          xSigPars <- t(as.matrix(xSigPars))
        }
        #======This is the third "if" condition================ 
        
        if(length(which(sy==sy[i]))==1)
        { 
          #--------------Delete beta, sigma, mu, pi---------
          if(is.matrix(beta_y)==T)
          {
            beta_y <- beta_y[-sy[i],]
          }else{
            beta_y <- beta_y[-sy[i]]
            beta_y <- t(as.matrix(beta_y))
          }
          if(is.matrix(sig2_y)==T)
          {
            sig2_y <- sig2_y[-sy[i],]
          }else{
            sig2_y <- sig2_y[-sy[i]]
            sig2_y <- t(as.matrix(sig2_y))
          }
          if(length(which(beta_m[,2]==sy[i])) > 0)
          {
            beta_m <- beta_m[-which(beta_m[,2]==sy[i]),]
          }
          if(length(which(sig_m[,2]==sy[i])) > 0)
          {
            sig_m <- sig_m[-which(sig_m[,2]==sy[i]),]
          }
          if(length(which(xPiPars[,1]==sy[i])) > 0)
          {
            xPiPars <- xPiPars[-which(xPiPars[,1]==sy[i]),]
          }
          if(length(which(xMuPars[,1]==sy[i])) > 0)
          {
            xMuPars <- xMuPars[-which(xMuPars[,1]==sy[i]),]
          }
          if(length(which(xSigPars[,1]==sy[i])) > 0)
          {
            xSigPars <- xSigPars[-which(xSigPars[,1]==sy[i]),]
          }
          #--------- Making sure everything is matrix before the next "if" condition-----
          if(is.matrix(beta_y)==F)
          {
            beta_y <- t(as.matrix(beta_y))
          }
          if(is.matrix(sig2_y)==F)
          {
            sig2_y <- t(as.matrix(sig2_y))
          }
          if(is.matrix(beta_m)==F)
          {
            beta_m <- t(as.matrix(beta_m))
          }
          if(is.matrix(sig_m)==F)
          {
            sig_m <- t(as.matrix(sig_m))
          }
          if(is.matrix(xPiPars)==F)
          {
            xPiPars <- t(as.matrix(xPiPars))
          }
          if(is.matrix(xMuPars)==F)
          {
            xMuPars <- t(as.matrix(xMuPars))
          }
          if(is.matrix(xSigPars)==F)
          {
            xSigPars <- t(as.matrix(xSigPars))
          }
        }
        # end of third "if" condition=====
      }
      #end of the second "if" condition
      
      #==If none of the 2nd and 3rd "if" holds, we jump to the following step directly===
      #Please note that these "if conditions" need to be checked !
      
      if(length(which(xPiPars[,1]==sy[i] & xPiPars[,2]==sm[i] & xPiPars[,3]==sx[i])) > 0)
      {
        xPiPars <- xPiPars[-which(xPiPars[,1]==sy[i] & xPiPars[,2]==sm[i] & xPiPars[,3]==sx[i]),]
      }
      if(length(which(xMuPars[,1]==sy[i] & xMuPars[,2]==sm[i] & xMuPars[,3]==sx[i])) > 0)
      {
        xMuPars <- xMuPars[-which(xMuPars[,1]==sy[i] & xMuPars[,2]==sm[i] & xMuPars[,3]==sx[i]),]
      }
      if (length(which(xSigPars[,1]==sy[i] & xSigPars[,2]==sm[i] & xSigPars[,3]==sx[i])) > 0)
      {
        xSigPars <- xSigPars[-which(xSigPars[,1]==sy[i] & xSigPars[,2]==sm[i] & xSigPars[,3]==sx[i]),]
      }
      
      
      #======== Re-label x clusters==============
      #Note : Remember we are just Re-labeling in this step. We will delete later
      
      for (j in 1:length(sx))
      {
        if(sy[j]==sy[i] & sm[j]==sm[i] & sx[j] > sx[i])
        {
          sx[j] = sx[j] -1
        }
      }
      for ( j in 1: nrow(uniqueS))
      {
        if(uniqueS[j,1]==sy[i]& uniqueS[j,2]==sm[i] & uniqueS[j,3] > sx[i])
        {
          uniqueS[j,3] = uniqueS[j,3] - 1
        }
      }
      
      #3--------------Re-Label M cluster---------------------------
      
      if(length(which(sy==sy[i] & sm == sm[i]))==1)
      {
        for (j in 1:length(sm))
        {
          if(sy[j]==sy[i] & sm[j] > sm[i])
          {
            sm[j] = sm[j] -1
          }
        }
        for ( j in 1: nrow(uniqueS))
        {
          if(uniqueS[j,1]==sy[i]& uniqueS[j,2]>sm[i])
          {
            uniqueS[j,2] = uniqueS[j,2] - 1
          }
        }
      }
      
      #3--------------Re-Label Y cluster---------------------------
      if(length(which(sy==sy[i]))==1)
      {
        for (j in 1:length(sy))
        {
          if(sy[j]>sy[i])
          {
            sy[j] = sy[j] -1
          }
        }
        for ( j in 1: nrow(uniqueS))
        {
          if(uniqueS[j,1]>sy[i])
          {
            uniqueS[j,1] = uniqueS[j,1] - 1
          }
        }
      }
      #=============Now Get rid of the rows of uniqueS===========
      if(length(drop_index)>0)
      {
        if(nrow(uniqueS)==2)
        {
          uniqueS <- t(as.matrix(uniqueS[-drop_index,]))
        }else{
          uniqueS <- uniqueS[-drop_index,]   
        } 
      }
    }
    #end of first if condition
    
    #======Please note that the followings are needed before we proceed !======
    #------ This is just to attach the cluster numbers to the parameters------
    #print(i)
    if(nrow(uniqueS)==1)
    {
      xMuPars <- cbind(uniqueS,t(as.matrix(xMuPars[-c(1:3)])))
      xSigPars <- cbind(uniqueS,t(as.matrix(xSigPars[-c(1:3)])))
      xPiPars <- cbind(uniqueS,t(as.matrix(xPiPars[-c(1:3)])))
    }else{
      xMuPars <- cbind(uniqueS,xMuPars[,-c(1:3)])
      xSigPars <- cbind(uniqueS,xSigPars[,-c(1:3)])
      xPiPars <- cbind(uniqueS,xPiPars[,-c(1:3)])
    }
    
    if(length(unique(uniqueS[,1]))==1)
    {
      beta_y <- cbind(unique(uniqueS[,1]),t(as.matrix(beta_y[-1])))
      sig2_y <- cbind(unique(uniqueS[,1]),t(as.matrix(sig2_y[-1])))
    }else{
      beta_y <- cbind(unique(uniqueS[,1]),beta_y[,-1])
      sig2_y <- cbind(unique(uniqueS[,1]), sig2_y[,-1])
    }
    
    if(nrow(uniqueS)==1)
    {
      select <- uniqueS[,c(1,2)]
      select1 <- cbind(rep(select[1],q),rep(select[2],q))
    }else{
      select <- unique(uniqueS[,c(1,2)])
      select1 <- cbind(rep(select[,1],q),rep(select[,2],q))
    }
    #print(beta_m)
    #print(select1)
    beta_m[,c(2,3)] <- select1
    sig_m[,c(2,3)] <- select1
    
    #===============Now delete the 'i'-th row of sy, sm and sx=====
    sy <- sy[-i]
    sm <- sm[-i]
    sx <- sx[-i]
    
    #==== Re-calculate number of unique clusters==========
    
    # Note that numy, numm and ntot can be obtained from uniqueS also
    numy <- nrow(beta_y)
    numm <- length(which(beta_m[,1]==1))
    ntot <- nrow(xMuPars)
    totalposs <- numy + numm + ntot + 1
    probs <- numeric(totalposs)
    
    #====== Fill in probs for existing clusters============
    count <- 1
    for ( j in 1: numy)
    {
      
      sub11 <- uniqueS[which(uniqueS[,1]==j),]
      if(nrow(as.matrix(t(sub11)))==1)
      {
        nummj <- length(unique(sub11[2]))
      }else{
        nummj <- length(unique(sub11[,2]))
      }
      nj_woi <- length(which(sy==j))
      likeregy <- dnorm(y[i], dot(full_mat[i,], beta_y[j,-c(1)]),sqrt(sig2_y[j,2]))
      for ( k in 1: nummj)
      {
        sub1 <- uniqueS[which(uniqueS[,1]==j & uniqueS[,2]==k),]
        if(nrow(as.matrix(t(sub1)))==1)
        {
          numx_kj <- length(unique(sub1[3]))
        }else{
          numx_kj <- length(unique(sub1[,3]))
        }
        nkj_woi <- length(which(sy==j & sm ==k))
        likeregm <- 1
        for ( r in 1:q)
        {
          temp <- dnorm(M[i,r], dot(matX[i,], beta_m[which(beta_m[,1]==r & beta_m[,2]==j & beta_m[,3]==k),-c(1:3)]),sqrt(sig_m[which(sig_m[,1]==r & sig_m[,2]==j & sig_m[,3]==k),4]))
          likeregm <- likeregm*temp
        }
        
        for ( l in 1: numx_kj)
        {
          nlkj_woi <- length(which(sy==j & sm ==k & sx==l))
          prodx1 <- 1
          prodx2 <- 1
          for ( t in 1:(ptx+p1))
          {
            prodx1 <- prodx1 * dbinom(X[i,t],1,xPiPars[count,(3+t)])
          }
          for ( t in 1:p2)
          {
            prodx2 <- prodx2 * dnorm(X[i,ptx+p1+t],xMuPars[count,(3+t)],sqrt(xSigPars[count,(3+t)]))
          }
          probs[count] <- (nj_woi/(alpha+n))* (nkj_woi/(mu+nj_woi)) * (nlkj_woi/(psai+nkj_woi))*likeregy*likeregm*prodx1*prodx2
          count <- count + 1
        }
      }
    }
    #====== Fill in probs for old Y and M clusters, but new x cluster============
    count1 <- 1
    for ( j in 1: numy)
    {
      sub11 <- uniqueS[which(uniqueS[,1]==j),]
      if(nrow(as.matrix(t(sub11)))==1)
      {
        nummj <- length(unique(sub11[2]))
      }else{
        nummj <- length(unique(sub11[,2]))
      }
      nj_woi <- length(which(sy==j))
      likeregy <- dnorm(y[i], dot(full_mat[i,], beta_y[j,-c(1)]),sqrt(sig2_y[j,2]))
      for ( k in 1: nummj)
      {
        #numx_kj <- length(unique((uniqueS[which(uniqueS[,1]==j & uniqueS[,2]==k),])[,3]))
        nkj_woi <- length(which(sy==j & sm ==k))
        likeregm <- 1
        for ( r in 1:q)
        {
          temp <- dnorm(M[i,r], dot(matX[i,], beta_m[which(beta_m[,1]==r & beta_m[,2]==j & beta_m[,3]==k),-c(1:3)]),sqrt(sig_m[which(sig_m[,1]==r & sig_m[,2]==j & sig_m[,3]==k),4]))
          likeregm <- likeregm*temp
        }
        
        probs[ntot+count1] <- (nj_woi/(alpha+n))* (nkj_woi/(mu+nj_woi)) * (psai/(psai+nkj_woi))*likeregy*likeregm*h0i[i]
        count1 <- count1 + 1
      }
    }
    #====== Fill in probs for old Y but new M clusters and new x cluster============
    count2 <- 1
    for ( j in 1: numy)
    {
      sub11 <- uniqueS[which(uniqueS[,1]==j),]
      if(nrow(as.matrix(t(sub11)))==1)
      {
        nummj <- length(unique(sub11[2]))
      }else{
        nummj <- length(unique(sub11[,2]))
      }
      nj_woi <- length(which(sy==j))
      likeregy <- dnorm(y[i], dot(full_mat[i,], beta_y[j,-c(1)]),sqrt(sig2_y[j,2]))
      probs[ntot+numm+count2] <- (nj_woi/(alpha+n))* (mu/(mu+nj_woi))*likeregy*h0m[i]*h0i[i]
      count2 <- count2+1
    }
    #====== Fill in probs for new Y cluster============
    probs[totalposs] <- (alpha/(alpha+n)) * h0y[i] * h0m[i] * h0i[i]
    
    #======choose new cluster==========================
    
    #Note: frame, frame2 and frame3 will help in choosing cluster
    
    frame <- numeric()
    for ( j in 1: length(unique(uniqueS[,1])))
    {
      sub11 <- uniqueS[which(uniqueS[,1]==j),]
      if(nrow(as.matrix(t(sub11)))==1)
      {
        temp <- length(unique(sub11[2]))
      }else{
        temp <- length(unique(sub11[,2]))
      }
      temp1 <- c(j,temp)
      frame <- rbind(frame,temp1)
    }
    frame <- data.frame(frame)
    colnames(frame) <- c("c", "sc_Number")  
    frame$cum <- cumsum(frame$sc_Number)
    
    frame2 <- numeric()
    for ( j in 1: numy)
    {
      sub11 <- uniqueS[which(uniqueS[,1]==j),]
      if(nrow(as.matrix(t(sub11)))==1)
      {
        nummj <- length(unique(sub11[2]))
      }else{
        nummj <- length(unique(sub11[,2]))
      }
      for ( k in 1: nummj)
      {
        sub22 <- uniqueS[which(uniqueS[,1]==j & uniqueS[,2]==k),]
        if(nrow(as.matrix(t(sub22)))==1)
        {
          numx_kj <- length(unique(sub22[3]))
        }else{
          numx_kj <- length(unique(sub22[,3]))
        }
        temp4 <- c(j,k, numx_kj)
        frame2 <- rbind(frame2, temp4)
      }
    }
    frame2 <- data.frame(frame2)
    colnames(frame2) <- c("c", "s", "ss")  
    
    frame3 <- numeric()
    for ( j in 1: numy)
    {
      sub11 <- uniqueS[which(uniqueS[,1]==j),]
      if(nrow(as.matrix(t(sub11)))==1)
      {
        nummj <- length(unique(sub11[2]))
      }else{
        nummj <- length(unique(sub11[,2]))
      }
      temp5 <- c(j,nummj)
      frame3 <- rbind(frame3,temp5)
    }
    frame3 <- data.frame(frame3)
    colnames(frame3) <- c("c", "s")  
    
    
    pro <- exp(log(probs, base=exp(1)) - max(log(probs,base=exp(1))))
    newcluster <- rmultinomf(pro)
    #print(newcluster)
    #============ Existing clusters ( To update sy, sm and sx)============
    if(newcluster < ( ntot+1))
    {
      sy <- insertvec( sy, i,uniqueS[newcluster,1])
      sm <- insertvec(sm, i, uniqueS[newcluster,2])
      sx <- insertvec( sx, i, uniqueS[newcluster,3])
    }  
    #============= Old y and old m cluster, but new x cluster (To update sy,sm, sx, xmu, xsig, xpi and uniqueS)
    if (ntot < newcluster & newcluster < (ntot+numm+1))
    {
      t <- newcluster - ntot
      tt <- frame$cum
      sy <- insertvec(sy, i, which((frame$cum)>=t)[1])
      if(sy[i]==1)
      {
        sm <- insertvec(sm, i, t)
      }else{
        sm <- insertvec(sm, i, t - tt[sy[i]-1])
      }
      sx <- insertvec( sx, i, frame2$ss[which(frame2$c==sy[i] & frame2$s==sm[i])] + 1 )
      newsigpars <- numeric()
      newmupars <- numeric()
      newpipars <- numeric()
      for(j in 1:(p1+ptx))
      {
        temp8 <- rbeta(1,X[i,j]+a0, 1-X[i,j]+b0)
        newpipars <- rbind(newpipars, temp8)
      }
      for ( j in 1:p2)
      {
        temp6 <- updtvar(X[i,p1+ptx+j], nu0,tau0,c0,mu0)
        temp7 <- updtmean(X[i,p1+ptx+j],temp6,c0,mu0)
        newsigpars <- rbind(newsigpars,temp6)
        newmupars <- rbind(newmupars,temp7)
      }
      position <- max(which(uniqueS[,1]==sy[i] & uniqueS[,2]==sm[i])) + 1
      xMuPars <- insertmat(xMuPars, position, c(sy[i], sm[i], sx[i],t(newmupars)))
      xSigPars <- insertmat(xSigPars, position, c(sy[i], sm[i], sx[i],t(newsigpars)))
      xPiPars <- insertmat(xPiPars, position, c(sy[i], sm[i], sx[i],t(newpipars)))
      uniqueS <- insertmat(uniqueS, position, c(sy[i], sm[i], sx[i]))
    }  
    #=============Old y cluster, new M and new x cluster (To update sy, sm, sx, beta_m, sig_m, xmu, xsig, xpi, uniqueS)  
    if ((ntot+numm) < newcluster & newcluster < (ntot+numm+numy+1)) 
    {
      sy <- insertvec(sy, i, newcluster - ntot - numm)
      sm <- insertvec(sm, i, frame3$s[which(frame3$c==sy[i])] + 1)
      sx <- insertvec(sx, i, 1)
      for ( r in 1:q)
      {
        newsig <- nwrgvar(matX[i,],M[i,r], betaa0, betab0,bs1[[r]])
        newbeta <- nwbet(matX[i,],M[i,r], newsig, bs1[[r]])
        to_include_sig <- c(r,sy[i], sm[i], newsig)
        to_include_beta <- c(r,sy[i], sm[i],newbeta)
        position <- max(which(sig_m[,1]==r & sig_m[,2]==sy[i]))+1 
        sig_m <- insertmat(sig_m,position, to_include_sig)
        beta_m <- insertmat(beta_m,position, to_include_beta)
      }
      newsigpars <- numeric()
      newmupars <- numeric()
      newpipars <- numeric()
      for(j in 1:(p1+ptx))
      {
        temp10 <- rbeta(1,X[i,j]+a0, 1-X[i,j]+b0)
        newpipars <- rbind(newpipars, temp10)
      }
      for ( j in 1:p2)
      {
        temp11 <- updtvar(X[i,p1+ptx+j], nu0,tau0,c0,mu0)
        temp12 <- updtmean(X[i,p1+ptx+j],temp11,c0,mu0)
        newsigpars <- rbind(newsigpars,temp11)
        newmupars <- rbind(newmupars,temp12)
      }
      position <- max(which(uniqueS[,1]==sy[i])) + 1
      xMuPars <- insertmat(xMuPars, position, c(sy[i], sm[i], sx[i],t(newmupars)))
      xSigPars <- insertmat(xSigPars, position, c(sy[i], sm[i], sx[i],t(newsigpars)))
      xPiPars <- insertmat(xPiPars, position, c(sy[i], sm[i], sx[i],t(newpipars)))
      uniqueS <- insertmat(uniqueS, position, c(sy[i], sm[i], sx[i]))
    }  
    #======================Completely new cluster=============================================  
    if (newcluster==ntot+numm+numy+1) 
    {
      sy <- insertvec(sy, i, numy +1)
      sm <- insertvec(sm, i, 1)
      sx <- insertvec(sx, i, 1)
      newsig <- nwrgvar(full_mat[i,],y[i], betaa0, betab0, bs)
      newbeta <- nwbet(full_mat[i,],y[i], newsig, bs)
      position <- nrow(beta_y) + 1
      to_include_sig <- c(sy[i],newsig)
      to_include_beta <- c(sy[i],newbeta)
      sig2_y <- insertmat(sig2_y, position, to_include_sig)
      beta_y <- insertmat(beta_y, position, to_include_beta)
      for ( r in 1:q)
      {
        newsig <- nwrgvar(matX[i,],M[i,r], betaa0, betab0, bs1[[r]])
        newbeta <- nwbet(matX[i,],M[i,r], newsig,bs1[[r]])
        to_include_sig <- c(r,sy[i], sm[i], newsig)
        to_include_beta <- c(r,sy[i], sm[i],newbeta)
        position <- max(which(sig_m[,1]==r))+1 
        sig_m <- insertmat(sig_m,position, to_include_sig)
        beta_m <- insertmat(beta_m,position, to_include_beta)
      }
      newsigpars <- numeric()
      newmupars <- numeric()
      newpipars <- numeric()
      for(j in 1:(p1+ptx))
      {
        temp10 <- rbeta(1,X[i,j]+a0, 1-X[i,j]+b0)
        newpipars <- rbind(newpipars, temp10)
      }
      for ( j in 1:p2)
      {
        temp11 <- updtvar(X[i,p1+ptx+j], nu0,tau0,c0,mu0)
        temp12 <- updtmean(X[i,p1+ptx+j],temp11,c0,mu0)
        newsigpars <- rbind(newsigpars,temp11)
        newmupars <- rbind(newmupars,temp12)
      }
      position <- nrow(uniqueS) + 1
      xMuPars <- insertmat(xMuPars, position, c(sy[i], sm[i], sx[i],t(newmupars)))
      xSigPars <- insertmat(xSigPars, position, c(sy[i], sm[i], sx[i],t(newsigpars)))
      xPiPars <- insertmat(xPiPars, position, c(sy[i], sm[i], sx[i],t(newpipars)))
      uniqueS <- insertmat(uniqueS, position, c(sy[i], sm[i], sx[i]))
      
    }
    
    multr <- c(multr, newcluster) 
    check[[i]] <- uniqueS 
    s_check[[i]] <- cbind(sy,sm,sx)
  }
  #====loop for "i"th person ends here=========
  cluster_res <- list(sy=sy, sm=sm, sx=sx,uniques=uniqueS,beta_y=beta_y, sig2_y=sig2_y, beta_m=beta_m, sig_m=sig_m,xmupars=xMuPars, xsigpars=xSigPars, xpipars=xPiPars, multr=multr, check =check, s_check = s_check)
  cluster_res
}

