# 1. For proposed estimator: beta_A;
# 2. No inverse probability, with augmentation using a single-index model
#    under MAR mechanism;
# 3. Need 'np' package;
# 4. Functions: generating data, point estimation, empirical/asymptotic 
#    SE calculation (Correction matrix applied), 95% coverages;
# 5. For simulations in the paper only. Fix dim(X)=1, dim(Z)=3;
# 6. Return a list:
#      (i). Beta: a MC*5 matrix storing the estimates of beta
#     (ii). bias: the averaged bias of MC simulations
#    (iii). epsd: empirical SD of MC simulations, based on Beta in (i)
#     (iv). asysd: sqrt of diagonal elements of the averaged Cov matrix over MC simulations
#      (v). cover: 95% coverages of beta

library(np)               

expit<-function(x){             
  y=exp(x)/(1+exp(x))
  return(y)
}

AWEE_MAR <- function(xd,inno,alpha,n,seed,trim){
  Beta=matrix(0,MC,5)
  AsySE=matrix(0,MC,5)
  Cov_beta=matrix(0,5,5)
  coverage<-rep(0,5)
  for (mc in 1:MC){
    set.seed(seed+mc)
    ww=matrix(0,nrow=n,ncol=5)
    xx=NULL
    xxt=NULL
    zz1=NULL
    zz2=NULL
    zz3=NULL
    rr=NULL
    yy=NULL
    P=NULL
    S1=0
    S2=0
    for (i in 1:n){
      X=(xd=="n")*rnorm(1)+(xd=="t")*rt(1,df=5)/sqrt(5/3)+(xd=="g")*(rgamma(1,shape=5,scale=1)-5)/sqrt(5)
      Z1=rnorm(1)
      Z2=rnorm(1)
      Z3=rnorm(1)
      W=c(1,X,Z1,Z2,Z3)
      innov=(inno=="n")*rnorm(1)+(inno=="t")*rt(1,df=5)/sqrt(5/3)+(inno=="g")*(rgamma(1,shape=5,scale=1)-5)/sqrt(5)
      y=t(W)%*%beta+innov
      p=expit(t(alpha)%*%c(1,y,Z1,Z2,Z3))
      R=rbinom(1,1,p)
      xxt=c(xxt,X)
      X=X*R-999*(1-R)
      W=c(1,X,Z1,Z2,Z3)
      ww[i,]=W
      yy=c(yy,y)
      xx=c(xx,X)
      zz1=c(zz1,Z1)
      zz2=c(zz2,Z2)
      zz3=c(zz3,Z3)
      rr=c(rr,R)
    }
    
    beta_hat=c(0,0,0,0,0)
    t=c(1,1,1,1,1)
    
    while (sum((beta_hat-t)^2)>1e-6){
      S1=0
      S2=0
      Phi=matrix(0,nrow=n,ncol=2)
      t=beta_hat
      theta=c(1,-beta_hat[3],-beta_hat[4],-beta_hat[5])
      Q=cbind(yy,zz1,zz2,zz3)
      sindex=Q%*%theta
      h=1*sd(sindex)*n^(-1/3)
      G1 <- npksum(txdat=sindex,tydat=rr,bws=h)$ksum
      G2 <- npksum(txdat=sindex,tydat=rr,bws=0.4*h)$ksum
      G1[G1==0] <- trim
      G2[G2==0] <- trim
      Phi[,1] <- npksum(txdat=sindex,tydat=rr*xx,bws=h)$ksum/G1
      Phi[,2] <- npksum(txdat=sindex,tydat=rr*(xx^2),bws=0.4*h)$ksum/G2
      for (i in 1:n){
        W=ww[i,]
        y=yy[i]
        R=rr[i]
        Z1=zz1[i]
        Z2=zz2[i]
        Z3=zz3[i]
        phi1=c(1,Phi[i,1],Z1,Z2,Z3)*y
        phi2=c(1,Phi[i,1],Z1,Z2,Z3)%*%t(c(1,Phi[i,1],Z1,Z2,Z3))
        phi2[2,2]=Phi[i,2]
        S1=S1+R*W%*%t(W)+(1-R)*phi2
        S2=S2+R*W*y+(1-R)*phi1
      } 
      beta_hat=solve(S1)%*%S2
    }
    Beta[mc,]=beta_hat
    
    #if (mc%%10==0) print(mc)
    
    theta=c(1,-beta_hat[3],-beta_hat[4],-beta_hat[5])
    Q=cbind(yy,zz1,zz2,zz3)
    sindex=Q%*%theta
    h=1*sd(sindex)*n^(-1/3)
    G1 <- npksum(txdat=sindex,tydat=rr,bws=h)$ksum
    G2 <- npksum(txdat=sindex,tydat=rr,bws=0.4*h)$ksum
    G1[G1==0] <- trim
    G2[G2==0] <- trim
    Phi[,1] <- npksum(txdat=sindex,tydat=rr*xx,bws=h)$ksum/G1
    Phi[,2] <- npksum(txdat=sindex,tydat=rr*(xx^2),bws=0.4*h)$ksum/G2
    P_nw <- npksum(txdat=sindex,tydat=rr,bws=h)$ksum/npksum(txdat=sindex,bws=h)$ksum
    P_nw[P_nw <= trim] <- trim
    eta_n <- (n*h^4+(n*h^2)^(-1))^(1/2)
    
    
    Gn=0
    Vn=0
    n_o=sum(rr)
    G3 <- npksum(txdat=sindex,tydat=1-rr,bws=h)$ksum
    G3[G3==0] <- trim
    kerwt <- npksum(txdat=sindex,bws=h,return.kernel.weights = TRUE)$kw
    phiz <- phizx <- matrix(0,nrow=n,ncol=3)
    for (i in 1:n){
      sum1 <- 0
      sum2 <- 0
      for (j in 1:n){
        W=ww[j,]
        X=xx[j]
        y=yy[j]
        R=rr[j]
        Z1=zz1[j]
        Z2=zz2[j]
        Z3=zz3[j]
        sum1 <- sum1 + (1-R)*c(Z1,Z2,Z3)*kerwt[i,j]
        sum2 <- sum2 + (1-R)*X*c(Z1,Z2,Z3)*kerwt[i,j]
      }
      phiz[i,] <- sum1/G3[i]
      phizx[i,] <- sum2/G3[i]
    }
    
    for (i in 1:n){
      W=ww[i,]
      y=yy[i]
      R=rr[i]
      X=xx[i]
      Z1=zz1[i]
      Z2=zz2[i]
      Z3=zz3[i]
      phi1=c(1,Phi[i,1],Z1,Z2,Z3)*y
      phi2=c(1,Phi[i,1],Z1,Z2,Z3)%*%t(c(1,Phi[i,1],Z1,Z2,Z3))
      phi2[2,2]=Phi[i,2]
      S=W*(y-t(W)%*%beta_hat)
      S_star=phi1-phi2%*%beta_hat
      Gn=Gn+(1/n)*(R*W%*%t(W)+(1-R)*phi2)
      S_stark=c(sindex[i]-beta_hat[1]-beta_hat[2]*Phi[i,1],Phi[i,1]*sindex[i]-Phi[i,1]*beta_hat[1]-beta_hat[2]*Phi[i,2],
                as.vector(phiz[i,])*(sindex[i]-beta_hat[1]-beta_hat[2]*Phi[i,1]))
      Sk=c(sindex[i]-beta_hat[1]-beta_hat[2]*X,X*(sindex[i]-beta_hat[1]-beta_hat[2]*X),as.vector(phiz[i,])*(sindex[i]-beta_hat[1]-beta_hat[2]*X))
      ck=(1/P_nw[i]-1)
      Un <- R*S + (1-R)*S_star + R*(Sk-S_stark)*ck
      Vn <- Vn + Un%*%t(Un)*(1/n)
    }
    #print(mc)
    
    fac <- 1-0.3*sum(1-rr)/n
    fac2 <- 1-0.7*sum(1-rr)/n
    M <- diag(c(1/fac,(1/fac2)*exp((500-n)/5000),rep(1/fac,3)))
    
    asycov <- M%*%solve(Gn)%*%(Vn)%*%t(solve(Gn))
    asyse <- sqrt(diag(asycov)/n)
    lower <- beta_hat-1.96*asyse
    upper <- beta_hat+1.96*asyse
    coverage <- coverage + (beta>=lower)*(beta<=upper)
    Cov_beta <- Cov_beta + asycov
  }
  asysd=sqrt(diag(Cov_beta/MC))
  bias=apply(Beta,2,mean)-beta
  epsd=apply(Beta,2,sd)/sqrt(MC)
  cover=coverage/MC
  return(list(Beta=Beta,bias=bias,epsd=epsd,asysd=asysd,cover=cover))
}



