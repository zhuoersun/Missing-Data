# 1. Application of the methodologies on 2010/2011 YSS data;
# 2. Focus on Asian students (Grade 6-8), with attributes esteem, sex, marks
#    smoking status, BMI. Other subjects and features are deleted;
# 3. The filtered data set is stored in "yss.csv", with columns:
#       (i). provID: province ID. Not in the model. Deleted in the code
#      (ii). GRADE: grade of the student. In range 6 to 8
#     (iii). SEX: gender of the student
#      (iv). GETHNCC1: race of the student. All = 1 (Asian). Deleted in the code
#       (v). OMARKSA1: marks during the last year
#      (vi). DVTY1ST: smoking status
#     (vii). DVSELF: self-esteem score
#    (viii). BMI: BMI
# 4. Full information of the attributes is in the codebook file;
# 5. The original complete data: "yss10_data_public_ver3_20120125.csv"
# 6. "yss.data" has 1715 rows. Valus=99 means missingness. We take a subset
#    of this data set with n=493, with 144 cases BMI missing, by using random
#    seed = 123+349 = 472;
# 7. The results include:
#       (i). beta_hat_: point estimate
#      (ii). SE_: standard errors by sand-wich formula with correction matrix
#     (iii). p.value_: p-values
#    with _ = 1.CC  2.PIP  3.PIPA  4.A (proposed estimator) 





library(np)
data <- read.csv("yss.csv",header=TRUE)
data <- data[,-c(1,4)]
data$SEX <- data$SEX - 1
colnames(data) <- c("grade","sex","marks","smoke","esteem","BMI")
data <- data[(data$marks!=99),]
data <- data[(data$esteem!=99),]
data$BMI[data$BMI==99] <- NA
data$smoke[data$smoke!=3] <- 1
data$smoke[data$smoke==3] <- 0
#data$marks[data$marks<=3] <- 0
#data$marks[data$marks>=4] <- 1
data$ind <- 0
data$ind[!is.na(data$BMI)] <- 1

set.seed(123+349)
mis <- data[is.na(data$BMI),]
obs <- data[!is.na(data$BMI),]
n <- 493
n_mis <- 144
n_obs <- n-n_mis
index1 <- sample(c(1:length(obs$sex)),n_obs)
index2 <- sample(c(1:length(mis$sex)),n_mis)
newdata <- rbind(mis[index2,],obs[index1,])
yy <- newdata$esteem
xx <- newdata$BMI
xx[is.na(xx)] = -999
zz1 <- newdata$sex
zz2 <- newdata$marks
zz3 <- newdata$smoke
ww <- cbind(rep(1,n),xx,zz1,zz2,zz3)
rr <- newdata$ind



###CC#######################################################
S1=0
S2=0
for (i in 1:n){
  W=ww[i,]
  y=yy[i]
  R=rr[i]
  S1=S1+R*W%*%t(W)
  S2=S2+R*W*y
}
beta_hat1 <- solve(S1)%*%S2
Ubar <- 0
for (i in 1:n){
  W=ww[i,]
  y=yy[i]
  R=rr[i]
  U <- R*W*(as.numeric(y-t(W)%*%beta_hat1))
  Ubar <- Ubar + U/n
}
Gn=0
Vn=0
for (i in 1:n){
  W=ww[i,]
  y=yy[i]
  R=rr[i]
  U <- R*W*(as.numeric(y-t(W)%*%beta_hat1))
  Gn=Gn+(1/n)*R*W%*%t(W)
  Vn=Vn+(1/n)*(U-Ubar)%*%t(U-Ubar)
}
Cov_beta1 <- solve(Gn)%*%Vn%*%solve(Gn)
SE1 <- sqrt(diag(Cov_beta1)/n)
p.value1 <- (1-pt(abs(beta_hat1)/SE1,df=n-5))*2



###PIP#######################################################
alpha.fit<-glm(rr~yy+zz1+zz2+zz3,family=binomial())
alpha_hat<-alpha.fit$coeff
P_hat<-predict(alpha.fit,type="response")
P_hat[P_hat <= 1e-4] <- 1e-4

S1=0
S2=0
for (i in 1:n){
  W=ww[i,]
  y=yy[i]
  R=rr[i]
  S1=S1+(R/P_hat[i])*W%*%t(W)
  S2=S2+(R/P_hat[i])*W*y
}
beta_hat2=solve(S1)%*%S2

Gn=0
Vn=0
En=0
An=0
for (i in 1:n){
  W=ww[i,]
  y=yy[i]
  R=rr[i]
  Q=c(1,y,zz1[i],zz2[i],zz3[i])
  Gn=Gn+(1/n)*(R/P_hat[i])*W%*%t(W)
  Vn=Vn+(1/n)*(R/(P_hat[i])^2)*W%*%t(W)*(as.numeric((y-t(W)%*%beta_hat2)^2))
  En=En+(1/n)*W%*%t(Q)*(R/P_hat[i]-R)*(as.numeric(y-t(W)%*%beta_hat2))
  An=An+(1/n)*Q%*%t(Q)*P_hat[i]*(1-P_hat[i])
}
Cov_beta2 <- solve(Gn)%*%(Vn-En%*%solve(An)%*%t(En))%*%t(solve(Gn))
SE2 <- sqrt(diag(Cov_beta2)/n)
p.value2 <- (1-pt(abs(beta_hat2)/SE2,df=n-5))*2



###PIPA#######################################################
alpha.fit<-glm(rr~yy+zz1+zz2+zz3,family=binomial())
alpha_hat<-alpha.fit$coeff
P_hat<-predict(alpha.fit,type="response")
P_hat[P_hat <= 1e-4] <- 1e-4

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
  G1[G1==0] <- 1e-4
  G2[G2==0] <- 1e-4
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
    S1=S1+(R/P_hat[i])*W%*%t(W)+(1-R/P_hat[i])*phi2
    S2=S2+(R/P_hat[i])*W*y+(1-R/P_hat[i])*phi1
  } 
  beta_hat=solve(S1)%*%S2
}

theta=c(1,-beta_hat[3],-beta_hat[4],-beta_hat[5])
Q=cbind(yy,zz1,zz2,zz3)
sindex=Q%*%theta
h=1*sd(sindex)*n^(-1/3)
G1 <- npksum(txdat=sindex,tydat=rr,bws=h)$ksum
G2 <- npksum(txdat=sindex,tydat=rr,bws=0.4*h)$ksum
G1[G1==0] <- 1e-4
G2[G2==0] <- 1e-4
Phi[,1] <- npksum(txdat=sindex,tydat=rr*xx,bws=h)$ksum/G1
Phi[,2] <- npksum(txdat=sindex,tydat=rr*(xx^2),bws=0.4*h)$ksum/G2
P_nw <- npksum(txdat=sindex,tydat=rr,bws=h)$ksum/npksum(txdat=sindex,bws=h)$ksum
P_nw[P_nw <= 1e-4] <- 1e-4

beta_hat3 = beta_hat

Gn=0
Vn=0
C_star=0
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
  S=W*(y-t(W)%*%beta_hat3)
  S_star=phi1-phi2%*%beta_hat3
  Gn=Gn+(1/n)*((R/P_hat[i])*W%*%t(W)+(1-R/P_hat[i])*phi2)
  Vn=Vn+(S-S_star)%*%t(S-S_star)*(R/(P_nw[i])^2)*(1/n)
  C_star=C_star+S_star%*%t(S_star)*(1/n)
}
Cov_beta3 <- solve(Gn)%*%(Vn+C_star)%*%t(solve(Gn))
SE3 <- sqrt(diag(Cov_beta3)/n)
p.value3 <- (1-pt(abs(beta_hat3)/SE3,df=n-5))*2



###A#######################################################
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
  G1[G1==0] <- 1e-4
  G2[G2==0] <- 1e-4
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

theta=c(1,-beta_hat[3],-beta_hat[4],-beta_hat[5])
Q=cbind(yy,zz1,zz2,zz3)
sindex=Q%*%theta
h=1*sd(sindex)*n^(-1/3)
G1 <- npksum(txdat=sindex,tydat=rr,bws=h)$ksum
G2 <- npksum(txdat=sindex,tydat=rr,bws=0.4*h)$ksum
G1[G1==0] <- 1e-4
G2[G2==0] <- 1e-4
Phi[,1] <- npksum(txdat=sindex,tydat=rr*xx,bws=h)$ksum/G1
Phi[,2] <- npksum(txdat=sindex,tydat=rr*(xx^2),bws=0.4*h)$ksum/G2
P_nw <- npksum(txdat=sindex,tydat=rr,bws=h)$ksum/npksum(txdat=sindex,bws=h)$ksum
P_nw[P_nw <= 1e-4] <- 1e-4

beta_hat4 = beta_hat

Gn=0
Vn=0
n_o=sum(rr)
G3 <- npksum(txdat=sindex,tydat=1-rr,bws=h)$ksum
G3[G3==0] <- 1e-4
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

Cov_beta4 <- M%*%solve(Gn)%*%(Vn)%*%t(solve(Gn))
SE4 <- sqrt(diag(Cov_beta4)/n)
p.value4 <- (1-pt(abs(beta_hat4)/SE4,df=n-5))*2


