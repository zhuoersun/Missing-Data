# 1. Simulations for normal data, 60% missing on average;
# 2. Sourced functions below should be in the same directory;
# 3. Sample size n=500, MC=1000 repeated simulations, random seed=1234;
# 4. X~standardized Gamma(5,1), Z~N_3(0,I_3), noise~standardized t(5);
# 5. Estimated seletion probabilities are bounded below from 1e-4 to avoid
#    invalid 'NaN' results;
# 6. Obtain:
#      (i). Beta_: a list storing all results of the called function
#     (ii). bias_: the averaged bias of MC simulations
#    (iii). epsd_: the empirical SE of MC simulations
#     (iv). asysd_: the averaged asymptotic SE over MC simulations
#      (v). cover_: 95% coverages of beta



source("full_MAR.R")
source("CC_MAR.R")
source("ip_noaug_MAR.R")
source("ip_aug_MAR.R")
source("AWEE_MAR.R")



beta=c(0,0.5,1,-1,-0.5)
MC=1000

#alpha<-c(3,-0.1,-0.1,0,0)
#alpha<-c(2.2,-0.9,-0.7,0,0)
#alpha<-c(0.5,-1,-0.5,0,0)
alpha<-c(-0.5,-0.5,-0.5,0,0)


####Full data analysis under MAR####
Beta0 <- full_MAR(xd="g",inno="t",alpha=alpha,n=500,seed=1234)
bias0 <- Beta0$bias
epsd0 <- Beta0$epsd
asysd0 <- (Beta0$asysd)/sqrt(MC*500)
cover0 <- Beta0$cover


####CC under MAR####
Beta1 <- CC_MAR(xd="g",inno="t",alpha=alpha,n=500,seed=1234)
bias1 <- Beta1$bias
epsd1 <- Beta1$epsd
asysd1 <- (Beta1$asysd)/sqrt(MC*500)
cover1 <- Beta1$cover


####IPW under MAR####
Beta2 <- ip_noaug_MAR(xd="g",inno="t",alpha=alpha,n=500,seed=1234,trim=1e-4)
bias2 <- Beta2$bias
epsd2 <- Beta2$epsd
asysd2 <- (Beta2$asysd)/sqrt(MC*500)
cover2 <- Beta2$cover


####PIPA under MAR####
Beta3 <- ip_aug_MAR(xd="g",inno="t",alpha=alpha,n=500,seed=1234,trim=1e-4)
bias3 <- Beta3$bias
epsd3 <- Beta3$epsd
cover3 <- Beta3$cover


####AWEE under MAR####
Beta5 <- AWEE_MAR(xd="g",inno="t",alpha=alpha,n=500,seed=1234,trim=1e-4)
bias5 <- Beta5$bias
epsd5 <- Beta5$epsd
asysd5 <- (Beta5$asysd)/sqrt(MC*500)
cover5 <- Beta5$cover




save.image(file="MC_nonnormal_n=500_4.RData")





