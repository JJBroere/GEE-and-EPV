####Joris Broere####
####Simulation script####
####Simulatiting correlated binomial data with gaussian predictors####
library(gee)
library(geepack)
#package 'mvtBinaryEP'not availeble via main CRAN website, can still be downloaded via CRAN archive
library(mvtBinaryEP)
library(mvtnorm)

set.seed(23)
Simdata <- function(Nsamp, Nsim, NP){
  
  #progressbar
  progressbar <- winProgressBar(title = "progress bar", min = 0, max = Nsim, width = 300)
  
  #Objects of storage
  MEANMU<-matrix(1,Nsim)
  MINMU<-matrix(1,Nsim)
  MAXMU<-matrix(1,Nsim)
  Datalist<- list(NA)
  
  #Used corrrlation matrix for binary outcomes
  r2<-matrix(c(1, 0.3, 0.3, 0.3, 0.3,
               0.3,  1, 0.3,  0.3, 0.3,
               0.3,  0.3, 1,  0.3, 0.3,
               0.3, 0.3, 0.3, 1, 0.3,
               0.3,  0.3, 0.3,	0.3, 1), ncol=5)
  
  #Corrrlation matrix for Predictor variables
#   times <- 1:NP
#   rho <- 0.5
#   sigma <- 1
#   H <- abs(outer(times, times, "-"))
#   V <- sigma * rho^H
#   p <- nrow(V)
#   V[cbind(1:p, 1:p)] <- V[cbind(1:p, 1:p)] * sigma
  
  V <- diag(4)
  
  
  #r1.1_mat <- matrix(c(1, 0.5, 0.25, 0.125, 
             #   0.5,  1, 0.5,  0.25, 
              #  0.25,  0.5, 1,  0.5,
               # 0.125, 0.25, 0.5, 1), ncol=4)

  
 # V[1:4,1:4] <- r1.1_mat
  #V[5:8,5:8] <- r1.1_mat
  #V[9:12,9:12] <- r1.1_mat
  
    #tart simulation
  simulation100 <- for(j in 1:Nsim){
    
    setWinProgressBar(progressbar, j, title=paste( round(j/Nsim*100, 0), "% complete"))
    
    #Objects of storage
    
    x1 <-matrix(1,Nsamp)
      
    Dataset <-matrix(1,Nsamp,5)
    
    #generating random data
    Xi <- matrix(rmvnorm(n=Nsamp,sigma = V), ncol=NP)
    Xijr <- matrix(rep(Xi, each=5), ncol=NP)
    
    #making a ID variable
    id <- rep(1:Nsamp, each =5)
    
    # Combining information
    x0 <- rep(1,Nsamp)
    Xij <- cbind(x0, Xi)
    
    #Defing population parameters
    #BETAij <- c(0,0.8,0.8,0.8,0.8,-0.5,-0.5,-0.5,-0.5,-0.3,-0.3,-0.3,-0.3,0,0,0,0,0,0,0,0)
    #BETAij <- c(0,0.5,0.5,0.5,-0.5,-0.5,-0.5,0.3,0.3,-0.3,-0.3,0,0,0,0,0,0,0,0,0,0)
    BETAij <- c(0,0.8,-0.5,-0.3,0)
    #Defing probabilities of scoring a one
    piij <- exp(Xij%*%BETAij)/(1+exp(Xij%*%BETAij))
    MU <- piij
    mean(MU)
    min(MU)
    max(MU)
    MEANMU[j]<-mean(MU)
    MINMU[j]<-min(MU)
    MAXMU[j]<-max(MU)
    
    #Generating binary outcomes, based on previous defined probabilities
    Binmatrix <- for(i in 1:Nsamp){
      
      #simulating outcomes, based on different probabilties for every simulated individual i 
      Mu3 <- c(MU[i], MU[i], MU[i], MU[i], MU[i])
      Simdat1 <- ep(mu=Mu3, R=r2, nRep=1, seed=NULL)
      Dataset[i,] <- Simdat1$y
    }
    
    outcome <- as.vector(t(Dataset))
    AA<-cbind(id, Xijr,outcome)
    BB <- as.data.frame(AA)
    Datalist[[j]] <- BB
    
  }
  close(progressbar)
  return(Datalist)
}
warnings()
Sim100N40P4u <- Simdata(40,100,4)
setwd("C:/Users/Joris/Google Drive/School/Methodology & Statistics/Master scriptie/Data archief/Gebruik map")
save(Sim100N40P4u,file="Sim100N40P4u.Rda")

summary(geeglm(outcome ~ V2+V3+V4+V5+V6+V7+V8+V9+V10+V11, id= id, data=Sim100N40P10[[16]], corstr="exchangeable", family = binomial))



setwd("C:/Users/Joris/Dropbox/School/Methodology & Statistics/Master scriptie/Data/Simulated data")
save(Sim107a0,file="Sim107a0.Rda")



Coefbeta0sim <-matrix(1,1000)
Coefbeta1sim <-matrix(1,1000)
Coefbeta2sim <-matrix(1,1000)
Coefbeta3sim <-matrix(1,1000)
Coefbeta4sim <-matrix(1,1000)
Coefbeta5sim <-matrix(1,1000)
Coefbeta6sim <-matrix(1,1000)
Coefbeta7sim <-matrix(1,1000)
Coefbeta8sim <-matrix(1,1000)

Coefbeta0simp<-matrix(1,1000)
Coefbeta1simp<-matrix(1,1000)
Coefbeta2simp<-matrix(1,1000)
Coefbeta3simp<-matrix(1,1000)
Coefbeta4simp<-matrix(1,1000)
Coefbeta5simp<-matrix(1,1000)
Coefbeta6simp<-matrix(1,1000)
Coefbeta7simp<-matrix(1,1000)
Coefbeta8simp<-matrix(1,1000)

Coefbeta0simpp <-matrix(1,1000)
Coefbeta1simpp <-matrix(1,1000)
Coefbeta2simpp <-matrix(1,1000)
Coefbeta3simpp <-matrix(1,1000)
Coefbeta4simpp <-matrix(1,1000)
Coefbeta5simpp <-matrix(1,1000)
Coefbeta6simpp <-matrix(1,1000)
Coefbeta7simpp <-matrix(1,1000)
Coefbeta8simpp <-matrix(1,1000)

COVE0<-matrix(1,1000)
COVE1<-matrix(1,1000)
COVE2<-matrix(1,1000)
COVE3<-matrix(1,1000)
COVE4<-matrix(1,1000)
COVE5<-matrix(1,1000)
COVE6<-matrix(1,1000)
COVE7<-matrix(1,1000)
COVE8<-matrix(1,1000)

for(j in 1:1000){
  
  ####Analysing simulated data####
  GEEX1coef<- summary(geeglm(outcome ~ V3+V4+V5+V6+V7+V8+V9+V10, id= id, data=Sim80a0[[j]], corstr="exchangeable", family = binomial))
  
  #Saving coefficients beta0 and beta1 for each replicated dataset
  #   Coefbeta0sim[j] <- GEEX1coef$coefficients[1,1]
  #   Coefbeta1sim[j] <- GEEX1coef$coefficients[2,1]
  #   Coefbeta2sim[j] <- GEEX1coef$coefficients[3,1]
  #   Coefbeta3sim[j] <- GEEX1coef$coefficients[4,1]
  #   Coefbeta4sim[j] <- GEEX1coef$coefficients[5,1]
  #   Coefbeta5sim[j] <- GEEX1coef$coefficients[6,1]
  #   Coefbeta6sim[j] <- GEEX1coef$coefficients[7,1]
  #   Coefbeta7sim[j] <- GEEX1coef$coefficients[8,1]  
  #   Coefbeta8sim[j] <- GEEX1coef$coefficients[9,1]
  #   
  #   Coefbeta0simp[j] <- GEEX1coef$coefficients[1,2]
  #   Coefbeta1simp[j] <- GEEX1coef$coefficients[2,2]
  #   Coefbeta2simp[j] <- GEEX1coef$coefficients[3,2]
  #   Coefbeta3simp[j] <- GEEX1coef$coefficients[4,2]
  #   Coefbeta4simp[j] <- GEEX1coef$coefficients[5,2]
  #   Coefbeta5simp[j] <- GEEX1coef$coefficients[6,2]
  #   Coefbeta6simp[j] <- GEEX1coef$coefficients[7,2]
  #   Coefbeta7simp[j] <- GEEX1coef$coefficients[8,2]  
  #   Coefbeta8simp[j] <- GEEX1coef$coefficients[9,2]
  #   
  # COVE0[j] <- (Coefbeta0sim[j] - 1.96*Coefbeta0simp[j] <= -1 & -1 >=  Coefbeta0sim[j] - 1.96*Coefbeta0simp[j])
  # COVE1[j] <- (Coefbeta1sim[j] - 1.96*Coefbeta1simp[j] <= 0.3 & 0.3 >=  Coefbeta1sim[j] - 1.96*Coefbeta1simp[j])
  # COVE2[j] <- (Coefbeta2sim[j] - 1.96*Coefbeta2simp[j] <= 0.5 & 0.5 >=  Coefbeta2sim[j] - 1.96*Coefbeta2simp[j])
  # COVE3[j] <- (Coefbeta3sim[j] - 1.96*Coefbeta3simp[j] <= -0.4 & -0.4 >=  Coefbeta3sim[j] - 1.96*Coefbeta3simp[j])
  # COVE4[j] <- (Coefbeta4sim[j] - 1.96*Coefbeta4simp[j] <= -0.3 & -0.3 >=  Coefbeta4sim[j] - 1.96*Coefbeta4simp[j])
  # COVE5[j] <- (Coefbeta5sim[j] - 1.96*Coefbeta5simp[j] <=  0 & 0 >=  Coefbeta5sim[j] - 1.96*Coefbeta5simp[j])
  # COVE6[j] <- (Coefbeta6sim[j] - 1.96*Coefbeta6simp[j] <= 0.2 & 0.2 >=  Coefbeta6sim[j] - 1.96*Coefbeta6simp[j])
  # COVE7[j] <- (Coefbeta7sim[j] - 1.96*Coefbeta7simp[j] <= 0 & 0 >=  Coefbeta7sim[j] - 1.96*Coefbeta7simp[j])
  # COVE8[j] <- (Coefbeta8sim[j] - 1.96*Coefbeta8simp[j] <= 0 & 0 >=  Coefbeta8sim[j] - 1.96*Coefbeta8simp[j])
  
  Coefbeta0simpp[j] <- GEEX1coef$coefficients[1,4]
  Coefbeta1simpp[j] <- GEEX1coef$coefficients[2,4]
  Coefbeta2simpp[j] <- GEEX1coef$coefficients[3,4]
  Coefbeta3simpp[j] <- GEEX1coef$coefficients[4,4]
  Coefbeta4simpp[j] <- GEEX1coef$coefficients[5,4]
  Coefbeta5simpp[j] <- GEEX1coef$coefficients[6,4]
  Coefbeta6simpp[j] <- GEEX1coef$coefficients[7,4]
  Coefbeta7simpp[j] <- GEEX1coef$coefficients[8,4]  
  Coefbeta8simpp[j] <- GEEX1coef$coefficients[9,4]
  
}

# Meansim100n200beta0 <- mean(Coefbeta0sim)
# Meansim100n200beta1 <- mean(Coefbeta1sim)
# Meansim100n200beta2 <- mean(Coefbeta2sim)
# Meansim100n200beta3 <- mean(Coefbeta3sim)
# Meansim100n200beta4 <- mean(Coefbeta4sim)
# Meansim100n200beta5 <- mean(Coefbeta5sim)
# Meansim100n200beta6 <- mean(Coefbeta6sim)
# Meansim100n200beta7 <- mean(Coefbeta7sim)
# Meansim100n200beta8 <- mean(Coefbeta8sim)
# 
# Meansim100n200beta0
# Meansim100n200beta1
# Meansim100n200beta2
# Meansim100n200beta3
# Meansim100n200beta4
# Meansim100n200beta5
# Meansim100n200beta6
# Meansim100n200beta7
# Meansim100n200beta8
# 
# mean(COVE0)
# mean(COVE1)
# mean(COVE2)
# mean(COVE3)
# mean(COVE4)
# mean(COVE5)
# mean(COVE6)
# mean(COVE7)
# mean(COVE8)


Tp1<-Coefbeta1simpp < 0.05
Tp2<-Coefbeta2simpp < 0.05
Tp3<-Coefbeta3simpp < 0.05
Tp4<-Coefbeta4simpp < 0.05
Fp1<-Coefbeta5simpp < 0.05
Tp5<-Coefbeta6simpp < 0.05
Fp2<-Coefbeta7simpp < 0.05
Fp3<-Coefbeta8simpp < 0.05

TP <- sum(Tp1, Tp2, Tp3, Tp4, Tp5)
Fp <- sum(Fp1, Fp2, Fp3)
TP5 <- sum(Tp5)
TP
Fp
TP5
