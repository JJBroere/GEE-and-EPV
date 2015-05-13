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
               0.3,  0.3, 0.3,  0.3, 1), ncol=5)
  
  V <- diag(25)
  
  #Start simulation
  simulation100 <- for(j in 1:Nsim){
    
    setWinProgressBar(progressbar, j, title=paste( round(j/Nsim*100, 0), "% complete"))
    
    #Objects of storage
    Nsamp <- 100
    NP<- 25
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
    BETAij <- c(0,0.8,-0.5,-0.3,-0.3,0.3,-0.3,0.3,-0.3,0.3,-0.3,0.3,0.3,0,0,0,0,0,0,0,0,0,0,0,0,0)
    
    #Defing probabilities of scoring a one
    piij <- exp(Xij%*%BETAij)/(1+exp(Xij%*%BETAij))
    MU <- piij
   
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

Sim100N100P25u <- Simdata(100,100,25)