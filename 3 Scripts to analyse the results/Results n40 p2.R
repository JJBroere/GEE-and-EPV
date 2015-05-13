###Calculating preformance of GEE, GEE_backwardselct and PGEE###
library(gee)
library(geepack)

#Obtaining estimates for GEE
EstimatesGEE <- matrix(1,2,100)
AA<-c(1,2)
for(i in 1:100){
  AAB<-summary(geeglm(outcome ~ V2+V3, id= id, data=Sim100N40P2u[[i]], corstr="exchangeable", family = binomial))
  EstimatesGEE[,i] <- AAB$coefficients[2:3,1]
}

#Obtaining estimates for GEE with backward selection
#Modified script, general script does not work because of models with zero selections, 
#general script does not accout for this to happen 
p<-2
ids <-matrix(NA,100,1)
for(i in 1:100){
  X.temp<<- as.matrix(Sim100N40P2u[[i]][,2:(p)])
  colnames(X.temp)<-c("V02")
  Y <<- as.matrix(Sim100N40P2u[[i]][,(p+2)])
  glmFull <- geeglm(Y ~ X.temp, id= id, Sim100N40P2u[[i]], corstr="exchangeable", family = binomial)
  ids[i,] <- coef(summary(glmFull))[-1, 4]
}

BB<-which(ids<0.1)

EstimatesBACK <- matrix(NA,2,100)
for(i in BB){
  AAB<-summary(geeglm(outcome ~ V2, id= id, data=Sim100N40P2u[[i]], corstr="exchangeable", family = binomial))
  EstimatesBACK[,i] <- AAB$coefficients[2:3,1]
}
AA<-c(17,29,42,50,56,65,68,69,72,76)
for(i in AA){
  AAB<-summary(geeglm(outcome ~ V2+V3, id= id, data=Sim100N40P2u[[i]], corstr="exchangeable", family = binomial))
  EstimatesBACK[,i] <- AAB$coefficients[2:3,1]
}

BABA <- matrix(NA, 100,2)
for(i in AA){
  AAB<-summary(geeglm(outcome ~ V2+V3, id= id, data=Sim100N40P2u[[i]], corstr="exchangeable", family = binomial))
  BABA[i,] <- coef(AAB)[-1, 4]
}
EstimatesBACK[1,17] <-NA
EstimatesBACK[1,72] <-NA
EstimatesBACK[1,76] <-NA

EstimatesBACK[is.na(EstimatesBACK)]<- 0

#Calculating MSE for GEE
MEGEE <- matrix(1,100)
BetaTrue <- c(0.5,0)
for(i in 1:100){
  MEGEE[i] <- t(EstimatesGEE[,i]-BetaTrue)%*%(EstimatesGEE[,i]-BetaTrue)
}

#Calculating MSE GEE_Backward selection
MEGEEBACK <- matrix(1,100)
BetaTrue <- c(0.5,0)
for(i in 1:100){
  MEGEEBACK[i] <- t(EstimatesBACK[,i]-BetaTrue)%*%(EstimatesBACK[,i]-BetaTrue)
}

#Calculating MSE PGEE
MESCAD <- matrix(1,100)
BetaTrue <- c(0.5,0)
for(i in 1:100){
  MESCAD[i] <- t(EstimatesSCAD4002[,i]-BetaTrue)%*%(EstimatesSCAD4002[,i]-BetaTrue)
}

#Mean MSE
mean(MEGEE)
mean(MEGEEBACK)
mean(MESCAD)
sd(MEGEE)/sqrt(100)
sd(MEGEEBACK)/sqrt(100)
sd(MESCAD)/sqrt(100)

#Underselction
UnderselectBACK <- mean(EstimatesBACK[1,]==0)
UnderselectPGEE <- mean(EstimatesSCAD4002[1,]==0)

#Overselection
OverselectBACK <- 1-mean(EstimatesBACK[2,]==0)
OverselectPGEE <- 1-mean(EstimatesSCAD4002[2,]==0)
UnderselectBACK
UnderselectPGEE
OverselectBACK
OverselectPGEE

#Calculating parameter Bias for GEE
RM <-rowMeans(EstimatesGEE)
Bias <- (RM-BetaTrue)/BetaTrue
Bias

#Calculating parameter Bias for GEEBACK
EstimatesGEEBackWithout <- EstimatesBACK
EstimatesGEEBackWithout[EstimatesGEEBackWithout==0] <- NA
RM <-rowMeans(EstimatesGEEBackWithout, na.rm=T)
Bias <- (RM-BetaTrue)/BetaTrue
Bias

#Calculating parameter Bias for PGEE
EstimatesPGEEWithout <- EstimatesSCAD4002
EstimatesPGEEWithout[EstimatesPGEEWithout==0] <- NA
RM <-rowMeans(EstimatesPGEEWithout, na.rm =T)
Bias <- (RM-BetaTrue)/BetaTrue
Bias



###Detemine the covarage of GEE###
Coefbeta1sim <-matrix(1,100)
Coefbeta2sim <-matrix(1,100)

Coefbeta1simp<-matrix(1,100)
Coefbeta2simp<-matrix(1,100)

COVE1<-matrix(1,100)
COVE2<-matrix(1,100)

for(j in 1:100){
  
  ####Analysing simulated data####
  GEEX1coef<- summary(geeglm(outcome ~ V2+V3, id= id, data=Sim100N40P2u[[j]], corstr="exchangeable", family = binomial))
  
  #Saving coefficients beta0 and beta1 for each replicated dataset
  
  Coefbeta1sim[j] <- GEEX1coef$coefficients[2,1]
  Coefbeta2sim[j] <- GEEX1coef$coefficients[3,1]
  
  
  Coefbeta1simp[j] <- GEEX1coef$coefficients[2,2]
  Coefbeta2simp[j] <- GEEX1coef$coefficients[3,2]
    
  COVE1[j] <- (Coefbeta1sim[j] - 1.96*Coefbeta1simp[j] <= 0.5 & 0.5 >=  Coefbeta1sim[j] - 1.96*Coefbeta1simp[j])
  COVE2[j] <- (Coefbeta2sim[j] - 1.96*Coefbeta2simp[j] <= 0 & 0 >=  Coefbeta2sim[j] - 1.96*Coefbeta2simp[j])
  
}


mean(COVE1)
mean(COVE2)


