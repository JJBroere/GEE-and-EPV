###Calculating preformance of GEE, GEE_backwardselct and PGEE###
library(gee)
library(geepack)

#Obtaining estimates for GEE
EstimatesGEE <- matrix(1,5,100)
for(i in 1:100){
  AAB<-summary(geeglm(outcome ~ V2+V3+V4+V5+V6, id= id, data=Sim100N100P5u[[i]], corstr="exchangeable", family = binomial))
  EstimatesGEE[,i] <- AAB$coefficients[2:6,1]
}

#Backward selecting GEE
p<-5
EstimatesB <- matrix(NA,5,100)
EstimatesBack <- matrix(NA,5,100)
row.names(EstimatesBack) <-c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05", "X.tempV06")
row.names(EstimatesB) <- c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05", "X.tempV06")
for(i in 1:100){
  X.temp<<- as.matrix(Sim100N100P5u[[i]][,2:(p+1)])
  colnames(X.temp)<-c("V02", "V03", "V04", "V05", "V06")
  Y <<- as.matrix(Sim100N100P5u[[i]][,(p+2)])
  glmFull <- geeglm(Y ~ X.temp, id= id, Sim100N100P5u[[i]], corstr="independence", family = binomial)
  ids <- which(coef(summary(glmFull))[-1, 4] < 0.5)
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N100P5u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.3)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N100P5u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.2)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N100P5u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.1)]
  X.temp <- X.temp[,ids]
  GB<- geeglm(Y ~ X.temp, id= id, Sim100N100P5u[[i]], corstr="independence", family = binomial)
  GBA <- GB$coefficients
  BBAA <-merge(EstimatesB[,i], GBA, by = "row.names", all.x = T,sort=T)
  EstimatesBack[,i] <- BBAA$y
}

EstimatesBack[is.na(EstimatesBack)]<- 0#Setting notselected varables to zere

#Calculating MSE for GEE
MEGEE <- matrix(1,100)
BetaTrue <- c(0.8,-0.5,-0.3,0,0)
for(i in 1:100){
  MEGEE[i] <- t(EstimatesGEE[,i]-BetaTrue)%*%(EstimatesGEE[,i]-BetaTrue)
}

#Calculating MSE GEE_Backward selection
MEGEEBACK <- matrix(1,100)
BetaTrue <- c(0.8,-0.5,-0.3,0,0)
for(i in 1:100){
  MEGEEBACK[i] <- t(EstimatesBack[,i]-BetaTrue)%*%(EstimatesBack[,i]-BetaTrue)
}

#Calculating MSE PGEE
MESCAD <- matrix(1,100)
BetaTrue <- c(0.8,-0.5,-0.3,0,0)
for(i in 1:100){
  MESCAD[i] <- t(EstimatesSCAD1005u[,i]-BetaTrue)%*%(EstimatesSCAD1005u[,i]-BetaTrue)
}

#Mean MSE
mean(MEGEE)
mean(MEGEEBACK)
mean(MESCAD)
sd(MEGEE)/sqrt(100)
sd(MEGEEBACK)/sqrt(100)
sd(MESCAD)/sqrt(100)

#Underselction
UnderselectBACK <- mean(EstimatesBack[1:3,]==0)
UnderselectPGEE <- mean(EstimatesSCAD1005u[1:3,]==0)
#Overselection
OverselectBACK <- 1-mean(EstimatesBack[4:5,]==0)
OverselectPGEE <- 1-mean(EstimatesSCAD1005u[4:5,]==0)
UnderselectBACK
UnderselectPGEE
OverselectBACK
OverselectPGEE

#Calculating parameter Bias for GEE
RM <-rowMeans(EstimatesGEE)
Bias <- (RM-BetaTrue)/BetaTrue
Bias

#Calculating parameter Bias for GEEBACK
EstimatesGEEBackWithout <- EstimatesBack
EstimatesGEEBackWithout[EstimatesGEEBackWithout==0] <- NA
RM <-rowMeans(EstimatesGEEBackWithout, na.rm=T)
Bias <- (RM-BetaTrue)/BetaTrue
Bias

#Calculating parameter Bias for PGEE
EstimatesPGEEWithout <- EstimatesSCAD1005u
EstimatesPGEEWithout[EstimatesPGEEWithout==0] <- NA
RM <-rowMeans(EstimatesPGEEWithout, na.rm =T)
Bias <- (RM-BetaTrue)/BetaTrue
Bias

###Detemine the covarage of GEE###
Coefbeta1sim <-matrix(1,100)
Coefbeta2sim <-matrix(1,100)
Coefbeta3sim <-matrix(1,100)

Coefbeta1simp<-matrix(1,100)
Coefbeta2simp<-matrix(1,100)
Coefbeta3simp<-matrix(1,100)

COVE1<-matrix(1,100)
COVE2<-matrix(1,100)
COVE3<-matrix(1,100)
COVE4<-matrix(1,100)
COVE5<-matrix(1,100)

for(j in 1:100){
  
  ####Analysing simulated data####
  GEEX1coef<- summary(geeglm(outcome ~ V2+V3+V4+V5+V6, id= id, data=Sim100N100P5u[[j]], corstr="exchangeable", family = binomial))
  #Saving coefficients beta0 and beta1 for each replicated dataset
  
  Coefbeta1sim[j] <- GEEX1coef$coefficients[2,1]
  Coefbeta2sim[j] <- GEEX1coef$coefficients[3,1]
  Coefbeta3sim[j] <- GEEX1coef$coefficients[4,1]
  Coefbeta4sim[j] <- GEEX1coef$coefficients[5,1]
  Coefbeta5sim[j] <- GEEX1coef$coefficients[6,1]

  
  Coefbeta1simp[j] <- GEEX1coef$coefficients[2,2]
  Coefbeta2simp[j] <- GEEX1coef$coefficients[3,2]
  Coefbeta3simp[j] <- GEEX1coef$coefficients[4,2]
  Coefbeta4simp[j] <- GEEX1coef$coefficients[5,2]
  Coefbeta5simp[j] <- GEEX1coef$coefficients[6,2]
  
  
  COVE1[j] <- (Coefbeta1sim[j] - 1.96*Coefbeta1simp[j] <= 0.8 & 0.8 >=  Coefbeta1sim[j] - 1.96*Coefbeta1simp[j])
  COVE2[j] <- (Coefbeta2sim[j] - 1.96*Coefbeta2simp[j] <= -0.5 & -0.5 >=  Coefbeta2sim[j] - 1.96*Coefbeta2simp[j])
  COVE3[j] <- (Coefbeta3sim[j] - 1.96*Coefbeta3simp[j] <= -0.3 & -0.3 >=  Coefbeta3sim[j] - 1.96*Coefbeta3simp[j])
  COVE4[j] <- (Coefbeta4sim[j] - 1.96*Coefbeta4simp[j] <= 0 & 0 >=  Coefbeta4sim[j] - 1.96*Coefbeta4simp[j])
  COVE5[j] <- (Coefbeta5sim[j] - 1.96*Coefbeta5simp[j] <=  0 & 0 >=  Coefbeta5sim[j] - 1.96*Coefbeta5simp[j])
  
}

mean(COVE1)
mean(COVE2)
mean(COVE3)
mean(COVE4)
mean(COVE5)

