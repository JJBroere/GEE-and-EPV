###Calculating preformance of GEE, GEE_backwardselct and PGEE###
library(gee)
library(geepack)

#Obtaining estimates for GEE
EstimatesGEE <- matrix(1,10,100)
for(i in 1:100){
  AAB<-summary(geeglm(outcome ~ V2+V3+V4+V5+V6+V7+V8+V9+V10+V11, id= id, data=Sim100N40P10u[[i]], corstr="exchangeable", family = binomial))
  EstimatesGEE[,i] <- AAB$coefficients[2:11,1]
}

#Backward selecting GEE
p<-10
EstimatesB <- matrix(NA,10,100)
EstimatesBack <- matrix(NA,10,100)
row.names(EstimatesBack) <-c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05", "X.tempV06", "X.tempV07",
                             "X.tempV08", "X.tempV09", "X.tempV10", "X.tempV11")
row.names(EstimatesB) <- c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05", "X.tempV06", "X.tempV07",
                           "X.tempV08", "X.tempV09", "X.tempV10", "X.tempV11")
for(i in 1:100){
  X.temp<<- as.matrix(Sim100N40P10u[[i]][,2:(p+1)])
  colnames(X.temp)<-c("V02", "V03", "V04", "V05", "V06", "V07",
                      "V08", "V09", "V10", "V11")
  Y <<- as.matrix(Sim100N40P10u[[i]][,(p+2)])
  glmFull <- geeglm(Y ~ X.temp, id= id, Sim100N40P10u[[i]], corstr="independence", family = binomial)
  ids <- which(coef(summary(glmFull))[-1, 4] < 0.5)
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N40P10u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.3)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N40P10u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.2)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N40P10u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.1)]
  X.temp <- X.temp[,ids]
  GB<- geeglm(Y ~ X.temp, id= id, Sim100N40P10u[[i]], corstr="independence", family = binomial)
  GBA <- GB$coefficients
  BBAA <-merge(EstimatesB[,i], GBA, by = "row.names", all.x = T,sort=T)
  EstimatesBack[,i] <- BBAA$y
}

EstimatesBack[is.na(EstimatesBack)]<- 0#Setting notselected varables to zere

#Calculating MSE for GEE
MEGEE <- matrix(1,100)
BetaTrue <- c(0.8,0.8,-0.5,-0.5,-0.3,-0.3,0,0,0,0)
for(i in 1:100){
  MEGEE[i] <- t(EstimatesGEE[,i]-BetaTrue)%*%(EstimatesGEE[,i]-BetaTrue)
}

#Calculating MSE GEE_Backward selection
MEGEEBACK <- matrix(1,100)
BetaTrue <- c(0.8,0.8,-0.5,-0.5,-0.3,-0.3,0,0,0,0)
for(i in 1:100){
  MEGEEBACK[i] <- t(EstimatesBack[,i]-BetaTrue)%*%(EstimatesBack[,i]-BetaTrue)
}

#Calculating MSE PGEE
MESCAD <- matrix(1,100)
BetaTrue <- c(0.8,0.8,-0.5,-0.5,-0.3,-0.3,0,0,0,0)
for(i in 1:100){
  MESCAD[i] <- t(EstimatesSCAD4010u[,i]-BetaTrue)%*%(EstimatesSCAD4010u[,i]-BetaTrue)
}
 
#Mean MSE
mean(MEGEE)
mean(MEGEEBACK)
mean(MESCAD)
sd(MEGEE)/sqrt(100)
sd(MEGEEBACK)/sqrt(100)
sd(MESCAD)/sqrt(100)

#Underselction
UnderselectBACK <- mean(EstimatesBack[1:6,]==0)
UnderselectPGEE <- mean(EstimatesSCAD4010u[1:6,]==0)
#Overselection
OverselectBACK <- 1-mean(EstimatesBack[7:10,]==0)
OverselectPGEE <- 1-mean(EstimatesSCAD4010u[7:10,]==0)
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
EstimatesPGEEWithout <- EstimatesSCAD4010u
EstimatesPGEEWithout[EstimatesPGEEWithout==0] <- NA
RM <-rowMeans(EstimatesPGEEWithout, na.rm =T)
Bias <- (RM-BetaTrue)/BetaTrue
Bias


###Detemine the covarage of GEE###
Coefbeta1sim <-matrix(1,100)
Coefbeta2sim <-matrix(1,100)
Coefbeta3sim <-matrix(1,100)
Coefbeta4sim <-matrix(1,100)
Coefbeta5sim <-matrix(1,100)
Coefbeta6sim <-matrix(1,100)
Coefbeta7sim <-matrix(1,100)
Coefbeta8sim <-matrix(1,100)
Coefbeta9sim <-matrix(1,100)
Coefbeta10sim <-matrix(1,100)


Coefbeta1simp<-matrix(1,100)
Coefbeta2simp<-matrix(1,100)
Coefbeta3simp<-matrix(1,100)
Coefbeta4simp<-matrix(1,100)
Coefbeta5simp<-matrix(1,100)
Coefbeta6simp<-matrix(1,100)
Coefbeta7simp<-matrix(1,100)
Coefbeta8simp<-matrix(1,100)
Coefbeta9simp<-matrix(1,100)
Coefbeta10simp<-matrix(1,100)

COVE1<-matrix(1,100)
COVE2<-matrix(1,100)
COVE3<-matrix(1,100)
COVE4<-matrix(1,100)
COVE5<-matrix(1,100)
COVE6<-matrix(1,100)
COVE7<-matrix(1,100)
COVE8<-matrix(1,100)
COVE9<-matrix(1,100)
COVE10<-matrix(1,100)

for(j in 1:100){
  
  ####Analysing simulated data####
  GEEX1coef<- summary(geeglm(outcome ~ V2+V3+V4+V5+V6+V7+V8+V9+V10+V11, id= id, data=Sim100N40P10u[[j]], corstr="exchangeable", family = binomial))
  
  #Saving coefficients beta0 and beta1 for each replicated dataset
     Coefbeta0sim[j] <- GEEX1coef$coefficients[1,1]
     Coefbeta1sim[j] <- GEEX1coef$coefficients[2,1]
     Coefbeta2sim[j] <- GEEX1coef$coefficients[3,1]
     Coefbeta3sim[j] <- GEEX1coef$coefficients[4,1]
     Coefbeta4sim[j] <- GEEX1coef$coefficients[5,1]
     Coefbeta5sim[j] <- GEEX1coef$coefficients[6,1]
     Coefbeta6sim[j] <- GEEX1coef$coefficients[7,1]
     Coefbeta7sim[j] <- GEEX1coef$coefficients[8,1]  
     Coefbeta8sim[j] <- GEEX1coef$coefficients[9,1]
     Coefbeta9sim[j] <- GEEX1coef$coefficients[10,1]
    Coefbeta10sim[j] <- GEEX1coef$coefficients[11,1]
  #   
     Coefbeta0simp[j] <- GEEX1coef$coefficients[1,2]
     Coefbeta1simp[j] <- GEEX1coef$coefficients[2,2]
     Coefbeta2simp[j] <- GEEX1coef$coefficients[3,2]
     Coefbeta3simp[j] <- GEEX1coef$coefficients[4,2]
     Coefbeta4simp[j] <- GEEX1coef$coefficients[5,2]
     Coefbeta5simp[j] <- GEEX1coef$coefficients[6,2]
     Coefbeta6simp[j] <- GEEX1coef$coefficients[7,2]
     Coefbeta7simp[j] <- GEEX1coef$coefficients[8,2]  
     Coefbeta8simp[j] <- GEEX1coef$coefficients[9,2]
  Coefbeta9simp[j] <- GEEX1coef$coefficients[10,2]
  Coefbeta10simp[j] <- GEEX1coef$coefficients[11,2]
  #   
   
   COVE1[j] <- (Coefbeta1sim[j] - 1.96*Coefbeta1simp[j] <= 0.8 & 0.8 >=  Coefbeta1sim[j] - 1.96*Coefbeta1simp[j])
   COVE2[j] <- (Coefbeta2sim[j] - 1.96*Coefbeta2simp[j] <= 0.8 & 0.8 >=  Coefbeta2sim[j] - 1.96*Coefbeta2simp[j])
   COVE3[j] <- (Coefbeta3sim[j] - 1.96*Coefbeta3simp[j] <= -0.5 & -0.5 >=  Coefbeta3sim[j] - 1.96*Coefbeta3simp[j])
   COVE4[j] <- (Coefbeta4sim[j] - 1.96*Coefbeta4simp[j] <= -0.5 & -0.5 >=  Coefbeta4sim[j] - 1.96*Coefbeta4simp[j])
   COVE5[j] <- (Coefbeta5sim[j] - 1.96*Coefbeta5simp[j] <=  -0.3 & -0.3 >=  Coefbeta5sim[j] - 1.96*Coefbeta5simp[j])
   COVE6[j] <- (Coefbeta6sim[j] - 1.96*Coefbeta6simp[j] <= -0.3 & -0.3 >=  Coefbeta6sim[j] - 1.96*Coefbeta6simp[j])
 COVE7[j] <- (Coefbeta7sim[j] - 1.96*Coefbeta7simp[j] <= 0 & 0 >=  Coefbeta7sim[j] - 1.96*Coefbeta7simp[j])
   COVE8[j] <- (Coefbeta8sim[j] - 1.96*Coefbeta8simp[j] <= 0 & 0 >=  Coefbeta8sim[j] - 1.96*Coefbeta8simp[j])
  COVE9[j] <- (Coefbeta9sim[j] - 1.96*Coefbeta9simp[j] <= 0 & 0 >=  Coefbeta9sim[j] - 1.96*Coefbeta9simp[j])
 COVE10[j] <- (Coefbeta10sim[j] - 1.96*Coefbeta10simp[j] <= 0 & 0 >=  Coefbeta10sim[j] - 1.96*Coefbeta10simp[j])
  
  
}


 mean(COVE1)
 mean(COVE2)
 mean(COVE3)
 mean(COVE4)
 mean(COVE5)
 mean(COVE6)
 mean(COVE7)
 mean(COVE8)
mean(COVE9)
mean(COVE10)




