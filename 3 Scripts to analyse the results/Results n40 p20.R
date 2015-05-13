###Calculating preformance of GEE, GEE_backwardselct and PGEE###
library(gee)
library(geepack)

#Obtaining estimates for GEE
EstimatesGEE <- matrix(1,20,100)
for(i in 1:100){
  AAB<-summary(geeglm(outcome ~ V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18+V19+V20+V21, id= id, data=Sim100N40P20uu[[i]], corstr="exchangeable", family = binomial))
  EstimatesGEE[,i] <- AAB$coefficients[2:21,1]
}

#Backward selecting GEE
p<-20
EstimatesB <- matrix(NA,20,100)
EstimatesBack <- matrix(NA,20,100)
row.names(EstimatesBack) <-c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05", "X.tempV06", "X.tempV07",
                             "X.tempV08", "X.tempV09", "X.tempV10", "X.tempV11","X.tempV12", "X.tempV13", 
                             "X.tempV14", "X.tempV15", "X.tempV16", "X.tempV17",
                             "X.tempV18", "X.tempV19", "X.tempV20", "X.tempV21")
row.names(EstimatesB) <- c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05", "X.tempV06", "X.tempV07",
                           "X.tempV08", "X.tempV09", "X.tempV10", "X.tempV11","X.tempV12", "X.tempV13", 
                           "X.tempV14", "X.tempV15", "X.tempV16", "X.tempV17",
                           "X.tempV18", "X.tempV19", "X.tempV20", "X.tempV21")
for(i in 1:100){
  X.temp<<- as.matrix(Sim100N40P20uu[[i]][,2:(p+1)])
  colnames(X.temp)<-c("V02", "V03", "V04", "V05", "V06", "V07",
                      "V08", "V09", "V10", "V11","V12", "V13", "V14", "V15", "V16", "V17",
                      "V18", "V19", "V20", "V21")
  Y <<- as.matrix(Sim100N40P20uu[[i]][,(p+2)])
  glmFull <- geeglm(Y ~ X.temp, id= id, Sim100N40P20uu[[i]], corstr="exchangeable", family = binomial)
  ids <- which(coef(summary(glmFull))[-1, 4] < 0.5)
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N40P20uu[[i]], corstr="exchangeable", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.3)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N40P20uu[[i]], corstr="exchangeable", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.2)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N40P20uu[[i]], corstr="exchangeable", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.1)]
  X.temp <- X.temp[,ids]
  GB<- geeglm(Y ~ X.temp, id= id, Sim100N40P20uu[[i]], corstr="exchangeable", family = binomial)
  GBA <- GB$coefficients
  BBAA <-merge(EstimatesB[,i], GBA, by = "row.names", all.x = T,sort=T)
  EstimatesBack[,i] <- BBAA$y
}
EstimatesBack[is.na(EstimatesBack)]<- 0

#Calculating MSE for GEE
MEGEE <- matrix(1,100)
BetaTrue <- c(0.8,-0.5,-0.3,0.3,-0.3,0.3,-0.3,0,0,0,0,0,0,0,0,0,0,0,0,0)
for(i in 1:100){
  MEGEE[i] <- t(EstimatesGEE[,i]-BetaTrue)%*%(EstimatesGEE[,i]-BetaTrue)
}

#Calculating MSE GEE_Backward selection
MEGEEBACK <- matrix(1,100)
BetaTrue <- c(0.8,-0.5,-0.3,0.3,-0.3,0.3,-0.3,0,0,0,0,0,0,0,0,0,0,0,0,0)
for(i in 1:100){
  MEGEEBACK[i] <- t(EstimatesBack[,i]-BetaTrue)%*%(EstimatesBack[,i]-BetaTrue)
}

#Calculating MSE PGEE
MESCAD <- matrix(1,100)
BetaTrue <- c(0.8,-0.5,-0.3,0.3,-0.3,0.3,-0.3,0,0,0,0,0,0,0,0,0,0,0,0,0)
#BetaTrue <- c(0.5,0.3,-0.3,-0.3,0.3,0.3,-0.3,-0.3,0.3,-0.3,0,0,0,0,0,0,0,0,0,0)
#BetaTrue <- c(0.8,-0.5,-0.3,-0.3,0.3,-0.3,0.3,-0.3,0.3,0,0,0,0,0,0,0,0,0,0,0)
for(i in 1:100){
  MESCAD[i] <- t(EstimatesSCAD4020u[,i]-BetaTrue)%*%(EstimatesSCAD4020u[,i]-BetaTrue)
}

#Mean MSE
mean(MEGEE)
mean(MEGEEBACK)
mean(MESCAD)
sd(MEGEE)/sqrt(100)
sd(MEGEEBACK)/sqrt(100)
sd(MESCAD)/sqrt(100)


#Underselction
UnderselectBACK <- mean(EstimatesBack[1:7,]==0)
UnderselectPGEE <- mean(EstimatesSCAD4020u[1:7,]==0)

#Overselection
OverselectBACK <- 1-mean(EstimatesBack[8:20,]==0)
OverselectPGEE <- 1-mean(EstimatesSCAD4020u[8:20,]==0)
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
EstimatesPGEEWithout <- EstimatesSCAD4020u
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
Coefbeta11sim <-matrix(1,100)
Coefbeta12sim <-matrix(1,100)
Coefbeta13sim <-matrix(1,100)
Coefbeta14sim <-matrix(1,100)
Coefbeta15sim <-matrix(1,100)
Coefbeta16sim <-matrix(1,100)
Coefbeta17sim <-matrix(1,100)
Coefbeta18sim <-matrix(1,100)
Coefbeta19sim <-matrix(1,100)
Coefbeta20sim <-matrix(1,100)


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
Coefbeta11simp<-matrix(1,100)
Coefbeta12simp<-matrix(1,100)
Coefbeta13simp<-matrix(1,100)
Coefbeta14simp<-matrix(1,100)
Coefbeta15simp<-matrix(1,100)
Coefbeta16simp<-matrix(1,100)
Coefbeta17simp<-matrix(1,100)
Coefbeta18simp<-matrix(1,100)
Coefbeta19simp<-matrix(1,100)
Coefbeta20simp<-matrix(1,100)

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
COVE11<-matrix(1,100)
COVE12<-matrix(1,100)
COVE13<-matrix(1,100)
COVE14<-matrix(1,100)
COVE15<-matrix(1,100)
COVE16<-matrix(1,100)
COVE17<-matrix(1,100)
COVE18<-matrix(1,100)
COVE19<-matrix(1,100)
COVE20<-matrix(1,100)

for(j in 1:100){
  
  ####Analysing simulated data####
  GEEX1coef<- summary(geeglm(outcome ~ V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18+V19+V20+V21, id= id, data=Sim100N40P20uu[[j]], corstr="exchangeable", family = binomial))
  
  #Saving coefficients beta0 and beta1 for each replicated dataset
 
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
  Coefbeta11sim[j] <- GEEX1coef$coefficients[12,1]
  Coefbeta12sim[j] <- GEEX1coef$coefficients[13,1]
  Coefbeta13sim[j] <- GEEX1coef$coefficients[14,1]
  Coefbeta14sim[j] <- GEEX1coef$coefficients[15,1]
  Coefbeta15sim[j] <- GEEX1coef$coefficients[16,1]
  Coefbeta16sim[j] <- GEEX1coef$coefficients[17,1]
  Coefbeta17sim[j] <- GEEX1coef$coefficients[18,1]  
  Coefbeta18sim[j] <- GEEX1coef$coefficients[19,1]
  Coefbeta19sim[j] <- GEEX1coef$coefficients[20,1]
  Coefbeta20sim[j] <- GEEX1coef$coefficients[21,1]
  #   
  
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
  Coefbeta11simp[j] <- GEEX1coef$coefficients[12,2]
  Coefbeta12simp[j] <- GEEX1coef$coefficients[13,2]
  Coefbeta13simp[j] <- GEEX1coef$coefficients[14,2]
  Coefbeta14simp[j] <- GEEX1coef$coefficients[15,2]
  Coefbeta15simp[j] <- GEEX1coef$coefficients[16,2]
  Coefbeta16simp[j] <- GEEX1coef$coefficients[17,2]
  Coefbeta17simp[j] <- GEEX1coef$coefficients[18,2]  
  Coefbeta18simp[j] <- GEEX1coef$coefficients[19,2]
  Coefbeta19simp[j] <- GEEX1coef$coefficients[20,2]
  Coefbeta20simp[j] <- GEEX1coef$coefficients[21,2]
  #   
  
  COVE1[j] <- (Coefbeta1sim[j] - 1.96*Coefbeta1simp[j] <= 0.8 & 0.8 >=  Coefbeta1sim[j] - 1.96*Coefbeta1simp[j])
  COVE2[j] <- (Coefbeta2sim[j] - 1.96*Coefbeta2simp[j] <= -0.5 & -0.5 >=  Coefbeta2sim[j] - 1.96*Coefbeta2simp[j])
  COVE3[j] <- (Coefbeta3sim[j] - 1.96*Coefbeta3simp[j] <= -0.3 & -0.3 >=  Coefbeta3sim[j] - 1.96*Coefbeta3simp[j])
  COVE4[j] <- (Coefbeta4sim[j] - 1.96*Coefbeta4simp[j] <= 0.3 & 0.3 >=  Coefbeta4sim[j] - 1.96*Coefbeta4simp[j])
  COVE5[j] <- (Coefbeta5sim[j] - 1.96*Coefbeta5simp[j] <=  -0.3 & -0.3 >=  Coefbeta5sim[j] - 1.96*Coefbeta5simp[j])
  COVE6[j] <- (Coefbeta6sim[j] - 1.96*Coefbeta6simp[j] <= 0.3 & 0.3 >=  Coefbeta6sim[j] - 1.96*Coefbeta6simp[j])
  COVE7[j] <- (Coefbeta7sim[j] - 1.96*Coefbeta7simp[j] <= -0.3 & -0.3 >=  Coefbeta7sim[j] - 1.96*Coefbeta7simp[j])
  COVE8[j] <- (Coefbeta8sim[j] - 1.96*Coefbeta8simp[j] <= 0 & 0 >=  Coefbeta8sim[j] - 1.96*Coefbeta8simp[j])
  COVE9[j] <- (Coefbeta9sim[j] - 1.96*Coefbeta9simp[j] <= 0 & 0 >=  Coefbeta9sim[j] - 1.96*Coefbeta9simp[j])
  COVE10[j] <- (Coefbeta10sim[j] - 1.96*Coefbeta10simp[j] <= 0 & 0 >=  Coefbeta10sim[j] - 1.96*Coefbeta10simp[j])
  COVE11[j] <- (Coefbeta11sim[j] - 1.96*Coefbeta11simp[j] <= 0 & 0 >=  Coefbeta11sim[j] - 1.96*Coefbeta11simp[j])
  COVE12[j] <- (Coefbeta12sim[j] - 1.96*Coefbeta12simp[j] <= 0 & 0 >=  Coefbeta12sim[j] - 1.96*Coefbeta12simp[j])
  COVE13[j] <- (Coefbeta13sim[j] - 1.96*Coefbeta13simp[j] <= 0 & 0 >=  Coefbeta13sim[j] - 1.96*Coefbeta13simp[j])
  COVE14[j] <- (Coefbeta14sim[j] - 1.96*Coefbeta14simp[j] <= 0 & 0 >=  Coefbeta14sim[j] - 1.96*Coefbeta14simp[j])
  COVE15[j] <- (Coefbeta15sim[j] - 1.96*Coefbeta15simp[j] <=  0 & 0 >=  Coefbeta15sim[j] - 1.96*Coefbeta15simp[j])
  COVE16[j] <- (Coefbeta16sim[j] - 1.96*Coefbeta16simp[j] <= 0 & 0 >=  Coefbeta16sim[j] - 1.96*Coefbeta16simp[j])
  COVE17[j] <- (Coefbeta17sim[j] - 1.96*Coefbeta17simp[j] <= 0 & 0 >=  Coefbeta17sim[j] - 1.96*Coefbeta17simp[j])
  COVE18[j] <- (Coefbeta18sim[j] - 1.96*Coefbeta18simp[j] <= 0 & 0 >=  Coefbeta18sim[j] - 1.96*Coefbeta18simp[j])
  COVE19[j] <- (Coefbeta19sim[j] - 1.96*Coefbeta19simp[j] <= 0 & 0 >=  Coefbeta19sim[j] - 1.96*Coefbeta19simp[j])
  COVE20[j] <- (Coefbeta20sim[j] - 1.96*Coefbeta20simp[j] <= 0 & 0 >=  Coefbeta20sim[j] - 1.96*Coefbeta20simp[j])
  
  
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
mean(COVE11)
mean(COVE12)
mean(COVE13)
mean(COVE14)
mean(COVE15)
mean(COVE16)
mean(COVE17)
mean(COVE18)
mean(COVE19)
mean(COVE20)

