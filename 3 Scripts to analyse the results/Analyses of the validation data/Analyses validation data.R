###Script for analysing the validation datasets###
###In order to run the full script, the simulated data files and the PGEE estimates need 
#to be loaded###
library(geepack)

###Prediction validation N40 P20###

#Calculating GEE coefficients
EstimatesGEE <- matrix(1,20,100)
for(i in 1:100){
  AAB<-summary(geeglm(outcome ~ V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18+V19+V20+V21, id= id, data=Sim100N40P20uu[[i]], corstr="exchangeable", family = binomial))
  EstimatesGEE[,i] <- AAB$coefficients[2:21,1]
}

#Calculting model error for GEE
CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp20uu[3:22])%*% as.vector(EstimatesGEE[,i])
  PP<- exp(PA)/(1+exp(PA))

ClassCheck<-cbind(as.matrix(Valp20uu[,24]),PP)

CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])

}
mean(CheckCL)
sd(CheckCL)

#Calculating coeffiecient of GEE with backward selection
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

#Calculting model error for GEE with backward selection
CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp20uu[3:22])%*% as.vector(EstimatesBack[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp20uu[,24]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)

#Calculting model error for penalized GEE
CheckCL2 <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp20uu[3:22])%*% as.vector(EstimatesSCAD4020u[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp20uu[,24]),PP)
  
  CheckCL2[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)


####Prediction validation N40 P10, using the same steps as above
EstimatesGEE <- matrix(1,10,100)
for(i in 1:100){
  AAB<-summary(geeglm(outcome ~ V2+V3+V4+V5+V6+V7+V8+V9+V10+V11, id= id, data=Sim100N40P10u[[i]], corstr="exchangeable", family = binomial))
  EstimatesGEE[,i] <- AAB$coefficients[2:11,1]
}
CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp10uu[3:12])%*% as.vector(EstimatesGEE[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp10uu[,14]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)

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

EstimatesBack[is.na(EstimatesBack)]<- 0

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp10uu[3:12])%*% as.vector(EstimatesBack[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp10uu[,14]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)


CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp10uu[3:12])%*% as.vector(EstimatesSCAD4010u[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp10uu[,14]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)


####Validation N40 P5
EstimatesGEE <- matrix(1,5,100)
for(i in 1:100){
  AAB<-summary(geeglm(outcome ~ V2+V3+V4+V5+V6, id= id, data=Sim100N40P5u[[i]], corstr="exchangeable", family = binomial))
  EstimatesGEE[,i] <- AAB$coefficients[2:6,1]
}

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp5uu[3:7])%*% as.vector(EstimatesGEE[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp5uu[,9]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)
#!!!Does not select anything in dataset 34. Run from 1:33 and 35:100, other wise procedure stops!!!
p<-5
EstimatesB <- matrix(NA,5,100)
EstimatesBack <- matrix(NA,5,100)
row.names(EstimatesBack) <-c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05", "X.tempV06")
row.names(EstimatesB) <- c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05", "X.tempV06")
for(i in 35:100){
  X.temp<<- as.matrix(Sim100N40P5u[[i]][,2:(p+1)])
  colnames(X.temp)<-c("V02", "V03", "V04", "V05", "V06")
  Y <<- as.matrix(Sim100N40P5u[[i]][,(p+2)])
  glmFull <- geeglm(Y ~ X.temp, id= id, Sim100N40P5u[[i]], corstr="exchangeable", family = binomial)
  ids <- which(coef(summary(glmFull))[-1, 4] < 0.5)
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N40P5u[[i]], corstr="exchangeable", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.3)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N40P5u[[i]], corstr="exchangeable", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.2)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N40P5u[[i]], corstr="exchangeable", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.1)]
  X.temp <- X.temp[,ids]
  GB<- geeglm(Y ~ X.temp, id= id, Sim100N40P5u[[i]], corstr="exchangeable", family = binomial)
  GBA <- GB$coefficients
  BBAA <-merge(EstimatesB[,i], GBA, by = "row.names", all.x = T,sort=T)
  EstimatesBack[,i] <- BBAA$y
}
EstimatesBack[is.na(EstimatesBack)]<- 0

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp5uu[3:7])%*% as.vector(EstimatesBack[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp5uu[,9]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)


CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp5uu[3:7])%*% as.vector(EstimatesSCAD4005[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp5uu[,9]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)

####Prediction validation N40 P4
EstimatesGEE <- matrix(1,4,100)
for(i in 1:100){
  AAB<-summary(geeglm(outcome ~ V2+V3+V4+V5, id= id, data=Sim100N40P4u[[i]], corstr="exchangeable", family = binomial))
  EstimatesGEE[,i] <- AAB$coefficients[2:5,1]
}

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp4uu[3:6])%*% as.vector(EstimatesGEE[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp4uu[,8]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)

p<-4
EstimatesB <- matrix(NA,4,100)
EstimatesBack <- matrix(NA,4,100)
row.names(EstimatesBack) <-c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05")
row.names(EstimatesB) <- c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05")
for(i in 1:100){
  X.temp<<- as.matrix(Sim100N40P4u[[i]][,2:(p+1)])
  colnames(X.temp)<-c("V02", "V03", "V04", "V05")
  Y <<- as.matrix(Sim100N40P4u[[i]][,(p+2)])
  glmFull <- geeglm(Y ~ X.temp, id= id, Sim100N40P4u[[i]], corstr="exchangeable", family = binomial)
  ids <- which(coef(summary(glmFull))[-1, 4] < 0.5)
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N40P4u[[i]], corstr="exchangeable", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.3)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N40P4u[[i]], corstr="exchangeable", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.2)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N40P4u[[i]], corstr="exchangeable", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.1)]
  X.temp <- X.temp[,ids]
  GB<- geeglm(Y ~ X.temp, id= id, Sim100N40P4u[[i]], corstr="exchangeable", family = binomial)
  GBA <- GB$coefficients
  BBAA <-merge(EstimatesB[,i], GBA, by = "row.names", all.x = T,sort=T)
  EstimatesBack[,i] <- BBAA$y
}
EstimatesBack[is.na(EstimatesBack)]<- 0

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp4uu[3:6])%*% as.vector(EstimatesBack[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp4uu[,8]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp4uu[3:6])%*% as.vector(EstimatesSCAD4004[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp4uu[,8]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)

####Prediction validation N40 P2
EstimatesGEE <- matrix(1,2,100)
AA<-c(1,2)
for(i in 1:100){
  AAB<-summary(geeglm(outcome ~ V2+V3, id= id, data=Sim100N40P2u[[i]], corstr="exchangeable", family = binomial))
  EstimatesGEE[,i] <- AAB$coefficients[2:3,1]
}

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp2uu[3:4])%*% as.vector(EstimatesGEE[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp2uu[,6]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)

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

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp2uu[3:4])%*% as.vector(EstimatesBACK[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp2uu[,6]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)


CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp2uu[3:4])%*% as.vector(EstimatesSCAD4002[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp2uu[,6]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)

####Prediction validation N100 P50
EstimatesGEE <- matrix(1,50,100)
for(i in 1:100){
  AAB<-summary(geeglm(outcome ~ V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18+V19+V20+V21+V22+V23+V24+V25+V26+
                        V27+V28+V29+V30+V31+V32+V33+V34+V35+V36+V37+V38+V39+V40+V41+V42+V43+V44+V45+V46+V47+V48+V49+V50+V51, id= id, data=Sim100N100P50u[[i]], corstr="exchangeable", family = binomial))
  EstimatesGEE[,i] <- AAB$coefficients[2:51,1]
}
CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp50uu[3:52])%*% as.vector(EstimatesGEE[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp50uu[,54]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)


#Backward selecting GEE
p<-50
EstimatesB <- matrix(NA,50,100)
EstimatesBack <- matrix(NA,50,100)
row.names(EstimatesBack) <-c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05", "X.tempV06", "X.tempV07",
                             "X.tempV08", "X.tempV09", "X.tempV10", "X.tempV11",
                             "X.tempV12", "X.tempV13", "X.tempV14", "X.tempV15", "X.tempV16", "X.tempV17",
                             "X.tempV18", "X.tempV19", "X.tempV20", "X.tempV21","X.tempV22", "X.tempV23", 
                             "X.tempV24", "X.tempV25", "X.tempV26","X.tempV27",
                             "X.tempV28", "X.tempV29", "X.tempV30", "X.tempV31","X.tempV32", "X.tempV33", 
                             "X.tempV34", "X.tempV35", "X.tempV36","X.tempV37",
                             "X.tempV38", "X.tempV39", "X.tempV40", "X.tempV41","X.tempV42", "X.tempV43", 
                             "X.tempV44", "X.tempV45", "X.tempV46", "X.tempV47",
                             "X.tempV48", "X.tempV49", "X.tempV50","X.tempV51")
row.names(EstimatesB) <- c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05", "X.tempV06", "X.tempV07",
                           "X.tempV08", "X.tempV09", "X.tempV10", "X.tempV11",
                           "X.tempV12", "X.tempV13", "X.tempV14", "X.tempV15", "X.tempV16", "X.tempV17",
                           "X.tempV18", "X.tempV19", "X.tempV20", "X.tempV21","X.tempV22", "X.tempV23", 
                           "X.tempV24", "X.tempV25", "X.tempV26","X.tempV27",
                           "X.tempV28", "X.tempV29", "X.tempV30", "X.tempV31","X.tempV32", "X.tempV33", 
                           "X.tempV34", "X.tempV35", "X.tempV36","X.tempV37",
                           "X.tempV38", "X.tempV39", "X.tempV40", "X.tempV41","X.tempV42", "X.tempV43", 
                           "X.tempV44", "X.tempV45", "X.tempV46", "X.tempV47",
                           "X.tempV48", "X.tempV49", "X.tempV50","X.tempV51")
for(i in 1:100){
  X.temp<<- as.matrix(Sim100N100P50u[[i]][,2:(p+1)])
  colnames(X.temp)<-c("V02", "V03", "V04", "V05", "V06", "V07",
                      "V08", "V09", "V10", "V11","V12", "V13", "V14", "V15", "V16", "V17",
                      "V18", "V19", "V20", "V21","V22", "V23", "V24", "V25", "V26", "V27",
                      "V28", "V29", "V30", "V31","V32", "V33", "V34", "V35", "V36","V37",
                      "V38", "V39", "V40", "V41","V42", "V43", "V44", "V45", "V46", "V47",
                      "V48", "V49", "V50","V51")
  Y <<- as.matrix(Sim100N100P50u[[i]][,(p+2)])
  glmFull <- geeglm(Y ~ X.temp, id= id, Sim100N100P50u[[i]], corstr="independence", family = binomial)
  ids <- which(coef(summary(glmFull))[-1, 4] < 0.5)
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N100P50u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.3)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N100P50u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.2)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N100P50u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.1)]
  X.temp <- X.temp[,ids]
  GB<- geeglm(Y ~ X.temp, id= id, Sim100N100P50u[[i]], corstr="independence", family = binomial)
  GBA <- GB$coefficients
  BBAA <-merge(EstimatesB[,i], GBA, by = "row.names", all.x = T,sort=T)
  EstimatesBack[,i] <- BBAA$y
}

EstimatesBack[is.na(EstimatesBack)]<- 0#Setting notselected varables to zere

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp50uu[3:52])%*% as.vector(EstimatesBack[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp50uu[,54]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)


CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp50uu[3:52])%*% as.vector(EstimatesSCAD10050u[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp50uu[,54]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)


####Prediction validation N100 P25
EstimatesGEE <- matrix(1,25,100)
for(i in 1:100){
  AAB<-summary(geeglm(outcome ~ V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18+V19+V20+V21+V22+V23+V24+V25+V26, id= id, data=Sim100N100P25u[[i]], corstr="exchangeable", family = binomial))
  EstimatesGEE[,i] <- AAB$coefficients[2:26,1]
}
CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp25uu[3:27])%*% as.vector(EstimatesGEE[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp25uu[,29]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)

#Backward selecting GEE
p<-25
EstimatesB <- matrix(NA,25,100)
EstimatesBack <- matrix(NA,25,100)
row.names(EstimatesBack) <-c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05", "X.tempV06", "X.tempV07",
                             "X.tempV08", "X.tempV09", "X.tempV10", "X.tempV11",
                             "X.tempV12", "X.tempV13", "X.tempV14", "X.tempV15", "X.tempV16", "X.tempV17",
                             "X.tempV18", "X.tempV19", "X.tempV20", "X.tempV21","X.tempV22", "X.tempV23", 
                             "X.tempV24", "X.tempV25", "X.tempV26")
row.names(EstimatesB) <- c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05", "X.tempV06", "X.tempV07",
                           "X.tempV08", "X.tempV09", "X.tempV10", "X.tempV11",
                           "X.tempV12", "X.tempV13", "X.tempV14", "X.tempV15", "X.tempV16", "X.tempV17",
                           "X.tempV18", "X.tempV19", "X.tempV20", "X.tempV21","X.tempV22", "X.tempV23", 
                           "X.tempV24", "X.tempV25", "X.tempV26")
for(i in 1:100){
  X.temp<<- as.matrix(Sim100N100P25u[[i]][,2:(p+1)])
  colnames(X.temp)<-c("V02", "V03", "V04", "V05", "V06", "V07",
                      "V08", "V09", "V10", "V11","V12", "V13", "V14", "V15", "V16", "V17",
                      "V18", "V19", "V20", "V21","V22", "V23", "V24", "V25", "V26")
  Y <<- as.matrix(Sim100N100P25u[[i]][,(p+2)])
  glmFull <- geeglm(Y ~ X.temp, id= id, Sim100N100P25u[[i]], corstr="independence", family = binomial)
  ids <- which(coef(summary(glmFull))[-1, 4] < 0.5)
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N100P25u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.3)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N100P25u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.2)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N100P25u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.1)]
  X.temp <- X.temp[,ids]
  GB<- geeglm(Y ~ X.temp, id= id, Sim100N100P25u[[i]], corstr="independence", family = binomial)
  GBA <- GB$coefficients
  BBAA <-merge(EstimatesB[,i], GBA, by = "row.names", all.x = T,sort=T)
  EstimatesBack[,i] <- BBAA$y
}

EstimatesBack[is.na(EstimatesBack)]<- 0#Setting notselected varables to zere

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp25uu[3:27])%*% as.vector(EstimatesBack[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp25uu[,29]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)




CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp25uu[3:27])%*% as.vector(EstimatesSCAD10025u[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp25uu[,29]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)


####Prediction validation N100 P13
EstimatesGEE <- matrix(1,13,100)
for(i in 1:100){
  AAB<-summary(geeglm(outcome ~ V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14, id= id, data=Sim100N100P13u[[i]], corstr="exchangeable", family = binomial))
  EstimatesGEE[,i] <- AAB$coefficients[2:14,1]
}

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp13uu[3:15])%*% as.vector(EstimatesGEE[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp13uu[,17]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)


#Backward selecting GEE
p<-13
EstimatesB <- matrix(NA,13,100)
EstimatesBack <- matrix(NA,13,100)
row.names(EstimatesBack) <-c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05", "X.tempV06", "X.tempV07",
                             "X.tempV08", "X.tempV09", "X.tempV10", "X.tempV11",
                             "X.tempV12", "X.tempV13", "X.tempV14")
row.names(EstimatesB) <- c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05", "X.tempV06", "X.tempV07",
                           "X.tempV08", "X.tempV09", "X.tempV10", "X.tempV11",
                           "X.tempV12", "X.tempV13", "X.tempV14")
for(i in 1:100){
  X.temp<<- as.matrix(Sim100N100P13u[[i]][,2:(p+1)])
  colnames(X.temp)<-c("V02", "V03", "V04", "V05", "V06", "V07",
                      "V08", "V09", "V10", "V11","V12", "V13", "V14")
  Y <<- as.matrix(Sim100N100P13u[[i]][,(p+2)])
  glmFull <- geeglm(Y ~ X.temp, id= id, Sim100N100P13u[[i]], corstr="independence", family = binomial)
  ids <- which(coef(summary(glmFull))[-1, 4] < 0.5)
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N100P13u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.3)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N100P13u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.2)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N100P13u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.1)]
  X.temp <- X.temp[,ids]
  GB<- geeglm(Y ~ X.temp, id= id, Sim100N100P13u[[i]], corstr="independence", family = binomial)
  GBA <- GB$coefficients
  BBAA <-merge(EstimatesB[,i], GBA, by = "row.names", all.x = T,sort=T)
  EstimatesBack[,i] <- BBAA$y
}

EstimatesBack[is.na(EstimatesBack)]<- 0#Setting notselected varables to zere

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp13uu[3:15])%*% as.vector(EstimatesBack[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp13uu[,17]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp13uu[3:15])%*% as.vector(EstimatesSCAD10013u[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp13uu[,17]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)

####Prediction validation N100 P10
EstimatesGEE <- matrix(1,10,100)
for(i in 1:100){
  AAB<-summary(geeglm(outcome ~ V2+V3+V4+V5+V6+V7+V8+V9+V10+V11, id= id, data=Sim100N100P10u[[i]], corstr="exchangeable", family = binomial))
  EstimatesGEE[,i] <- AAB$coefficients[2:11,1]
}
CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp10uu[3:12])%*% as.vector(EstimatesGEE[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp10uu[,14]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)

#Backward selecting GEE
p<-10
EstimatesB <- matrix(NA,10,100)
EstimatesBack <- matrix(NA,10,100)
row.names(EstimatesBack) <-c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05", "X.tempV06", "X.tempV07",
                             "X.tempV08", "X.tempV09", "X.tempV10", "X.tempV11")
row.names(EstimatesB) <- c("X.tempV02", "X.tempV03", "X.tempV04", "X.tempV05", "X.tempV06", "X.tempV07",
                           "X.tempV08", "X.tempV09", "X.tempV10", "X.tempV11")
for(i in 1:100){
  X.temp<<- as.matrix(Sim100N100P10u[[i]][,2:(p+1)])
  colnames(X.temp)<-c("V02", "V03", "V04", "V05", "V06", "V07",
                      "V08", "V09", "V10", "V11")
  Y <<- as.matrix(Sim100N100P10u[[i]][,(p+2)])
  glmFull <- geeglm(Y ~ X.temp, id= id, Sim100N100P10u[[i]], corstr="independence", family = binomial)
  ids <- which(coef(summary(glmFull))[-1, 4] < 0.5)
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N100P10u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.3)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N100P10u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.2)]
  glmFull <- geeglm(Y ~ X.temp[,ids], id= id, Sim100N100P10u[[i]], corstr="independence", family = binomial)
  ids <- ids[which(coef(summary(glmFull))[-1, 4] < 0.1)]
  X.temp <- X.temp[,ids]
  GB<- geeglm(Y ~ X.temp, id= id, Sim100N100P10u[[i]], corstr="independence", family = binomial)
  GBA <- GB$coefficients
  BBAA <-merge(EstimatesB[,i], GBA, by = "row.names", all.x = T,sort=T)
  EstimatesBack[,i] <- BBAA$y
}

EstimatesBack[is.na(EstimatesBack)]<- 0#Setting notselected varables to zere

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp10uu[3:12])%*% as.vector(EstimatesBack[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp10uu[,14]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp10uu[3:12])%*% as.vector(EstimatesSCAD10010u[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp10uu[,14]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)

####Prediction validation N100 P5
EstimatesGEE <- matrix(1,5,100)
for(i in 1:100){
  AAB<-summary(geeglm(outcome ~ V2+V3+V4+V5+V6, id= id, data=Sim100N100P5u[[i]], corstr="exchangeable", family = binomial))
  EstimatesGEE[,i] <- AAB$coefficients[2:6,1]
}

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp5uu[3:7])%*% as.vector(EstimatesGEE[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp5uu[,9]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)

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

CheckCL <- matrix(1,100,10000)
for(i in 1:100){
  PA <-as.matrix(Valp5uu[3:7])%*% as.vector(EstimatesBack[,i])
  PP<- exp(PA)/(1+exp(PA))
  
  ClassCheck<-cbind(as.matrix(Valp5uu[,9]),PP)
  
  CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
  
}
mean(CheckCL)
sd(CheckCL)

  CheckCL <- matrix(1,100,10000)
  for(i in 1:100){
    PA <-as.matrix(Valp5uu[3:7])%*% as.vector(EstimatesSCAD1005u[,i])
    PP<- exp(PA)/(1+exp(PA))
    
    ClassCheck<-cbind(as.matrix(Valp5uu[,9]),PP)
    
    CheckCL[i,]<-abs(ClassCheck[,1]-ClassCheck[,2])
    
  }
  mean(CheckCL)
  sd(CheckCL)



