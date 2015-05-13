########################################################################
########################################################################
########## ALGORITME FAN en LI implementatie voor elastic net penalty
########################################################################
########################################################################




#objects needed 
#data #dataframe met [id,X,Y]
#method <- 'EN', 'BRIDGE' or 'SCAD'
#alpha #parameter 1 = lasso/SCAD, 0 =ridge
#lambda
#V
#M_it 
#tresh_beta <-  #for stability if a betaj closer than this treshold, set to zero
#M_diff <- (convergence indicator= precision)
#k (parameter van initiële RIDGE/DECORRELATIE FIT, put to zero or a small value if it doesnt work increase)
#option <- 'adaptive'  to get adaptive LASSO/EN  and anything else to get nonadaptive version which is recommended
#distr <- "normal"/"binomial"



PEE_LI2 <- function(data,method)
{
  
  #initial beta = OLS-fit
  X_mat <<- as.matrix(data[,2:(p+1)])
  Y_mat <<- as.matrix(data[,(p+2)])
  ID_vec <<- as.vector(data[,1])
  n_obs <<- length(Y_mat)
  n_subj <<- max(data[,1])   #AANPASSING 13/04
  
  beta_decor_init <<- solve(t(X_mat)%*%X_mat + k*diag(p)) %*%(t(X_mat)%*%Y_mat) *(1 + k/length(X_mat[,1]))#aanpassing: ridge initialisation
  #initial gaussian case is ridge
  if(distr=='normal')
  {
    beta <- beta_decor_init #initialize beta
  }
  #initial binomial is glm fit
  if(distr=='binomial')
  {
    fit_glm <- glm(Y_mat~X_mat, family=binomial(link="logit")) 
    beta <- fit_glm$coefficients[2:(p+1)]
  }
  
  
  beta_init <- beta
  ad_weights <<- AD_WEIGHTS(beta_init)
  #
  
  #set to small beta vaules equal to zero
  beta[abs(beta)<tresh_beta] = 0
  
  #iteration of the fitting algorithm
  it <- 1 
  diff_beta <- M_diff + 10 ;# initialiseren anders geen iteraties
  while((it < M_it) & (diff_beta > M_diff))
  { 
    it <- it + 1
    #calculate objects
    est_eq <- EST_EQ_LI2(beta)
    S <- est_eq$S
    H_mat <- est_eq$H_mat
    P_mat <- P_MAT(beta,method)
    U_mat <- P_mat%*%beta     #Hier geen problemen want diagonaal matrix
    
    #Hessiaan berekeningen voor newton raphson: beta waarden nul uithalen, achteraf nul rijen en kollommen toevoegen	
    Hessian_NR <- (H_mat - P_mat)
    Hessian_NR_red <-  Hessian_NR[abs(beta) > tresh_beta,abs(beta) > tresh_beta]
    inv_H_red <- matrix()
    if (sum(abs(beta)) != 0)
    {
      inv_H_red <- solve(Hessian_NR_red)
    }
    
    
    inv_H <- matrix(0,nrow=p,ncol=p)
    inv_H[abs(beta) > tresh_beta,abs(beta) > tresh_beta]  <- inv_H_red 
    
    
    
    
    #beta updata
    beta1 <- beta - inv_H%*%(S - U_mat)    #aanpassing 17/03/2012 geen n_obs meer ervoor 
    #set to small beta values equal to zero
    beta1[abs(beta1) < tresh_beta] = 0
    
    diff_beta <- t(beta-beta1)%*%(beta-beta1)
    beta <- beta1
    
    #print('it')
    #print(it)
    #print('beta')
    #print(beta)
    #print('diff_beta')
    #print(diff_beta)
    #print('#########################################')
  }
  #	t4 <- proc.time();
  
  #convergence?
  if(diff_beta < M_diff) {print('converged')} else {print('no convergence')}
  
  #output/from naÃ”ve elastic net to elastic net
  
  if(method=='EN')
  {
    beta_EN <- beta*(1 + lambda*(1-alpha))
    output <- list(beta,beta_EN,beta_init,diff_beta) 
    names(output) <- list('beta naive elnet','beta elastic net','beta init','diff_beta')
    #print("method is Elastic net regularisation")                                       
  }
  
  if(method =='BRIDGE')  
  {
    output <- list(beta,beta_init)
    names(output) <- list('beta Bridge','beta_init')
    #print("method is Bridge regularisation") 
  }                                                    
  
  if(method =='SCAD')  
  {
    beta_SCAD <- beta*(1 + lambda*(1-alpha))
    output <- list(beta,beta_SCAD, beta_init,diff_beta)
    names(output) <- list('beta naive SCAD', 'beta SCAD','beta_init','diff_beta')
    #print("method is SCAD regularisation") 
  }
  
  if((method != 'EN') & (method != 'BRIDGE') & (method != 'SCAD'))
  {
    output = 'ERROR METHOD FILL IN METHOD="BIDGE" or "EN" or "SCAD"'
  }
  output  
}


EST_EQ_LI2 <- function(beta)    
{
  if(distr == 'normal')
  {S <- t(X_mat)%*%(Y_mat-X_mat%*%beta)
   H <- -t(X_mat)%*%X_mat;
  }
  
  if(distr == 'binomial')
  { S <- as.vector(c(0,0))
    H <- matrix(0,p,p)
    #for(i in 1:n_subj)
    #{
    X_mat_i <- X_mat # X_mat[ID_vec#==i,]   change if you want to do a loop with genearl working corrleation
    Y_mat_i  <-  Y_mat #Y_mat[ID_vec==i,]
    n_i = length(Y_mat_i)
    lin_pred_i <- X_mat_i%*%beta
    #avoid to large a value for the linear predictorl
    lin_pred_i <- apply(lin_pred_i,1,min_100)
    mu_i <- exp(lin_pred_i)/(1+exp(lin_pred_i))
    #numeriek probleem: infinity of exp function --> exact 1 invullen in predictie
    
    
    
    #U_i <- as.vector(mu_i*(1-mu_i))
    #D_i <- diag(mu_i*(1-mu_i))%*%(X_mat_i)
    #V_i  <- diag(U_i^(1/2))%*%diag(n_i)%*%diag(U_i^(1/2)) #Working independence is used, other working correlation not considered;
    #S_i <- t(D_i)%*%solve(V_i) %*%(Y_mat_i - mu_i)
    S_i <- (1/n_obs)*t(X_mat_i)%*%(Y_mat_i - mu_i)
    S <- S + S_i 
    #H_i <- t(D_i)%*%solve(V_i)%*%(D_i)
    H_i <- - (1/n_obs)*t(X_mat_i)%*%diag(as.vector(mu_i*(1-mu_i)))%*%(X_mat_i)   
    H <- H + H_i
    #}
  }
  S[abs(beta) < tresh_beta]=0	
  S  <- S 
  H_mat <- H
  output_est_eq <- list(S,H_mat)
  names(output_est_eq) <- c("S", "H_mat")
  output_est_eq
}


################################################################################

min_100 <- function(x)
{
  output <- min(x,100)
  output
}

################################################################################

##########
###Penalty functions
#########


#approximation elastic net penalty
P_MAT_EN <- function(beta)
{
  p_mat <- matrix(0,nrow=p,ncol=p)
  for(j in 1:p)
  {
    if (abs(beta[j]) < tresh_beta) {p_mat[j,j] = 0}
    else {p_mat[j,j] = lambda*(alpha*(ad_weights[j])  + 2*(1-alpha)*abs(beta[j]))/abs(beta[j])} 
  }
  p_mat
}

#Weights function adaptive EN
AD_WEIGHTS <- function(beta_init)
{
  weights <- 1/beta_init
  weights <- abs(weights)/mean(abs(weights))
  weights*(option == 'adaptive')  + 1*(option != 'adaptive')
}




#approximation bridge penalty 
P_MAT_BRIDGE <- function(beta)
{
  p_mat <- matrix(0,nrow=p,ncol=p)
  for(j in 1:p) 
  {
    if (abs(beta[j]) < tresh_beta) {p_mat[j,j] = 0}
    else {p_mat[j,j] = (lambda*gamma*(abs(beta[j]))^(gamma-1))/abs(beta[j])} 
  }
  p_mat
}

#approximation SCAD penalty
P_MAT_SCAD <- function(beta)
{
  p_mat <- matrix(0,nrow=p,ncol=p)
  for(j in 1:p)
  {
    lambda_SCAD <- lambda*alpha*(ad_weights[j]) #adaptive SCAD part
    A <- (abs(beta[j]) <= lambda_SCAD)*1 ;#aanpassing alpha*
    B <- (a*lambda_SCAD-abs(beta[j])) *  ((a*lambda_SCAD-abs(beta[j])) > 0) #aanpassing lambda_SCAD
    
    
    if (abs(beta[j]) < tresh_beta) {p_mat[j,j] = 0}
    else {p_mat[j,j] = ((lambda_SCAD*(A + B / ((a-1)*lambda_SCAD) * (1-A))     
                         + 2*lambda*(1-alpha)*abs(beta[j]))     )/abs(beta[j])} 
  }
  p_mat
}


#choose method
P_MAT <- function(beta,method='EN')
{
  if (method == 'EN') {p_mat = P_MAT_EN(beta)}
  
  if(method == 'BRIDGE') {p_mat = P_MAT_BRIDGE(beta)}
  
  if(method == 'SCAD')  {p_mat = P_MAT_SCAD(beta)}
  
  if((method != 'EN') & (method != 'BRIDGE') & (method != 'SCAD')) 
  {P_mat=''; print("incorrect method choose between 'EN' 'BRIDGE' 'SCAD' ")}
  
  p_mat
}


################################################################################

#DATA_LIST function puts data into another structure
DATA_LIST <- function(data)
{
  N <- max(data[,1])
  p <- length(data[1,])-2
  X_list  <- as.list(1:N)
  Y_list <- as.list(1:N)
  id <- data[,1]
  
  for (i in 1:N)
  {
    X_list[[i]] <- as.matrix(data[id==i,2:(p+1)])
    Y_list[[i]] <- as.matrix(data[id==i,length(data[1,])])
  }
  data_list <- list(X_list,Y_list,N,p)
  data_list
}

#############################################################################

##########
###Tuning parameter selection
#########


#Effective number of parameters (Hessian based)
#library(psych)
EF_PAR <- function(beta,method)
{
  H_mat <- H_MAT(beta)
  P_mat <- P_MAT(beta,method)
  
  #0 componenten eruit kegelen
  indic <- (abs(beta) > tresh_beta) 
  H_mat_red <- H_mat[indic==1,indic==1]
  P_mat_red <- P_mat[indic==1,indic==1]
  if(sum(indic)==0) {e <- 0}
  else  {  e <-  tr(solve(H_mat_red + n_obs*P_mat_red)%*%(H_mat_red ))  }    
  e
}

#Residuals  (pearson)
PEARS_RES <- function(beta)
{
  res <- Y_mat - X_mat%*%beta
  res
}
# probleem van standaardisatie: noodzakelijk bij BIC? anders groter Y variantie, invloed op penalty

#  BIC

BIC_DZIAK <- function(result)
{
  beta_N <- result[[1]]
  beta <- result[[2]]
  res <- PEARS_RES(beta)
  RSS <- t(res)%*%res
  dev <- log(RSS/n_obs)/n_obs
  e <- EF_PAR(beta_N,method);#aangepast hiervoor was beta i.p.v. beta N
  bic <- log(RSS/n_obs)*n_obs + log(n_obs)*e
  out <- list(bic,e,dev)
  names(out) <- list('BIC','e','dev')
  out
}


#  AIC

AIC_DZIAK <- function(result)
{
  beta_N <- result[[1]]
  beta <- result[[2]]
  res <- PEARS_RES(beta)
  RSS <- t(res)%*%res
  dev <- log(RSS/n_obs)/n_obs
  e <- EF_PAR(beta_N,method);#aangepast hiervoor was beta i.p.v. beta N
  aic <- log(RSS/n_obs)*n_obs + 2*e
  out <- list(aic,e,dev)
  names(out) <- list('AIC','e','dev')
  out
}

#fout in bovenstaande conde: e niet berekend op naïve par --> kan resultaten opnieuw verwerken

#aanpassing nodig: df op basis van naïve, res op basis van niet naïve 


#QGCV Fu , als tuning selection, use empircal correlation matrix, hessian based number of parameters
QGCV_FU <- function(result)
  
{
  beta_N <- result[[1]]
  beta <- result[[2]]
  res <- PEARS_RES(beta)
  W_mat <- diag(length(res));# implicatie working independence, wsch geen goede schatter
  Dev <- t(res)%*%W_mat%*%res
  e <-  EF_PAR(beta_N,method)
  QGCV <- Dev/(N*(1-e/n_obs))^2
  output <- list(QGCV,e,Dev) 
  names(output) <- c('QGCV','eff par','RSS')
  output
}



#############################################################################
##############################################################################

##########
###Cross_validatie: "Leave one subject out"   
##########




CV.PEE_LI <- function(data)
{
  p <- length(data[1,])-2
  names(data) <- c('id',paste('X',1:p,sep=""),'Y')  
  n_subj <- max(data[,1])
  
  beta_cv <- array(dim=c(n_subj/10,p))
  conv_cv <- array(dim = c(n_subj/10))
  RSS_cv <- array(dim=c(n_subj/10))
  
  Start <- 1
  Till <- n_subj/10
  
  for (ident in 1:(n_subj/10))
  {
    
    Till2 <- Till*ident
    ding <- seq(Start,Till2)
    #Trow one subject out
    data_min_ident <- subset(data,  !(id %in% ding))
    
    Start <- 1 +Till2
    #id opnieuw nummeren 
    greater <- (data_min_ident$id > ident)
    data_min_ident$id[greater] <-  data_min_ident$id[greater] - (n_subj/10)
    
    #Fit model
    result <- PEE_LI2(data_min_ident,method)  #aanpassing PEE_LI2 ipv PEE_LI
    
    
    #save results
    beta_cv[ident,] <- result[[2]]
    conv_cv[ident] <- result[[4]]
    
    data_ident <- data[data$id==ident,]
    #res <- data_ident[,p+2] -  as.matrix(data_ident[,2:(p+1)])%*%as.vector(beta_cv[ident,])
    PA <- as.matrix(data_ident[,2:(p+1)])%*%as.vector(beta_cv[ident,])
    PP <- exp(PA)/(1+exp(PA))
    #PD <- as.matrix(exp(data_ident[,2:(p+1)]%*%as.vector(beta_cv[ident,])))
    #PP <- as.matrix(exp(data_ident[,2:(p+1)]%*%as.vector(beta_cv[ident,]))/(1+exp(data_ident[,2:(p+1)])%*%as.vector(beta_cv[ident,])))
    #data_ident[,p+2][data_ident[,p+2]==0] <- -1
    #res <-  max(0, 1-data_ident[,p+2]*PP) 
    res <- data_ident[,p+2]-PP
    RSS_cv[ident] <- t(res)%*%res/length(res) #ok is scaled by number of obs for that subject
  }
  #output elements
  out <- list(beta_cv,conv_cv,RSS_cv)
  names(out) <- c('beta_cv','conv_cv','RSS_cv')
  prest <- mean(RSS_cv)
  prest_sde <- sd(as.vector(RSS_cv))/n_subj
  prest_Median <-  median(RSS_cv)    #meer robust tegen outliers
  output <- c(prest,prest_sde,prest_Median,out)
  names(output) <-  c('MSE_CV','sde','MedianSE_CV','beta_cv','conv_cv','RSS_cv')
  output  
}

#CV only suited as programmed here for the non-naive EN or SCAD_L2 estimator

##########
###Remarks
##########

#if a beta component drops below tresh_beta, the value of beta, its score and penalty part are put to zero (so no reappearance when a beta is put to zero)
#this not immediatly indicated in the paper of Li, possibly an addapation is possible(sort sensitivty) to reput a beta into the model

