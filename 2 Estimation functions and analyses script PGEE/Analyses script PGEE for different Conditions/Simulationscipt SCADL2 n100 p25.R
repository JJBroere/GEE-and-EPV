#########################################################################################################################################################

source("Estimation and crosvalidation functions PGEE.R")
load('Sim100N100P25u.Rda')#simulated data


#generate grid of penalization parameters to evaluate
library(matlab)
n_lambda <-40
n_alpha <- 5
nsim <- 100
lambda_vec <- linspace(0.0025,0.5,n=n_lambda)
alpha_vec <- rev(linspace (0.001,1,n_alpha)) 
grid <- expand.grid(lambda_vec,alpha_vec)
names(grid) <- c('lambda' , 'alpha')
n_grid <- length(grid[,1])



#instellingen

M_it <- 1000 #(genoeg om convergentie te garanderen, indien convergentie stopt het sneller)
tresh_beta <- 0.00001 #for stability small enough, lots of coeff if a betaj closer than this treshold, set to zero (aanpassing: M_diff <<tresh_beta)
M_diff <- 0.0000000001 #(convergence indicator= precision)
method <- 'SCAD'

p <- 25 #the number of covariates
k <- 1 #ridge initialisation
a <- 3.7 #additional parameter, keep fixed as recommended buy Fan and Li
k <- 0 #initialisation parameter for first fit, if numarical problems with initialization use a larger value (parameter van initi�le RIDGE/DECORRELATIE FIT, put to zero or a small value if it doesnt work increase)
option <- '' #'adaptive'  to get adaptive LASSO/EN  and anything else to get nonadaptive version which is recommended
distr <- "binomial" #indicates the distribution used, options: "normal" or "binomial"


#objecten
MSE_CV <- array(dim=c(n_grid))
sde_CV <-  array(dim=c(n_grid))
MedianSE_CV <-  array(dim=c(n_grid))
n_subj <- 100  #aanpassen aan aantal subjecten
beta_cv <- array(dim=c(n_grid,n_subj,p))
conv_array <- array(dim=c(n_grid,n_subj))


EstimatesSCAD10025u <- matrix(1,25,nsim)
Optimal10025u <- matrix(1,nsim,2)

for(jj in 1:nsim){
  for (i_grid in 1:n_grid){
    lambda <- grid[i_grid,1]
    alpha <- grid[i_grid,2]
    
    resultaat <-  CV.PEE_LI(data=Sim100N100P25u[[jj]])
    MSE_CV[i_grid] <- resultaat$MSE_CV
    MedianSE_CV[i_grid] <- resultaat$MedianSE_CV
    sde_CV[i_grid] <- resultaat$sde
    beta_cv[i_grid,,] <- resultaat$beta_cv
    conv_array[i_grid,] <- resultaat$conv_cv
  }
  LOW_MSE <-which(MSE_CV %in% sort(MSE_CV)[1])
  
  Optimal10025u[jj,] <- as.numeric(grid[LOW_MSE[1],])
}
save(Optimal10025u, file="Optimal10025u.Rda")

for(jj in 1:nsim){
  alpha <- Optimal10025u[jj,2] #parameter 1 = lasso/SCAD, 0 =ridge controls the amount of penalty devided by sparceness and grouping effect producing penalties
  lambda <-Optimal10025u[jj,1]
  AAAA <-  PEE_LI2(data=Sim100N100P25u[[jj]],method='SCAD')
  
  EstimatesSCAD10025u[,jj] <- as.matrix(AAAA[[2]])
}

save(EstimatesSCAD10025u, file="EstimatesSCAD10025u.Rda")
rowMeans(EstimatesSCAD10025u)
ABBA<- EstimatesSCAD10025u
ABBA[ABBA==0]<-NA
rowMeans(ABBA, na.rm=T)