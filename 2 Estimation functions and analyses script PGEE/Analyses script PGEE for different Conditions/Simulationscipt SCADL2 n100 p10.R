#########################################################################################################################################################
setwd("C:/Users/Joris Broere/Desktop/Gebruik map 100")
source("Estimation and crosvalidation functions PGEE.R")
load('Sim100N100P10u.Rda')#simulated data

#generate grid of penalization parameters to evaluate
library(matlab)
    n_lambda <-40
    n_alpha <- 5
    nsim <- 100
    lambda_vec <- linspace(0.0025,0.5,n=n_lambda)
    alpha_vec <- rev(linspace (0.001,1,n_alpha)) #beginnen waarmee meest problemen verwacht
    grid <- expand.grid(lambda_vec,alpha_vec)
    names(grid) <- c('lambda' , 'alpha')
    n_grid <- length(grid[,1])
    
    
    
    #instellingen
    
    M_it <- 1000 #(genoeg om convergentie te garanderen, indien convergentie stopt het sneller)
    tresh_beta <- 0.00001 #for stability small enough, lots of coeff if a betaj closer than this treshold, set to zero (aanpassing: M_diff <<tresh_beta)
    M_diff <- 0.0000000001 #(convergence indicator= precision)
    method <- 'SCAD'
    
    p <- 10 #the number of covariates
    k <- 1 #ridge initialisation
    a <- 3.7 #additional parameter, keep fixed as recommended buy Fan and Li
    #V
    # numaric constants to stop the algorithm is the norm of the difference in estimated beta vector is smaller than this constant(convergence indicator= precision)
    k <- 0 #initialisation parameter for first fit, if numarical problems with initialization use a larger value (parameter van initiële RIDGE/DECORRELATIE FIT, put to zero or a small value if it doesnt work increase)
    option <- '' #'adaptive'  to get adaptive LASSO/EN  and anything else to get nonadaptive version which is recommended
    distr <- "binomial" #indicates the distribution used, options: "normal" or "binomial"
    
    
    #objecten
    MSE_CV <- array(dim=c(n_grid))
    sde_CV <-  array(dim=c(n_grid))
    MedianSE_CV <-  array(dim=c(n_grid))
    n_subj <- 100  #aanpassen aan aantal subjecten
    beta_cv <- array(dim=c(n_grid,n_subj,p))
    conv_array <- array(dim=c(n_grid,n_subj))
    

EstimatesSCAD10010u <- matrix(1,10,nsim)
  Optimal10010u <- matrix(1,nsim,2)
      
for(jj in 1:nsim){
    for (i_grid in 1:n_grid){
      lambda <- grid[i_grid,1]
      alpha <- grid[i_grid,2]
      
      resultaat <-  CV.PEE_LI(data=Sim100N100P10u[[jj]])
      MSE_CV[i_grid] <- resultaat$MSE_CV
      MedianSE_CV[i_grid] <- resultaat$MedianSE_CV
      sde_CV[i_grid] <- resultaat$sde
      beta_cv[i_grid,,] <- resultaat$beta_cv
      conv_array[i_grid,] <- resultaat$conv_cv
    }
    LOW_MSE <-which(MSE_CV %in% sort(MSE_CV)[1])
    
    Optimal10010u[jj,] <- as.numeric(grid[LOW_MSE[1],])
    }
save(Optimal10010u, file="Optimal10010u.Rda")
   
  for(jj in 1:nsim){
    alpha <- Optimal10010u[jj,2] #parameter 1 = lasso/SCAD, 0 =ridge controls the amount of penalty devided by sparceness and grouping effect producing penalties
    lambda <-Optimal10010u[jj,1]
    AAAA <-  PEE_LI2(data=Sim100N100P10u[[jj]],method='SCAD')
    
    EstimatesSCAD10010u[,jj] <- as.matrix(AAAA[[2]])
  }

save(EstimatesSCAD10010u, file="EstimatesSCAD10010u.Rda")
rowMeans(EstimatesSCAD10010u)
ABA <-EstimatesSCAD10010u
ABA[ABA== 0] <- NA
rowMeans(ABA, na.rm=T)