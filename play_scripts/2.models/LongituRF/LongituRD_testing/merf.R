#' Mixed Effects Random Forest (MERF) for Longitudinal Data Analysis
#' 
#' This function implements a Mixed Effects Random Forest model for longitudinal data,
#' combining random forest with mixed effects modeling. It handles both fixed and random effects,
#' and can incorporate different types of stochastic processes (fractional Brownian motion,
#' exponential correlation, or none) to model temporal correlations.
#'
MERF <- function(X,Y,id,Z,iter=100,mtry=ceiling(ncol(X)/3),ntree=500, time, sto, delta = 0.001){
  # Initialize dimensions and parameters
  q <- dim(Z)[2]
  nind <- length(unique(id))
  btilde <- matrix(0,nind,q) # Random effects matrix - one row per individual
  sigmahat <- 1              # Initial error variance
  Btilde <- diag(rep(1,q))  # Initial random effects covariance matrix
  epsilonhat <- rep(0,length(Y))  # Residuals
  id_btilde <- unique(id)
  Tiime <- sort(unique(time))
  omega <- rep(0,length(Y))  # Stochastic process values
  sigma2 <- 1                # Initial stochastic process variance
  Vrai <- NULL              # Store likelihood values
  inc <- 1                  # Convergence increment
  OOB <- NULL              # Store Out-of-bag errors
  
  # Handle different types of stochastic processes
  if (class(sto)=="character"){
    # Case 1: Fractional Brownian Motion
    if (sto=="fbm"){
      # Create mapping matrix for individuals to timepoints
      id_omega <- matrix(0,nind,length(unique(time)))
      for (i in 1:length(unique(id))){
        w <- which(id ==id_btilde[i])
        time11 <- time[w]
        where <- NULL
        for (j in 1:length(time11)){
          where <- c(where,which(Tiime==time11[j]))
        }
        id_omega[i,where] <- 1
      }
      
      omega <- matrix(0,nind,length(unique(time)))
      omega2 <- rep(0,length(Y))
      # Optimize Hurst parameter for FBM
      h <- opti.FBM(X,Y,id,Z,iter, mtry,ntree,time)
      
      # Main iteration loop
      for (i in 1:iter){
        # Remove random effects and stochastic process from response
        ystar <- rep(0,length(Y))
        for (k in 1:nind){
          indiv <- which(id==unique(id)[k])
          ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
        }
        
        # Fit random forest to adjusted response
        forest <- randomForest(X,ystar,mtry=mtry,ntree=ntree, 
                               importance = TRUE)
        fhat <- predict(forest)
        OOB[i] <- forest$mse[ntree]
        
        # Update random effects and stochastic process for each individual
        for (k in 1:nind){
          indiv <- which(id==unique(id)[k])
          K <- cov.fbm(time[indiv], h)  # FBM covariance matrix
          # Total covariance matrix
          V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+
                 diag(as.numeric(sigmahat),length(indiv),
                      length(indiv))+ sigma2*K
          # Update random effects
          btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          # Update stochastic process values
          omega[k,which(id_omega[k,]==1)] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          omega2[indiv] <- omega[k,which(id_omega[k,]==1)]
          # Calculate residuals
          epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[k,which(id_omega[k,]==1)]
        }
        
        # Store current parameter values
        sigm <- sigmahat
        B <- Btilde
        
        # Update variance components
        sigmahat <- sig.fbm(Y,sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,h)
        Btilde  <- bay.fbm(btilde,Btilde,Z,id,sigm, time, sigma2,h)
        sigma2 <- gam_fbm(Y,sigm,id,Z,B,time,sigma2,omega,id_omega,h)
        
        # Calculate likelihood and check convergence
        Vrai <- c(Vrai, logV.fbm(Y,fhat,Z[,,drop=FALSE],
                                 time,
                                 id,
                                 Btilde,
                                 sigma2,
                                 sigmahat,h))
        if (i>1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
        if (inc < delta) {
          print(paste0("stopped after ", i, " iterations."))
          sortie <- list(forest=forest,random_effects=btilde,
                         var_random_effects=Btilde,sigma=sigmahat,
                         sigma_sto=sigma2, id_btilde=unique(id), sto= sto, 
                         vraisemblance = Vrai,id=id, time =time, 
                         Hurst=h, OOB =OOB, omega=omega2)
          class(sortie)<-"longituRF"
          return(sortie)
        }
      }
      
      # Return final results if max iterations reached
      sortie <- list(forest=forest,random_effects=btilde,
                     var_random_effects=Btilde,sigma=sigmahat, 
                     id_btilde=unique(id),sigma_sto=sigma2,omega=omega2, 
                     sigma_sto =sigma2, time = time, sto= sto, Hurst =h, id=id, 
                     Vraisemblance=Vrai, OOB =OOB)
      class(sortie) <- "longituRF"
      return(sortie)
    }
    
    # Case 2: Exponential correlation structure
    if (sto=="exp"){
      # Similar structure to FBM case but with exponential correlation
      # Create mapping matrix for individuals to timepoints
      id_omega <- matrix(0,nind,length(unique(time)))
      for (i in 1:length(unique(id))){
        w <- which(id ==id_btilde[i])
        time11 <- time[w]
        where <- NULL
        for (j in 1:length(time11)){
          where <- c(where,which(Tiime==time11[j]))
        }
        id_omega[i,where] <- 1
      }
      
      omega <- matrix(0,nind,length(unique(time)))
      omega2 <- rep(0,length(Y))
      # Optimize alpha parameter for exponential correlation
      alpha <- opti.exp(X,Y,id,Z,iter, mtry,ntree,time)
      
      # Main iteration loop - similar structure to FBM case
      for (i in 1:iter){
        ystar <- rep(0,length(Y))
        for (k in 1:nind){
          indiv <- which(id==unique(id)[k])
          # Editted below 
          ystar[indiv] <- Y[indiv]- as.vector(Z[indiv,, drop=FALSE] %*% matrix(btilde[k,], ncol=1))- omega[indiv]
        }
        
        forest <- randomForest(X,ystar,mtry=mtry,ntree=ntree, importance = TRUE)
        fhat <- predict(forest)
        OOB[i] <- forest$mse[ntree]
        
        for (k in 1:nind){
          indiv <- which(id==unique(id)[k])
          K <- cov.exp(time[indiv], alpha)
          V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+
                diag(as.numeric(sigmahat),length(indiv),length(indiv))+ sigma2*K
          btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          omega[k,which(id_omega[k,]==1)] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          omega2[indiv] <- omega[k,which(id_omega[k,]==1)]
          epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
        }
        
        sigm <- sigmahat
        B <- Btilde
        sigmahat <- sig.exp(Y,
                            sigmahat,id, 
                            Z, 
                            epsilonhat, 
                            Btilde, 
                            time, 
                            sigma2,
                            alpha)
        Btilde  <- bay.exp(btilde,Btilde,Z,id,sigm, time, sigma2,alpha)
        sigma2 <- gam_exp(Y,sigm,id,Z,B,time,
                          sigma2,omega,id_omega,alpha)
        
        Vrai <- c(Vrai,logV.exp(Y,fhat,Z[,,drop=FALSE],
                                time,id,Btilde,sigma2,
                                sigmahat,alpha))
        if (i>1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
        if (inc < delta) {
          print(paste0("stopped after ", i, " iterations."))
          sortie <- list(forest=forest,
                         random_effects=btilde,
                         var_random_effects=Btilde,
                         sigma=sigmahat, 
                         id_btilde=unique(id), 
                         sto= sto, vraisemblance = 
                           Vrai,id=id, time=time, 
                         alpha = alpha, OOB =OOB, 
                         omega=omega2)
          class(sortie) <- "longituRF"
          return(sortie)
        }
      }
      
      sortie <- list(forest=forest,
                     random_effects=btilde,
                     var_random_effects=Btilde,
                     sigma=sigmahat, id_btilde=unique(id), 
                     omega=omega2, sigma_sto =sigma2, 
                     time = time, sto= sto, alpha=alpha, 
                     id=id, Vraisemblance=Vrai, OOB =OOB)
      class(sortie) <- "longituRF"
      return(sortie)
    }
    
    # Case 3: No stochastic process
    if ( sto=="none"){
      # Simpler iteration loop without stochastic process components
      for (i in 1:iter){
        ystar <- rep(NA,length(Y))
        for (k in 1:nind){
          indiv <- which(id==unique(id)[k])
          ystar[indiv] <- Y[indiv]- Z[indiv,,drop=FALSE]%*%btilde[k,]
        }
        
        forest <- randomForest(X,ystar,mtry=mtry,ntree=ntree, importance = TRUE)
        fhat <- predict(forest)
        OOB[i] <- forest$mse[ntree]
        
        for (k in 1:nind){
          indiv <- which(id==unique(id)[k])
          V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+
            diag(as.numeric(sigmahat),length(indiv),length(indiv))
          btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]
        }
        
        sigm <- sigmahat
        sigmahat <- sig(sigma = sigmahat,id = id, Z = Z, epsilon = epsilonhat,Btilde = Btilde)
        Btilde  <- bay(bhat = btilde,Bhat = Btilde,Z = Z,id = id,sigmahat = sigm)
        Vrai <- c(Vrai, logV(Y,fhat,Z,time,id,Btilde,0,sigmahat,sto))
        if (i>1) inc <- abs((Vrai[i-1]-Vrai[i])/Vrai[i-1])
        
        if (inc < delta) {
          print(paste0("stopped after ", i, " iterations."))
          sortie <- list(forest=forest,
                         random_effects=btilde,
                         var_random_effects=Btilde,sigma=sigmahat, 
                         id_btilde=unique(id), sto= sto, 
                         vraisemblance = Vrai,id=id, time=time, OOB =OOB)
          class(sortie) <- "longituRF"
          return(sortie)
        }
      }
      
      sortie <- list(forest=forest,
                     random_effects=btilde,
                     var_random_effects=Btilde,
                     sigma=sigmahat, 
                     id_btilde=unique(id), sto= sto, 
                     vraisemblance=Vrai,id=id, time=time, OOB =OOB)
      class(sortie) <- "longituRF"
      return(sortie)
    }
  }
  
  # Default case: Custom stochastic process
  for (i in 1:iter){
    ystar <- rep(0,length(Y))
    for (k in 1:nind){
      indiv <- which(id==unique(id)[k])
      ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
    }
    
    forest <- randomForest(X,ystar,mtry=mtry,ntree=ntree, importance=TRUE)
    fhat <- predict(forest)
    OOB[i] <- forest$mse[ntree]
    
    for (k in 1:nind){
      indiv <- which(id==unique(id)[k])
      K <- sto_analysis(sto,time[indiv])
      V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+
             diag(as.numeric(sigmahat),
                  length(indiv),
                  length(indiv))+ sigma2*K
      btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      omega[indiv] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
    }
    
    sigm <- sigmahat
    B <- Btilde
    sigmahat <- sig_sto(sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,sto)
    Btilde  <- bay_sto(btilde,Btilde,Z,id,sigm, time, sigma2,sto)
    sigma2 <- gam_sto(sigm,id,Z,B,time,sigma2,sto,omega)
    
    Vrai <- c(Vrai, logV(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,sto))
    if (i>1) inc <- abs((Vrai[i-1]-Vrai[i])/Vrai[i-1])
    if (inc < delta) {
      print(paste0("stopped after ", i, " iterations."))
      sortie <- list(forest=forest,random_effects=btilde,
                     var_random_effects=Btilde,sigma=sigmahat, 
                     id_btilde=unique(id), omega=omega, 
                     sigma_sto =sigma2, time = time, 
                     sto= sto,Vraisemblance=Vrai,id=id, OOB =OOB)
      class(sortie) <- "longituRF"
      return(sortie)
    }
  }
  
  sortie <- list(forest=forest,random_effects=btilde,
                 var_random_effects=Btilde,sigma=sigmahat, 
                 id_btilde=unique(id),omega=omega, 
                 sigma_sto =sigma2, time = time, 
                 sto= sto,Vraisemblance=Vrai,id=id, OOB =OOB)
  class(sortie) <- "longituRF"
  return(sortie)
}
