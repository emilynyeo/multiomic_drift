predict.longituRF <- function(object, X,Z,id,time,...){
  # This function predicts values for longitudinal data using a random forest model with random effects
  # It handles different types of stochastic processes (none, exponential, fractional Brownian motion, or other)
  # 
  # Args:
  #   object: A fitted longitudinal random forest model
  #   X: Matrix of predictor variables
  #   Z: Matrix of random effects covariates
  #   id: Vector of subject IDs
  #   time: Vector of time points
  #
  # Returns:
  #   Vector of predicted values incorporating fixed effects, random effects and stochastic processes
  
  # Get number of unique subjects
  n <- length(unique(id))
  id_btilde <- object$id_btilde
  # Get predictions from random forest for fixed effects
  f <- predict(object$forest,X)
  Time <- object$time
  id_btilde <- object$id_btilde
  # Initialize predictions vector
  Ypred <- rep(0,length(id))
  id.app=object$id
  
  # Case 1: No stochastic process
  if (object$sto=="none"){
    for (i in 1:length(unique(id))){
      # Get indices for current subject
      w <- which(id==unique(id)[i])
      k <- which(id_btilde==unique(id)[i])
      # Combine fixed and random effects
      Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,]
    }
    return(Ypred)
  }

  # Case 2: Exponential stochastic process
  if (object$sto=="exp"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      k <- which(id_btilde==unique(id)[i])
      om <- which(id.app==unique(id)[i])
      # Add exponential process prediction to fixed and random effects
      Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,] + predict.exp(object$omega[om],Time[om],time[w], object$alpha)
    }
    return(Ypred)
  }

  # Case 3: Fractional Brownian motion process
  if (object$sto=="fbm"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      k <- which(id_btilde==unique(id)[i])
      om <- which(id.app==unique(id)[i])
      # Add FBM process prediction to fixed and random effects
      Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,] + predict.fbm(object$omega[om],Time[om],time[w], object$Hurst)
    }
    return(Ypred)
  }

  # Case 4: Other stochastic process
  for (i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    k <- which(id_btilde==unique(id)[i])
    om <- which(id.app==unique(id)[i])
    # Add generic stochastic process prediction to fixed and random effects
    Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,] + predict.sto(object$omega[om],Time[om],time[w], object$sto)
  }
  return(Ypred)
}
