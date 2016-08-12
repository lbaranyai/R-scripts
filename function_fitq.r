## R statistics
## Report quality measures of curve fitting
## CC 3.0 - BY NC SA - by lbaranyai@github
#
## Function parameters are
## - observed: data vector
## - predicted: values produced by model
## - np: number of function coefficients (optional)
#
fitq <- function(observed,predicted,np=0)
{
  N <- length(observed)
  # Coefficient of determination
  tmp <- cor(observed,predicted)^2
  cat("R-squared:",tmp,"\n")
  # Adjusted value to balance model complexity
  if (np>0) {
    tmp <- 1 - (1-tmp)*(N-1)/(N-np)
    cat("adjusted R-squared:",tmp,"\n")
  }
  # ANOVA F value
  if (np>0) {
    my <- mean(predicted)
    sse <- sum( (y-predicted)^2 )
    ssm <- sum( (predicted-my)^2 )
    tmp <- (ssm/(np-1))/(sse/(N-np))
    cat("F value:",tmp,"for df",np-1,N-np,"\n")
  }
  diff <- observed - predicted
  # Mean Squared Error
  tmp <- mean(diff^2)
  cat("MSE:",tmp,"\n")
  # Root Mean Squared Error
  tmp <- sqrt(tmp)
  cat("RMSE:",tmp,"\n")
  # Akaike Information Criterion
  if (np>0) {
    tmp <- 2*(np-1) + N*log(mean(diff^2))
    cat("AIC:",tmp,"\n")
  }
  # Durbin-Watson test for autocorrelation
  tmp <- 0
  for (i in 2:N) {
    tmp <- tmp + (diff[i]-diff[i-1])^2
  }
  tmp <- tmp/sum(diff^2)
  cat("Durbin-Watson D:",tmp,"\n")
}
