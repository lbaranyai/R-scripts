## R statistics
## Perform Durbin-Watson test of autocorrelation
## CC 3.0 - BY NC SA - by lbaranyai@github
#
## Function parameters are
## - values: data vector of residuals, time serie, etc
## - offset: boolean, baseline correction if necessary (default=false)
#
DW.test <- function(values,offset=FALSE)
{
  N <- length(values)
  tmp <- values
  # Baseline correction upon request
  if (offset==TRUE) {
   tmp <- tmp - mean(tmp)
  }
  # Calculate Durbin-Watson value
  DW <- 0
  for (i in 2:N) {
    DW <- DW + (tmp[i]-tmp[i-1])^2
  }
  DW <- DW/sum(tmp^2)
  cat("Durbin-Watson D:",DW,"for df",N-1,"\n")
}
