#' @title random walk function
#' 
#' @description utility function not to be called by user
#' 
#' @importFrom CircStats rwrpcauchy
#' @importFrom stats rweibull
#' @export
#' 
rw <- function(n = 1, mu = 0, rho = 0.5, a = 3, b = 5){
  
  if (a > 0)
    st <- rweibull(n, a, b)
  else {
    st <- b
  }
  
  phi <- rwrpcauchy(n, mu, rho)

  cbind(sin(phi) * st, cos(phi) * st)
} 
