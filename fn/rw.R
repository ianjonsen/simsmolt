rw <- function(n = 1, mu = 0, rho = 0.5, a = 3, b = 5){
  
  s <- ifelse(a > 0, rweibull(n, a, b), b)
  
  phi <- CircStats::rwrpcauchy(n, mu, rho)
  
  cbind(sin(phi) * s, cos(phi) * s)
} 
