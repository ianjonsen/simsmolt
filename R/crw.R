#' @title correlated random walk function
#' 
#' @description utility function not to be called by user
#' 
#' @importFrom CircStats rwrpcauchy
#' @importFrom stats rweibull
#' @importFrom raster extract cellFromXY adjacent xyFromCell
#' @export
#' 
crw <- function(n = 1, data, xy = NULL, buffer = NULL, mu, rho, a, b){
  
  if (a > 0)
    st <- rweibull(n, a, b)
  else if( a == 0 & b > 0) {
    st <- b
  } else {
    st <- 0
  }
  
  d2l <- extract(data$land, rbind(xy))
  
  
  if(d2l > buffer[1] | xy[1] < 300) {
    # nothing to do as mu, rho are supplied as args
    
  } else if (xy[1] >= 300 & !all(xy[1] >= 950, 
                                 xy[1] <= 1065, 
                                 xy[2] >= 1200, 
                                 xy[2] <= 1305) & 
             d2l <= buffer[1]) {
    
    ## direct smolt to move eastward & parallel to shore (avoid land) only after passing through most of GoM
    mu <- (extract(data$land_dir, rbind(xy)) + 0.5 * pi) %% (2*pi)
    rho <- 0.9 # deviate from random walk when inside land buffer
    if(d2l <= 2) {
      mu <- mu + 0.5 * pi %% (2*pi) ## move in opposite direction of land if within 2km
      rho <- 0.95
    }
  } else if(xy[1] >= 300 & all(xy[1] >= 950, 
                               xy[1] <= 1065, 
                               xy[2] >= 1200, 
                               xy[2] <= 1305) & 
            d2l <= buffer[1]) {
    mu <- (extract(data$land_dir, rbind(xy))  + 0.5 * pi) %% (pi)
    rho <- 0.9 # deviate from random walk when inside land buffer
    if(d2l <= 2) {
      mu <- mu + 0.5 * pi %% (pi) ## move in opposite direction of land if within 2km
      rho <- 0.95
    }
  }
  
  phi <- rwrpcauchy(n, mu, rho)  
  
  new.xy <- c(xy[1] + sin(phi) * st, xy[2] + cos(phi) * st)
  new.d2l <- extract(data$land, rbind(new.xy))
  
  ## if new location on land (0) then adjust so it's in water
  if(new.d2l == 0) {
    ## find all nearby cells (16) & select the one furthest from land
    adj.cells <- adjacent(data$land, 
                          cellFromXY(data$land, new.xy), 
                          directions = 16, 
                          pairs = FALSE)
    the.cell <- adj.cells[which(data$land[adj.cells] == max(data$land[adj.cells]))][1]
    new.xy <- xyFromCell(data$land, the.cell)
    ## get final distance from land
    new.d2l <- extract(data$land, rbind(new.xy))
  }
  
  return(c(new.xy, phi, new.d2l))
}