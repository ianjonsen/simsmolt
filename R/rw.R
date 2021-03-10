#' @title random walk function
#' 
#' @description utility function not to be called by user
#' 
#' @importFrom CircStats rwrpcauchy
#' @importFrom stats rweibull
#' @importFrom raster extract xyFromCell
#' @export
#' 
rw <- function(n = 1, data, xy = NULL, buffer = NULL, rho, a, b, taxis, u, v){
  
  if (a > 0)
    st <- rweibull(n, a, b)
  else if( a == 0 & b > 0) {
    st <- b
  } else {
    st <- 0
  }
  
  d2l <- extract(data$land, rbind(xy))
  
  
  if(d2l > buffer[1] | xy[1] < 300) {
    mu <- runif(1, -pi, pi)
    
  } else if (xy[1] >= 300 & !all(xy[1] >= data$sobi.box[1], 
                                 xy[1] <= data$sobi.box[1], 
                                 xy[2] >= data$sobi.box[1], 
                                 xy[2] <= data$sobi.box[1]) & 
             d2l <= buffer[1]) {
    
    ## direct smolt to move eastward & parallel to shore (avoid land) only after passing through most of GoM
    mu <- (extract(data$land_dir, rbind(xy)) + 0.5 * pi) %% (2*pi)
    rho <- 0.9 # deviate from random walk when inside land buffer
    if(d2l <= 2) {
      mu <- mu + 0.5 * pi %% (2*pi) ## move in opposite direction of land if within 2km
      rho <- 0.95
    }
  } else if(xy[1] >= 300 & all(xy[1] >= data$sobi.box[1], 
                               xy[1] <= data$sobi.box[1], 
                               xy[2] >= data$sobi.box[1], 
                               xy[2] <= data$sobi.box[1]) & 
            d2l <= buffer[1]) {
    mu <- (extract(data$land_dir, rbind(xy))  + 0.5 * pi) %% (pi)
    rho <- 0.9 # deviate from random walk when inside land buffer
    if(d2l <= 2) {
      mu <- mu + 0.5 * pi %% (pi) ## move in opposite direction of land if within 2km
      rho <- 0.95
    }
  }
  
  phi0 <- rwrpcauchy(n, mu, rho)
  
  ## account for any rheotaxis by modifying the step bearing phi to orient against current (taxis = "p") or with current (taxis = "n")
  ## strength of rheotaxis (adherence to the current direction) is a function of the ratio of current magnitude to active movement step length ()
  if(!is.na(taxis)) {
    switch(taxis,
           p = {
             rhb <- atan2(u,v) # against current
             rhm <- sqrt(u^2+v^2)
             rh.wt <- ifelse(rhm / st > 1, 1, rhm / st)
             phi <- atan2(sum(sin(rhb * rh.wt), sin(phi0 * (1 - rh.wt))), sum(cos(rhb * rh.wt), cos(phi0 * (1 - rh.wt))))
           },
           n = {
             rhb <- atan2(v,u) # with current
             rhm <- sqrt(u^2+v^2)
             rh.wt <- ifelse(rhm / st > 1, 1, rhm / st)
             phi <- atan2(sum(sin(rhb * rh.wt), sin(phi0 * (1 - rh.wt))), sum(cos(rhb * rh.wt), cos(phi0 * (1 - rh.wt))))
           })
  } else {
    phi <- phi0
  }
  
  
  new.xy <- c(xy[1] + sin(phi) * st, xy[2] + cos(phi) * st)
  new.d2l <- extract(data$land, rbind(new.xy))
  
  ## if new location on land (0) then adjust so it's in water
  if(!is.na(new.d2l) & new.d2l == 0) {
    ## find all nearby cells within 3 km & select the one farthest from land
    cells <- extract(data$land, rbind(new.xy), buffer = 3, cellnumbers = TRUE, df = TRUE)
    cell.max <- cells[cells[, "layer"] == max(cells[, "layer"]), "cells"][1]
    new.xy <- xyFromCell(data$land, cell.max)
  
    return(new.xy)
    
  } else if(is.na(new.d2l)) {
    return(cbind(NA,NA))
    
  } else if(!is.na(new.d2l) & new.d2l > 0) {
    return(new.xy)
  }
}