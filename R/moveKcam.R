#' @title random walk movement kernel for Campbellton River kelts
#' 
#' @description utility function not to be called by user
#' 
#' @importFrom CircStats rwrpcauchy
#' @importFrom stats rweibull
#' @importFrom raster extract xyFromCell
#' @export
#' 
moveKcam <- function(data, xy = NULL, mpar, i, s, ts, w) {
  
  ## if xy S of 1100 then apply fixed movement to get away from NL Islands
  # if(xy[2] <= 1100 & xy[1] >= 1100 & xy[1] <= 1230) {
  #   phi <- rwrpcauchy(1, -5/180*pi, 0.9)
  #   new.xy <- c(xy[1] + sin(phi) * s, xy[2] + cos(phi) * s)
  # 
  # } else {
    ## regular movement
    ## calculate distance to land
    d2l <- terra::extract(data$land, rbind(xy))
    
    if (d2l > mpar$pars$buffer) {
      switch(mpar$scenario, 
             rs = { 
               if(i < mpar$pars$N/2) {
                 ## state 1: migration toward Greenland
                 ## if current xy not inside land buffer then employ biased migration
                 phi <- rwrpcauchy(1, mpar$pars$mdir, mpar$pars$rho)
                 
               } else if (i >= mpar$pars$N/2) {
                 mu <- atan2(mpar$pars$coa[1] - xy[1], mpar$pars$coa[2] - xy[2])
                 phi <- rwrpcauchy(1, mu, mpar$pars$rho)
               }
               })
      
      if(mpar$growth) {
        ## Temperature-dependent direction reversal (instantaneous)
        g.rng <- growth(w, seq(6, 20, l = 100), s)
        ts.mig <- seq(6, 20, l = 100)[which(g.rng >= w)] %>% min() * 0.75
        if(ts <= ts.mig) {
          phi <- ifelse(ts <= ts.mig, runif(1,pi-0.4,pi+0.4), phi)
          cat("\ncold water")
        }
      } else {
        ## Temperature-dependent travel rate
        s <- ifelse(ts <= 5, s * 0.1, s)
      }
      
    } else {
      ## direct kelt to move parallel to land
      mu <- (terra::extract(data$land_dir, rbind(xy)) + 0.5 * pi)
      if (d2l <= 2) {
        mu <-
          (mu + 0.5 * pi) ## move in opposite direction of land if within 2km
      }
      phi <- rwrpcauchy(1, mu, 0.95)
    }   
    
    new.xy <- c(xy[1] + sin(phi) * s, xy[2] + cos(phi) * s)
    if(mpar$shelf) {
      pv <- c(extract(data$shelf[[1]], rbind(new.xy))[1],
              extract(data$shelf[[2]], rbind(new.xy))[1])
      new.xy <- new.xy + pv * mpar$pars$beta
    }
    
    ## check if new xy close/on land
    new.d2l <- terra::extract(data$land, rbind(new.xy))
    
    ## if new location on land (0) then adjust so it's in water
    if(!is.na(new.d2l) & new.d2l == 0) {
      ## find all nearby cells within 3 km & select the one farthest from land
      cells <- terra::extract(data$land, rbind(new.xy), buffer = 3, cellnumbers = TRUE, df = TRUE)
      cell.max <- cells[cells[, 3] == max(cells[, 3], na.rm = TRUE), 2][1]
      new.xy <- xyFromCell(data$land, cell.max)
      
    } else if(is.na(new.d2l)) {
      new.xy <- c(NA,NA)
    } 
#  }
  
  cbind(new.xy[1], new.xy[2], phi %% (2*pi))
  
}