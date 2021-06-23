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
  
    ## regular movement
    ## calculate distance to land
    d2l <- terra::extract(data$land, rbind(xy))
    
    if (d2l > mpar$pars$buffer) {
      switch(mpar$scenario, 
             rs = { 
               if(i < mpar$pars$N/2) {
                 ## state 1: migration toward reconditioning area
                 ## if current xy not inside land buffer then employ biased migration
                 delta <- c(mpar$pars$coa[1,1] - xy[1], mpar$pars$coa[1,2] - xy[2])

                 mu <- atan2(delta[1], delta[2])
                 rho <- tanh(mpar$pars$r * sqrt(sum(delta^2)))
                 phi <- rwrpcauchy(1, mu, rho)
                 
               } else if (i >= mpar$pars$N/2) {
                 ## state 2: migration back to spawning river
                 mu <- atan2(mpar$pars$coa[2,1] - xy[1], mpar$pars$coa[2,2] - xy[2])
                 phi <- rwrpcauchy(1, mu, mpar$pars$rho)
               }
               },
             as = {
               if(i < round(mpar$pars$N * 0.85)) {
                 ## state 1: migration toward W Greenland
                 delta <- c(mpar$pars$coa[1,1] - xy[1], mpar$pars$coa[1,2] - xy[2])
                 
                 mu <- atan2(delta[1], delta[2])
                 rho <- tanh(mpar$pars$r * sqrt(sum(delta^2)))
                 phi <- rwrpcauchy(1, mu, rho)
               } else if (i >= round(mpar$pars$N * 0.85)) {
                 ## state 2: migration back to spawing river
                 mu <- atan2(mpar$pars$coa[2,1] - xy[1], mpar$pars$coa[2,2] - xy[2])
                 phi <- rwrpcauchy(1, mu, mpar$pars$rho)
               }
             })
      
#      if(mpar$growth) {
        ## Temperature-dependent direction reversal (instantaneous)
        g.rng <- growth(w, seq(6, 20, l = 100), s)
        ts.mig <- seq(6, 20, l = 100)[which(g.rng >= w)] %>% min() * 0.5 #0.75
        if(ts <= ts.mig) {
          phi <- ifelse(ts <= ts.mig, runif(1,pi-0.4,pi+0.4), phi)
        }
#      } else {
#        ## Temperature-dependent travel rate
#        s <- ifelse(ts <= 5, s * 0.1, s)
#      }
      
    } else {
      ## direct kelt to move parallel to land
      if(xy[2] < 1200) {
        mu <- (terra::extract(data$land_dir, rbind(xy)) + 0.5 * pi)
        if (d2l <= 2) {
          mu <-
            (mu + 0.5 * pi) ## move in opposite direction of land if within 2km
        }
      } else {
        ## if kelt near Greenland then move parallel to land in either direction
        sg <- sample(c(-1,1), size = 1)
        mu <- (terra::extract(data$land_dir, rbind(xy)) + sg * 0.5 * pi)
        if (d2l <= 2) {
          mu <-
            (mu + sg * 0.5 * pi) ## move in opposite direction of land if within 2km
        }
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