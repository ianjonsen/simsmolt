#' @title simulate movements & survival, using inputs from presim
#' 
#' @description simulates, movement & detection of acoustically-tagged animals
#' 
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#' 
#' @param N - number of 1-h time steps to simulate
#' @param tag - start location(s) of simulated animals
#' @param coa - optional Centre-Of-Attraction location(s) to provide movement bias(es)
#' @param params - movement & (otpional) survival parameters supplied as a list, see details
#' @param psim - a list of required data from \code{presim}
#' @importFrom raster extract
#' @importFrom dplyr %>%
#' @importFrom tibble as_tibble
#' @importFrom stats runif rbinom
#' @export
#' 
sim_move <-
  function(N = 1800,
           tag = c(150, 240),
           coa = c(900, 950),
           chgpt = c(800, 835),
           psim = NULL,
           param = list(
             a = 2,
             b = 0.864,
             rho = 0.8,
             ntries = 100,
             surv = 0.9936
           )) {
    
    if (is.null(psim))
      stop("A `presim` object must be supplied\n")
    if(class(psim$bathy)[1] != "RasterLayer") stop("bathymetry must be a RasterLayer")
    if(!is.null(psim$uv)) {
      if(class(psim$uv)[1] != "RasterBrick" || (!names(psim$uv) %in% c("u","v"))) 
        stop("current data must be supplied as a RasterBrick with u and v layers")
    }
    if(is.null(chgpt)) chgpt <- c(100000, 100000) # set to huge number to ensure all steps are before chgpt
    
    ## define location/survivourship, s matrix & start position
    X <- matrix(NA, N, 3)
    s <- matrix(NA, N, 2)
    X[1,] <- cbind(tag[1], tag[2], 1)
    delta <- theta <- ps <- c()
    thetaL <- runif(1, -20, 10)
    rhoL <- runif(1, 0.6, 0.8)
  
    ## recursion
    for (i in 2:N) {
      if(X[i - 1, 1] < chgpt[1] & X[i - 1, 2] < chgpt[2]) {
        delta[i] <- sqrt((coa[1] - X[i - 1, 1]) ^ 2 + (coa[2] - X[i - 1, 2]) ^ 2)
        theta[i] <- atan2(coa[1] - X[i - 1, 1], coa[2] - X[i - 1, 2])
      } else {
        ## if in Lab Sea then head NW
        delta[i] <- NA
        theta[i] <- thetaL / 180*pi
        rho <- rhoL
      }
    
      ## draw n proposal steps
      d <- rw(n = param$ntries, mu = theta[i], rho = param$rho, a = param$a, b = param$b) 
      
      ## add advection due to current, if raster is supplied
      if(is.null(psim$uv)) {
        tmp <- cbind(X[i - 1, 1] + d[, 1], X[i - 1, 2] + d[, 2])
      } else {
        uv <- extract(psim$uv, rbind(X[i-1, 1:2])) * 3600/1000 # to convert from m/s to km/h
        tmp <- cbind(X[i - 1, 1] + d[, 1] + uv[1], X[i - 1, 2] + d[, 2] + uv[2])
      }
     
      ## check if proposed update is on land or in deep ocean
      p.step <- extract(psim$bathy, tmp, method = "simple")
      idx <- which(p.step == max(p.step))[1]
      ps[i] <- p.step[p.step == max(p.step)][1]
      X[i, 1:2] <- tmp[idx,]
      X[i, 3] <- rbinom(1, 1, param$surv^(1/24)) # rescale daily survival to hourly prob. 
      if(X[i, 3] == 0) break
    }
 
    X <- data.frame(X)
    names(X) <- c("x", "y", "s")
    phi <- c(NA, diff(theta))
    fill <- rep(NA, nrow(X) - length(phi))
    
    sim <- data.frame(X, ps=c(ps,fill), delta=c(delta,fill), theta=c(theta,fill), phi=c(phi,fill)) %>%
      as_tibble()
    
    
    psim$tag <- tag
    psim$coa <- coa
    psim$chgpt <- chgpt
  
    out <- list(sim = sim, data = psim)  
    class(out) <- "simsmolt"
    
    return(out)
    
}
