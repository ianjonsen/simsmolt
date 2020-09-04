#' @title simulate simple movements & survival, using inputs from sim_setup
#' 
#' @description simulates movement
#' 
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#' 
#' @param id - identifier for simulation run (individual animal)
#' @param N - number of 1-h time steps to simulate
#' @param data - a list of required data from \code{presim}
#' @param mpar - simulation control parameters supplied as a list, see details
#' @param pb - use progress bar (logical)
#' @importFrom raster extract
#' @importFrom CircStats rwrpcauchy
#' @importFrom dplyr %>%
#' @importFrom tibble as_tibble
#' @importFrom stats runif rbinom
#' @export
#' 
sim_simple <-
  function(id=1, 
           N = 2760,
           data = NULL,
           mpar = list(),
           pb = TRUE
           ) {
    ## default move parameters
    
    mpar.full <- list(
      start = "sobi", # or "nf" - postsmolts start in SoBI or E of northern tip of NF
      move = "sc", # or "uc" - shelf-constrained or unconstrained movements
      bias = function(x) (cos(x / 350) + sin(x / 1933)) * 45, # coa back to a single attraction point (for now)
      a = 2,
      b = 0.864,
      rho = 0.8,
      ntries = 500,
      surv = 0.9936,
      buffer = c(10, 10),
      mindist = 8,
      maxdist = 850, 
      pdrf = c(3.017, -0.0139) #c(4.865, -0.0139) = p(0.5) @ 350 m (~ consistent w HFX line V9 @ high power) 
    )
    pnms <- names(mpar.full)
    
    if (is.null(data))
      stop("Can't find output from sim_setup()\n")
    if (class(data$land)[1] != "RasterLayer") stop("distance2land must be a RasterLayer")
    
    if (length(mpar)) {
      nms <- names(mpar)
      if (!is.list(mpar) || is.null(nms))
        stop("'mpar' argument must be a named list")
      pos <- pmatch(nms, pnms)
      if (any(nap <- is.na(pos))) {
        warning(sprintf(ngettext(length(nap), "unrecognized mpar element named %s ignored", 
                                 "unrecognized mpar elements named %s ignored"), 
                        paste(sQuote(nms[nap]), collapse = ", ")), domain = NA)
        pos <- pos[!nap]
        mpar <- mpar[!nap]
      }
      mpar.full[pos] <- mpar
    } 
    mpar <- mpar.full
    
 

    sloc <- switch(mpar$start, 
           sobi = c(560, 82.5),
           nf = c(sample(seq(675, 725, l=50), 1), 82.5)
    )
    
    ## define location/survivourship, sw matrix & start position
    X <- matrix(NA, N, 3)
    #sw <- matrix(NA, N, 2)
    X[1,] <- cbind(sloc[1], sloc[2], 1)
    theta_s <- rd <- rz <- c()
    theta_s[1] <- 0 / 180 * pi

    ## recursion
    for (i in 2:N) {
      if(i==2 && pb)  tpb <- txtProgressBar(min = 2, max = N, style = 3)
      
      ## get distance from land at first step
      if(i == 2) {
        rd[i-1] <- raster::extract(data$land, rbind(X[i - 1, 1:2]))
        rz[i-1] <- raster::extract(data$b900.dist, rbind(X[i - 1, 1:2]))
      }
      ## heading to move ~ parallel to coastline
      theta_s[i] <- mpar$bias(X[i - 1, 2]) / 180*pi
    
      ## move vectors
      ## draw n proposal steps for active swimming
      d_s <-
        rw(
          n = mpar$ntries,
          mu = theta_s[i],
          rho = mpar$rho,
          a = mpar$a,
          b = mpar$b
        ) 
      d_b <- cbind(0,0)
      d_c <- cbind(0,0)
      
      ## apply weighting based on current distance to coast, 900 m isobath & specified buffer distance
      ##  first buffer's influence starts after exiting SoBI and after influence of first coa
      if(rd[i-1] <= mpar$buffer[1] & X[i - 1, 2] > 100) {
        ## direction ~ parallel to coast (to East) at current location
        theta_c <- (raster::extract(data$land.dir, rbind(X[i - 1, 1:2]), na.rm=TRUE) + 0.6 * pi) %% (2*pi)
        ## draw n proposal steps for move deflection (to avoid land)
        d_c <-
          rw(
            n = mpar$ntries,
            mu = theta_c,
            rho = 1,
            a = mpar$a,
            b = mpar$b
          ) * max(c((1 - (rd[i-1] - mpar$mindist) / (mpar$buffer[1] - mpar$mindist)), 0))
        d_s <- d_s * min(c((rd[i-1] - mpar$mindist) / (mpar$buffer[1] - mpar$mindist), 1))
        
      } 
      
      if(!is.na(mpar$buffer[2]) && rz[i-1] < mpar$buffer[2]){
        ## direction opposite to -900 m isobath at current location
        theta_b <- (raster::extract(data$b900.dir, rbind(X[i - 1, 1:2]), na.rm=TRUE) - 0.6 * pi) %% (2*pi)
        ## draw n proposal steps for move deflection to stay in water < 900 m deep (~ on shelf)
        d_b <- rw(n = mpar$ntries,
                  mu = theta_b,
                  rho = 0.9,
                  a = mpar$a,
                  b = mpar$b
        ) #* (1 - rz[i-1] / mpar$buffer[2])
        d_s <- d_s * 0
        
      }
      
      ## apply envt'l forcing
      ## add advection due to current
      u <- extract(data$u, rbind(X[i -  1, 1:2])) * 3.6
      v <- extract(data$v, rbind(X[i -  1, 1:2])) * 3.6
      tmp <-
        cbind(X[i - 1, 1] + d_s[, 1] + d_c[, 1] + d_b[, 1] + u, 
              X[i - 1, 2] + d_s[, 2] + d_c[, 2] + d_b[, 2] + v)

      ## check if proposed update is on land
      dist <- extract(data$land, tmp)
      if(all(is.na(dist))) {
        if(i >= round(0.25 * N)) {
        X[i, 3] <- -1
        break
        } else {
          stop("\ntrack stuck at a coastal boundary")
        }
      }
      ## select first proposal that is > mindist (km) from land 
      idx <- which(dist > mpar$mindist)
      if(length(idx) == 0) {
        if(i >= round(0.2 * N)) {
        X[i, 3] <- -2
        break
        } else {
          stop("\ntrack stuck at a boundary")
        }
      }
      idx <- sample(idx, 1)
      
      rd[i] <- dist[idx]
      X[i, 1:2] <- tmp[idx, ]
      
      ## determine distance from current location to 900 m isobath
      rz[i] <- extract(data$b900.dist, rbind(X[i, 1:2]))
      
      ## determine survivourship
      X[i, 3] <- rbinom(1, 1, mpar$surv^(1/24)) # rescale daily survival to hourly prob. 
      
      if(pb){
        setTxtProgressBar(tpb, i)
        if(i==N | X[i, 3] == 0) close(tpb)
      }
      
      if(X[i, 3] == 0) break # stop if dead
    }
    ## truncate matrix if death occurs before i = N
    if(length(na.omit(X[, 3])) < N) {
      end <- which(X[, 3] < 1)
      X <- X[1:end, ]
      N <- nrow(X)
    }
    X <- data.frame(X)
    names(X) <- c("x", "y", "s")

    sim <- X %>% as_tibble() %>%
      mutate(id = rep(id, N)) %>%
      mutate(date = seq(ISOdatetime(2018,07,09,00,00,00), by = 3600, length.out = N)) %>%
      select(id, date, x, y, s)
  
    param <- mpar  
    out <- list(sim = sim, params = param)  
    class(out) <- "simsmolt"
    
    return(out)
    
}
