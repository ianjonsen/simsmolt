#' @title simulate movements & survival, using inputs from presim
#' 
#' @description simulates, movement & detection of acoustically-tagged animals
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
sim_move <-
  function(id=1, 
           N = 2760,
           data = NULL,
           mpar = list(),
           pb = TRUE
           ) {
    ## default move parameters
    
    mpar.full <- list(
      start = c(560, 82.5),
      coa = NULL,
      a = 2,
      b = 0.864,
      rho = 0.8,
      ntries = 100,
      surv = 0.9936,
      taxis = "no",
      buffer = c(10, 10),
      mindist = 8,
      maxdist = 850, 
      pdrf = c(3.017, -0.0139) #c(4.865, -0.0139) = p(0.5) @ 350 m (~ consistent w HFX line V9 @ high power) 
    )
    pnms <- names(mpar.full)
    
    if (is.null(data))
      stop("Can't find output from sim_setup()\n")
    if (class(data$land)[1] != "RasterLayer") stop("distance2land must be a RasterLayer")
    if (all(!is.null(data$u), !is.null(data$v))) {
      if (class(data$u)[1] != "RasterBrick" || (!names(data$u) %in% c("jul","aug","sep","oct","nov"))) 
        stop("zonal current (u) data must be supplied as a RasterBrick with `jul` to `nov` layers")
      if (class(data$v)[1] != "RasterBrick" || (!names(data$v) %in% c("jul","aug","sep","oct","nov"))) 
        stop("meridional current (v) data must be supplied as a RasterBrick with `jul` to `nov` layers")
    }
    
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
    
    if (is.null(mpar$coa)) {
      p <- runif(1)
      if (is.na(mpar$buffer[2])) {
        mpar$coa <-
          cbind(
            c(
              625 + 275 * p,
              480 + 520 * p,
              240 + 760 * p,
              120 + 880 * p,
              10  + 990 * p
            ),
            c(
              runif(1, 155, 195),
              runif(1, 400, 440),
              runif(1, 555, 595),
              runif(1, 805, 845),
              runif(1, 1155, 1195)
            )
          )
      } else if (!is.na(mpar$buffer[2])) {
        mpar$coa <-
          cbind(
            c(
              625 + 175 * p,
              480 + 280 * p,
              240 + 220 * p,
              120 + 170 * p,
              10  + 210 * p
            ),
            c(
              runif(1, 155, 195),
              runif(1, 400, 440),
              runif(1, 555, 595),
              runif(1, 805, 845),
              runif(1, 1155, 1195)
            )
          )
      }
    }
    
    ## define location/survivourship, sw matrix & start position
    X <- matrix(NA, N, 3)
    #sw <- matrix(NA, N, 2)
    X[1,] <- cbind(mpar$start[1], mpar$start[2], 1)
    theta_s <- rd <- rz <- c()
    theta_s[1] <- 0 / 180 * pi

    ## recursion
    for (i in 2:N) {
      if(i==2 && pb)  tpb <- txtProgressBar(min = 2, max = N, style = 3)
      
      ## determine month to sample from appropriate u,v raster layer (assumes 9 July start date)
      m.i <- ifelse(i < 552, 1, ifelse(i >=552 && i < 1296, 2, ifelse(i >= 1296 && i < 2016, 3, 4)))
      
      ## get distance from land at first step
      if(i == 2) {
        rd[i-1] <- extract(data$land, rbind(X[i - 1, 1:2]))
        rz[i-1] <- extract(data$b900.dist, rbind(X[i - 1, 1:2]))
      }
      ## direction to centre of attraction from current location
      if(all(!is.na(mpar$coa))) {
        ## get dist & dir to all coa's
        d2coa <- sqrt((mpar$coa[,1] - X[i - 1, 1])^2 + (mpar$coa[,2] - X[i - 1, 2])^2)
        if(i == 2) mdist <- c(d2coa[1], diff(d2coa))
        theta <- atan2(mpar$coa[, 1] - X[i - 1, 1], mpar$coa[, 2] - X[i - 1, 2])
      
        
        ## calculate weighted swimming direction
        if(X[i - 1, 2] <= mpar$coa[1,2]) {
          if(i == 2) mdist <- d2coa[1]
          wt <- d2coa[1]/mdist
          theta_s[i] <- wt * theta[1] + (1 - wt) * theta[2]
          j <- i

        } else if(X[i - 1, 2] > mpar$coa[1,2] & X[i - 1, 2] <= mpar$coa[2,2]) {
          if(i == j+1) mdist <- d2coa[2]
          wt <- d2coa[2]/mdist
          theta_s[i] <- wt * theta[2] + (1 - wt) * theta[3]
          k <- i

        }
        else if(X[i - 1, 2] > mpar$coa[2,2] & X[i - 1, 2] <= mpar$coa[3,2]) {
          if(i == k+1) mdist <- d2coa[3]
          wt <- d2coa[3]/mdist
          theta_s[i] <- wt * theta[3] + (1 - wt) * theta[4]
          q <- i

        } else if(X[i - 1, 2] > mpar$coa[3,2] & X[i - 1, 2] <= mpar$coa[4,2]) {
          if(i == q+1) mdist <- d2coa[4]
          wt <- d2coa[4]/mdist
          theta_s[i] <- wt * theta[4] + (1 - wt) * theta[5]
          h <- i

        } else if(X[i - 1, 2] > mpar$coa[4,2] & X[i - 1, 2] <= mpar$coa[5,2]) {
          if(i == h+1) mdist <- d2coa[5]
          wt <- d2coa[5]/mdist
          theta_s[i] <- wt * theta[5] + (1 - wt) * 0/180*pi

        } else if(X[i - 1, 2] > mpar$coa[5,2]) {
          theta_s[i] <- 15/180*pi
        }

      } else {
        ## do this if no coa exists
        theta_s[i] <- rwrpcauchy(1, theta_s[i-1], mpar$rho) %% (2*pi)
      }
    
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
      
      ## apply weighting based on current distance to coast, 700 m isobath & specified buffer distance
      ##  first buffer's influence starts after exiting SoBI and after influence of first coa
      if(rd[i-1] <= mpar$buffer[1] & X[i - 1, 2] > 850) {
        ## direction ~ parallel to coast (to East) at current location
        theta_c <- (extract(data$land.dir, rbind(X[i - 1, 1:2]), na.rm=TRUE) + 0.6 * pi) %% (2*pi)
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
        ## direction opposite to -700 m isobath at current location
        theta_b <- (extract(data$b900.dir, rbind(X[i - 1, 1:2]), na.rm=TRUE) - 0.6 * pi) %% (2*pi)
        ## draw n proposal steps for move deflection to stay in water < 700 m deep (~ on shelf)
        d_b <- rw(n = mpar$ntries,
                  mu = theta_b,
                  rho = 0.9,
                  a = mpar$a,
                  b = mpar$b
        ) #* (1 - rz[i-1] / mpar$buffer[2])
        d_s <- d_s * 0
        
      }
      
      ## apply behavioural rules - taxis
      switch(mpar$taxis, no = {
        ## biased RW, no taxis
        
        if (is.null(data$u) || is.null(data$v)) {
          ## no current advection 
          tmp <- cbind(X[i - 1, 1] + d_s[, 1] + d_c[, 1] + d_b[, 1], 
                       X[i - 1, 2] + d_s[, 2] + d_c[, 2] + d_b[, 2])
          
        } else {
          ## add advection due to current, if raster is supplied
         
          u <- extract(data$u[[m.i]], rbind(X[i - 1, 1:2])) * 3.6 # to convert from m/s to km/h
          v <- extract(data$v[[m.i]], rbind(X[i - 1, 1:2])) * 3.6
          tmp <-
            cbind(X[i - 1, 1] + d_s[, 1] + d_c[, 1] + d_b[, 1] + u, 
                  X[i - 1, 2] + d_s[, 2] + d_c[, 2] + d_b[, 2] + v)
        }
      },
      
      rp = {
        ## positive rheotaxis (swim against current)
        if (is.null(data$u) || is.null(data$v)) 
          stop("u & v current rasters must be supplied to simulate rheotaxis")
        u <- extract(data$u[[m.i]], rbind(X[i - 1, 1:2])) * 3.6
        v <- extract(data$v[[m.i]], rbind(X[i - 1, 1:2])) * 3.6
        mu <- (atan2(u, v) + pi) %% (2*pi) ## move dir is opposite current dir
        ## draw ntries proposal steps, after 15 d
        if(X[i - 1, 2] > 135)
        d_s <- rw(n = mpar$ntries, mu = mu, rho = 0.95, a = mpar$a, b = mpar$b)
        tmp <- cbind(X[i - 1, 1] + d_s[, 1] + d_c[, 1] + d_b[, 1] + u, 
                     X[i - 1, 2] + d_s[, 2] + d_c[, 2] + d_b[, 2] + v) 
      },
      
      rn = {
        ## negative rheotaxis (swim with current)
        if (is.null(data$u) || is.null(data$v)) 
          stop("u & v current rasters must be supplied to simulate rheotaxis")
        u <- extract(data$u[[m.i]], rbind(X[i - 1, 1:2])) * 3.6
        v <- extract(data$v[[m.i]], rbind(X[i - 1, 1:2])) * 3.6
        mu <- atan2(u, v) %% (2*pi) ## move dir is with current dir
        ## draw ntries proposal steps, after 15 d
        if(X[i - 1, 2] > 135)
        d_s <- rw(n = mpar$ntries, mu = mu, rho = 0.95, a = mpar$a, b = mpar$b)
        tmp <- cbind(X[i - 1, 1] + d_s[, 1] + d_c[, 1] + d_b[, 1] + u, 
                     X[i - 1, 2] + d_s[, 2] + d_c[, 2] + d_b[, 2] + v)
      })

      ## check if proposed update is on land
      dist <- extract(data$land, tmp)
      if(all(is.na(dist))) {
        if(taxis != "no" && i >= round(0.25 * N)) {
        X[i, 3] <- -1
        break
        } else {
          stop("\ntrack stuck at a coastal boundary")
        }
      }
      ## select first proposal that is > mindist (km) from land 
      idx <- which(dist > mpar$mindist)
      if(length(idx) == 0) {
        if(taxis != "no" && i >= round(0.25 * N)) {
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
