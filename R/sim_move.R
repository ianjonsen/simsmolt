#' @title simulate movements & survival, using inputs from presim
#' 
#' @description simulates, movement & detection of acoustically-tagged animals
#' 
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#' 
#' @param N - number of 1-h time steps to simulate
#' @param tag - start location(s) of simulated animals
#' @param coa - optional Centre-Of-Attraction location(s) to provide movement bias(es)
#' @param mpar - movement & (otpional) survival parameters supplied as a list, see details
#' @param data - a list of required data from \code{presim}
#' @importFrom raster extract
#' @importFrom dplyr %>%
#' @importFrom tibble as_tibble
#' @importFrom stats runif rbinom
#' @export
#' 
sim_move <-
  function(N = 1800,
           tag = c(712, 761),
           data = NULL,
           mpar = list(),
           pb = TRUE
           ) {
    ## default move parameters
    mpar.full <- list(
      coa = c(300,1700) + rnorm(2, 0, 100),
      a = 2,
      b = 0.864,
      rho_s = 0.9,
      ntries = 100,
      surv = 0.9936,
      taxis = "no",
      buffer = c(20, 30),
      mindist = 5,
      maxdist = 820
    )
    pnms <- names(mpar.full)
    
    if (is.null(data))
      stop("Can't find output from sim_setup()\n")
    if (class(data$land)[1] != "RasterLayer") stop("distance2land must be a RasterLayer")
    if (!is.null(data$uv)) {
      if (class(data$uv)[1] != "RasterBrick" || (!names(data$uv) %in% c("u","v"))) 
        stop("current data must be supplied as a RasterBrick with u and v layers")
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
    
    ## define location/survivourship, s matrix & start position
    X <- matrix(NA, N, 3)
    s <- matrix(NA, N, 2)
    X[1,] <- cbind(tag[1], tag[2], 1)
    delta <- theta_s <- theta_d <- theta_z <- rd <- rz <- c()

    ## recursion
    for (i in 2:N) {
      if(i==2 && pb)  tpb <- txtProgressBar(min = 2, max = N, style = 3)
      
      ## get distance from land at first step
      if(i == 2) {
        rd[i-1] <- extract(data$land, rbind(X[i - 1, 1:2]))
        rz[i-1] <- extract(data$b700.dist, rbind(X[i - 1, 1:2]))
      }
      ## direction to centre of attraction from current location
      theta_s[i] <- atan2(mpar$coa[1] - X[i - 1, 1], mpar$coa[2] - X[i - 1, 2])
      
      ## distance to centre of attraction from current location
      delta[i] <- sqrt((mpar$coa[1] - X[i - 1, 1]) ^ 2 + (mpar$coa[2] - X[i - 1, 2]) ^ 2)
      
      switch(mpar$taxis, no = {
        ## biased RW, no taxis
        ## draw n proposal steps for active swimming
        s <-
          rw(
            n = mpar$ntries,
            mu = theta_s[i],
            rho = mpar$rho_s,
            a = mpar$a,
            b = mpar$b
          ) 
        
        ## apply weighting based on current distance to coast, 700 m isobath & specified buffer distance
        if(rz[i-1] > mpar$buffer[2]) {
#          cat(cbind(i, rd[i-1], rz[i-1]),"\n")
          ## direction ~ parallel to coast (to East) at current location
          theta_d <- (extract(data$land.dir, rbind(X[i - 1, 1:2]), buffer=10, fun=mean, na.rm=TRUE) + 0.6 * pi) %% (2*pi)
          ## draw n proposal steps for move deflection (to avoid land)
          d <-
            rw(
              n = mpar$ntries,
              mu = theta_d,
              rho = 1,
              a = mpar$a,
              b = mpar$b
            ) * max(c((1 - (rd[i-1] - mpar$mindist) / (mpar$buffer[1] - mpar$mindist)), 0))
          s <- s * min(c((rd[i-1] - mpar$mindist) / (mpar$buffer[1] - mpar$mindist), 1))
          z <- cbind(0,0)
          
        } else if(rz[i-1] <= mpar$buffer[2]){
#          cat(cbind(i, rd[i-1], rz[i-1]),"\n")
          ## direction opposite to -700 m isobath at current location
          theta_z <- (extract(data$b700.dir, rbind(X[i - 1, 1:2]), na.rm=TRUE) - 0.75 * pi) %% (2*pi)
          ## draw n proposal steps for move deflection to stay in water < 700 m deep (~ on shelf)
          z <- rw(n = mpar$ntries,
                  mu = theta_z,
                  rho = 1,
                  a = mpar$a,
                  b = mpar$b
          ) * (1 - rz[i-1] / mpar$buffer[2])
          s <- s * rz[i-1] / mpar$buffer[2]
        } 
        
        if (is.null(data$uv)) {
          ## no current advection
          tmp <- cbind(X[i - 1, 1] + d[, 1] + s[, 1] + z[, 1], X[i - 1, 2] + d[, 2] + s[, 2] + z[, 2])
          
        } else {
          ## add advection due to current, if raster is supplied
          uv <-
            extract(data$uv, rbind(X[i - 1, 1:2])) * 3600 / 1000 # to convert from m/s to km/h
          tmp <-
            cbind(X[i - 1, 1] + d[, 1] + s[, 1] + uv[1], X[i - 1, 2] + d[, 2] + s[, 2] + uv[2])
        }
      },
      
      rp = {
        ## positive rheotaxis (swim against current)
        if (is.null(data$uv)) stop("a current raster must be supplied to simulate rheotaxis")
        uv <- extract(data$uv, rbind(X[i - 1, 1:2]), method = "bilinear") * 3600 / 1000
        mu <- (atan2(uv[1], uv[2]) + pi) %% (2*pi) ## move dir is opposite current dir
        ## draw ntries proposal steps
        d <- rw(n = mpar$ntries, mu = mu, rho = 0.95, a = mpar$a, b = mpar$b)
        tmp <- cbind(X[i - 1, 1] + d[, 1] , X[i - 1, 2] + d[, 2] ) # FIXME: remove uv advection as a test
      },
      
      rn = {
        ## negative rheotaxis (swim with current)
        if (is.null(data$uv)) stop("a current raster must be supplied to simulate rheotaxis")
        uv <- extract(data$uv, rbind(X[i - 1, 1:2]), method = "bilinear") * 3600 / 1000
        mu <- atan2(uv[1], uv[2]) %% (2*pi) ## move dir is with current dir
        ## draw ntries proposal steps
        d <- rw(n = mpar$ntries, mu = mu, rho = 0.95, a = mpar$a, b = mpar$b)
        tmp <- cbind(X[i - 1, 1] + d[, 1] + uv[1], X[i - 1, 2] + d[, 2] + uv[2])
      })

      ## check if proposed update is on land
      dist <- extract(data$land, tmp)
      if(all(is.na(dist))) {
        stop("track stuck at a coastal boundary, try increasing 'mpar$ntries'")
      }
      ## select first proposal that is > mindist (km) from land 
      idx <- which(dist > mpar$mindist)
      if(length(idx)==0) {
        X[i, 3] <- -1
        break
        #stop("track stuck at a coastal boundary")
      }
      idx <- sample(idx, 1)
      
      rd[i] <- dist[idx]
      X[i, 1:2] <- tmp[idx, ]
      
      ## determine distance from current location to 700 m isobath
      rz[i] <- extract(data$b700.dist, rbind(X[i, 1:2]))
      
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
    }
    X <- data.frame(X)
    names(X) <- c("x", "y", "s")

    sim <- X %>% as_tibble()
    
    
    data$tag <- tag
    data$coa <- mpar$coa
    data$y.cpts <- mpar$y.cpts
    data$taxis <- mpar$taxis
  
    out <- list(sim = sim, data = data)  
    class(out) <- "simsmolt"
    
    return(out)
    
}
