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
           N = 1800,
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
      pdrf = c(3.017, -0.0139)
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
    
    if (is.null(mpar$coa)) {
      if (is.na(mpar$buffer[2])) {
        mpar$coa <-
          cbind(
            c(
              runif(1, 625, 800),
              runif(1, 480, 1200),
              runif(1, 240, 1200),
              runif(1, 120, 1200),
              runif(1, 10, 1200)
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
              runif(1, 625, 800),
              runif(1, 480, 760),
              runif(1, 240, 460),
              runif(1, 120, 290),
              runif(1, 10, 220)
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
    
    ## define location/survivourship, s matrix & start position
    X <- matrix(NA, N, 3)
    s <- matrix(NA, N, 2)
    X[1,] <- cbind(mpar$start[1], mpar$start[2], 1)
    theta_s <- rd <- rz <- c()
    theta_s[1] <- 0 / 180 * pi
    
    ## recursion
    for (i in 2:N) {
      if(i==2 && pb)  tpb <- txtProgressBar(min = 2, max = N, style = 3)
      
      ## get distance from land at first step
      if(i == 2) {
        rd[i-1] <- extract(data$land, rbind(X[i - 1, 1:2]))
        rz[i-1] <- extract(data$b900.dist, rbind(X[i - 1, 1:2]))
      }
      ## direction to centre of attraction from current location
      if(all(!is.na(mpar$coa))) {
        ## get dist & dir to all coa's
        d2coa <- sqrt((mpar$coa[,1] - X[i - 1, 1])^2 + (mpar$coa[,2] - X[i - 1, 2])^2)
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
      
      switch(mpar$taxis, no = {
        ## biased RW, no taxis
        ## draw n proposal steps for active swimming
        s <-
          rw(
            n = mpar$ntries,
            mu = theta_s[i],
            rho = mpar$rho,
            a = mpar$a,
            b = mpar$b
          ) 
        z <- cbind(0,0)
        d <- cbind(0,0)
        
        ## apply weighting based on current distance to coast, 700 m isobath & specified buffer distance
        ##  first buffer's influence starts after exiting SoBI and after influence of first coa
         if(rd[i-1] <= mpar$buffer[1] & X[i - 1, 2] > 850) {
          ## direction ~ parallel to coast (to East) at current location
          theta_d <- (extract(data$land.dir, rbind(X[i - 1, 1:2]), na.rm=TRUE) + 0.6 * pi) %% (2*pi)
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

         } 

        if(!is.na(mpar$buffer[2]) && rz[i-1] < mpar$buffer[2]){
          ## direction opposite to -700 m isobath at current location
          theta_z <- (extract(data$b900.dir, rbind(X[i - 1, 1:2]), na.rm=TRUE) - 0.6 * pi) %% (2*pi)
          ## draw n proposal steps for move deflection to stay in water < 700 m deep (~ on shelf)
          z <- rw(n = mpar$ntries,
                  mu = theta_z,
                  rho = 0.9,
                  a = mpar$a,
                  b = mpar$b
          ) #* (1 - rz[i-1] / mpar$buffer[2])
          s <- s * 0
          
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
        stop("\ntrack stuck at a coastal boundary")
       # X[i, 3] <- -1
       # break
      }
      ## select first proposal that is > mindist (km) from land 
      idx <- which(dist > mpar$mindist)
      if(length(idx)==0) {
        stop("\ntrack stuck at a coastal boundary")
        #X[i, 3] <- -1
        #break
        #stop("track stuck at a coastal boundary")
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
      mutate(id = id) %>%
      mutate(date = seq(ISOdatetime(2018,07,09,00,00,00), by = 3600, length.out = N)) %>%
      select(id, date, x, y, s)
  
    param <- mpar  
    out <- list(sim = sim, params = param)  
    class(out) <- "simsmolt"
    
    return(out)
    
}
