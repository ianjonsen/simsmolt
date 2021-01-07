#' @title simulate movements with or without current advection, using inputs from sim_setup
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
#' @importFrom raster extract cellStats adjacent cellFromXY xyFromCell
#' @importFrom CircStats rwrpcauchy
#' @importFrom dplyr %>%
#' @importFrom tibble as_tibble
#' @importFrom stats runif rbinom
#' @export
#' 
sim_growth <-
  function(id=1, 
           N = 2760,
           data = NULL,
           mpar = list(),
           pb = TRUE
  ) {
    ## default move parameters
    
    mpar.full <- list(
      move = "brw",
      start = c(990, 1245),    #c(750, 200)
      coa = c(400, 2390),      #c(610, 394),
      rho = 0.8,
      ntries = 1,
      advect = TRUE,
      shelf = TRUE,
      growth = TRUE,
      buffer = c(5, 50),
      b = 1.6, ## assumed sustained travel speed of body-lengths / s
      a = 0, ## scale parameter of Weibull dist for move steps; bigger = less variable step lengths; 0 = no variability
      w0 = 47 ## starting mass g
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
    
    ## keep track of whether simulation hits a boundary or land (set to FALSE at start)
    mpar$boundary <- FALSE
    mpar$land <- FALSE
      
    
    ## define location matrix & initialise start position
    ## ds - active swimming displacements
    ## dl - displacements to deflect away from land
    xy <- matrix(NA, N, 2)
    xy[1,] <- cbind(mpar$start)       #cbind(sloc[1], sloc[2])
    ds <- matrix(NA, N, 4)
    ds[1,] <- c(NA, NA, 45/180*pi, NA)
    
    ## define other vectors
    u <- v <- vector("numeric", N)
    
    if(mpar$shelf) {
      rdb <- vector("numeric", N)
      rdb[1] <- extract(data$d2b900, rbind(mpar$start))
    }
    
    if(mpar$growth) {
      s <- ts <- w <- fl <- vector("numeric", N)
      w[1] <- mpar$w0
      fl[1] <- (w[1] / 8987.9) ^ (1 / 2.9639)
      s[1] <- fl[1] * mpar$b * 3.6 ## initial swim speed (fl * b body-lengths / s) - in km/h
    } else {
      fl <- (mpar$w0 / 8987.9) ^ (1 / 2.9639)
      s <- rep(fl * mpar$b * 3.6, N)
    }
    
    ## recursion
    for (i in 2:N) {
      if(i==2 && pb)  tpb <- txtProgressBar(min = 2, max = N, style = 3)
      
      ### Apply Energetics
      if(mpar$growth) {
      ## extract Temperature
      ts[i-1] <- extract(data$ts, rbind(xy[i - 1, ]))
      if(is.na(ts[i - 1])) ts[i - 1] <- ts[i - 2] #cellStats(data$ts, mean, na.rm = TRUE)
      
      ## calculate growth in current time step based on water temp at previous location, etc...
      w[i] <- growth(w[i-1], ts[i-1], s[i-1])
      ## convert W to fL - based on Byron et al 2014 (fig A1 in S1)
      fl[i] <- (w[i] / 8987.9) ^ (1 / 2.9639)
      
      ## determine size-based average move step for current time step
      ## assume avg swim speed of b bL/s
      s[i] <- fl[i] * mpar$b * 3.6 ## forklength * b m/s converted to km/h
      } 
      
      ### Current Advection
      if (mpar$advect) {
        ## determine envt'l forcing
        ## determine advection due to current, convert from m/s to km/h
        u[i] <- extract(data$u, rbind(xy[i -  1,])) * 3.6 
        v[i] <- extract(data$v, rbind(xy[i -  1,])) * 3.6
        u[i] <- ifelse(is.na(u[i]), 0, u[i])
        v[i] <- ifelse(is.na(v[i]), 0, v[i])
        
      } else {
        u[i] <- v[i] <- 0
      }

      ## Move Vectors
      ds[i, ] <- switch(mpar$move,
                        brw = {
                          biased_rw(n=1, 
                                    data, 
                                    xy = xy[i-1,], 
                                    coa = mpar$coa, 
                                    buffer = mpar$buffer, 
                                    rho = mpar$rho,
                                    a = mpar$a,
                                    b = s[i])
                        },
                        rw = {
                          rw(n=1,
                             data,
                             xy = xy[i-1,],
                             buffer = mpar$buffer,
                             a = mpar$a,
                             b = s[i])
                        },
                        crw = {
                          crw(n=1,
                              data,
                              xy = xy[i-1,],
                              buffer = mpar$buffer,
                              mu = ds[i-1, 3],
                              rho = mpar$rho,
                              a = mpar$a,
                              b = s[i])
                        })
      
      xy[i, 1:2] <- cbind(ds[i, 1] + u[i], 
                          ds[i, 2] + v[i])
      
      if(ds[i, 4] == 0) {
        mpar$land <- TRUE
        cat("\n stopping simulation: stuck on land")
        break
      }
      
      if(any(is.na(xy[i, ]))) {
        mpar$boundary <- TRUE
        cat("\n stopping simulation: hit a boundary")
        break
      } 
      
      if(pb){
        setTxtProgressBar(tpb, i)
        if(i==N) close(tpb)
      }
    }
    if(mpar$growth) {
    X <-
      data.frame(
        x = xy[, 1],
        y = xy[, 2],
        d2l = ds[, 4],
        u = u,
        v = v,
        ts = ts,
        w = w,
        fl = fl
      )
    } else if(!mpar$growth) {
      X <-
        data.frame(
          x = xy[, 1],
          y = xy[, 2],
          d2l = ds[, 4],
          u = u,
          v = v
        )
    }
    
    sim <- X %>% as_tibble() 
    sim <- sim %>%
      mutate(id = rep(id, nrow(sim))) %>%
      mutate(date = seq(ISOdatetime(2018,07,09,00,00,00), by = 3600, length.out = nrow(sim))) %>%
      select(id, date, everything())
    
    param <- mpar  
    out <- list(sim = sim, params = param)  
    class(out) <- "simsmolt"
    
    return(out)
    
  }
