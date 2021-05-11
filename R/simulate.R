#' @title simulate smolt migration 
#' 
#' @description simulates smolt migration with or without: 1) temperature-dependent growth; 2) temperature-dependent migration; 3) current advection, using inputs from sim_setup. Movements are simulated as either a biased random walk (brw) or a random walk (rw)
#' 
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#' 
#' @param id - identifier for simulation run (individual animal)
#' @param N - number of 1-h time steps to simulate
#' @param data - a list of required data from \code{presim}
#' @param mpar - simulation control parameters supplied as a list, see details
#' @param pb - use progress bar (logical)
#' @importFrom raster extract xyFromCell
#' @importFrom CircStats rwrpcauchy
#' @importFrom dplyr %>% mutate lag
#' @importFrom tibble as_tibble
#' @importFrom stats runif rbinom
#' @importFrom lubridate week yday
#' @importFrom stringr str_split
#' @export
#' 
simulate <-
  function(id=1, 
           N = 2760,
           data = NULL,
           mpar = sim_par(),
           pb = TRUE
  ) {
    
    if (is.null(data))
      stop("Can't find output from sim_setup()\n")
    if (class(data$land)[1] != "RasterLayer") stop("distance2land must be a RasterLayer")
    
    if(mpar$advect & !all(c("u","v") %in% names(data))) {
      cat("Turning off current-advected movement as data do not contain currents\n")
      mpar$advect <- FALSE
    }
      
    ## define location matrix & initialise start position
    ## ds - active swimming displacements
    ## dl - displacements to deflect away from land
    xy <- matrix(NA, N, 2)
    xy[1,] <- cbind(mpar$pars$start)       #cbind(sloc[1], sloc[2])
    ds <- matrix(NA, N, 2)
    ds[1,] <- c(NA, NA)
    
    ## define other vectors
    u <- v <- vector("numeric", N)
    
    if(mpar$growth) {
      s <- ts <- w <- fl <- vector("numeric", N)
     # move_dir <- rep(NA, N)
      w[1] <- mpar$pars$w0
      fl[1] <- (w[1] / 8987.9) ^ (1 / 2.9639)
      s[1] <- fl[1] * mpar$pars$b * 3.6 ## initial swim speed (fl * b body-lengths / s) - in km/h
    } else {
      fl <- (mpar$pars$w0 / 8987.9) ^ (1 / 2.9639)
      s <- rep(fl * mpar$pars$b * 3.6, N)
    }
    
    ## what is start week in env data
    if(data$ocean == "doy") {
      d1 <- as.numeric(str_split(names(data$ts)[1], "d", simplify = TRUE)[,2]) - 1
    }
      
    ## iterate movement
    for (i in 2:N) {
      if(i==2 && pb)  tpb <- txtProgressBar(min = 2, max = N, style = 3)
      
      ### Apply Energetics
      if(mpar$growth) {
      ## extract Temperature
        switch(data$ocean, 
               cl = {
                 ts[i-1] <- extract(data$ts, rbind(xy[i - 1, ])) - 273
                 if(is.na(ts[i-1])){
                   ## calc mean Temp within 2 km buffer of location @ time i-1
                   ts[i-1] <- extract(data$ts, rbind(xy[i - 1, ]), buffer = 2, df = TRUE)[,2] %>%
                     mean(., na.rm = TRUE) - 273
                 }
               },
               doy = {
                 ts[i-1] <- extract(data$ts[[(yday(mpar$pars$start.dt + i * 3600) - d1)]], rbind(xy[i - 1, ])) - 273
                 if(is.na(ts[i-1])) {
                   ## calc mean Temp within 2 km buffer of location @ time i-1
                   ts[i-1] <- extract(data$ts[[(yday(mpar$pars$start.dt + i * 3600) - d1)]], rbind(xy[i - 1, ]), buffer = 2, df = TRUE)[,2] %>%
                     mean(., na.rm = TRUE) - 273
                 }
               })
      
      ## calculate growth in current time step based on water temp at previous location, etc...
      if(ts[i-1] <= 0 & !is.na(ts[i-1])) {
        cat("\n stopping simulation: smolt has entered water <= 0 deg C")  
        break
      } else if(ts[i-1] > 0 & !is.na(ts[i-1])) {
        w[i] <- growth(w[i-1], ts[i-1], s[i-1])
        
      } else if(is.na(ts[i-1])) {
        cat("\n stopping simulation: NA value for temperature")
        break
      }
        
      ## convert W to fL - based on Byron et al 2014 (fig A1 in S1)
      fl[i] <- (w[i] / 8987.9) ^ (1 / 2.9639)
      ## smolts can't shrink their length... (this would unrealistically affect swim speed & energetics), so
      ## if change in weight implies reduction in forklength, then stick with last forklength  - mimics loss of condition
      if(fl[i] < fl[i-1]) fl[i] <- fl[i-1] 
      
      ## determine size-based average move step for current time step
      ## assume avg swim speed of b bL/s
      s[i] <- fl[i] * mpar$pars$b * 3.6 ## forklength * b m/s converted to km/h
      } 
      
      ### Current Advection
      if (mpar$advect) {
        ## determine envt'l forcing
        ## determine advection due to current, convert from m/s to km/h
        if(data$ocean == "doy") {
          u[i] <- extract(data$u[[(yday(mpar$pars$start.dt + i * 3600) - d1)]], rbind(xy[i - 1, ])) * 3.6 
          v[i] <- extract(data$v[[(yday(mpar$pars$start.dt + i * 3600) - d1)]], rbind(xy[i - 1, ])) * 3.6
          } else if(data$ocean == "cl") {
          u[i] <- extract(data$u, rbind(xy[i - 1, ])) * 3.6 
          v[i] <- extract(data$v, rbind(xy[i - 1, ])) * 3.6
          }
        
        # downscale advection effect over first 21 d then increase to 1 over 7 d
        u[i] <- ifelse(is.na(u[i]), 0, u[i]) * ifelse(i < 500, mpar$pars$psi, ifelse(i < 668, -2 + i * 0.004491, 1))
        v[i] <- ifelse(is.na(v[i]), 0, v[i]) * ifelse(i < 500, mpar$pars$psi, ifelse(i < 668, -2 + i * 0.004491, 1))
        
      } else if(!mpar$advect | all(xy[1] >= data$sobi.box[1], 
                                   xy[1] <= data$sobi.box[2], 
                                   xy[2] >= data$sobi.box[3], 
                                   xy[2] <= data$sobi.box[4])) {
        u[i] <- v[i] <- 0
      }

      ### Temperature-dependent movement
      if (mpar$temp) {
        
        ## reverse migration bias from mpar$pars$mdir to mpar$pars$mdir - pi, if smolt in < 5 C water for 12 h
        if(i > 12) { 
          
          ## movement direction influenced by ts spatial gradient of current ts 
          ##  implies 0 or -ve growth, for current mass(w[i]) and speed (s[i])
          g.rng <- growth(w[i], seq(6, 20, l = 100), s[i])
          ts.mig <- seq(6, 20, l=100)[which(g.rng >= w[i])] %>% min()
          ts.rng <- seq(6, 20, l = 100)[which(g.rng >= quantile(g.rng, mpar$pars$ts.q))]  %>% range()
          
          ## if smolt in optimal T range for growth then switch from biased RW to simple RW
          dir <- ifelse(all(ts[i - 1:12] <= ts.mig), mpar$pars$mdir[2], mpar$pars$mdir[1])
          move <- ifelse(ts[i-1] >= ts.rng[1] & ts[i-1] <= ts.rng[2], "rw", mpar$move)
            
        } else {
          dir <- mpar$pars$mdir[1]
          move <- mpar$move
        }
      
      } else if(!mpar$temp) {
        dir <- mpar$pars$mdir[1]
        move <- mpar$move
      }
        
        ## Movement
        ## First check if smolt is in SoBI & within 25 km of land, if so then move toward Lab Sea after which movement rules can be applied
      d2l <- extract(data$land, rbind(xy[i-1,]))
        if(all(xy[i-1,1] >= data$sobi.box[1], 
               xy[i-1,1] <= data$sobi.box[2], 
               xy[i-1,2] >= data$sobi.box[3], 
               xy[i-1,2] <= data$sobi.box[4]) & d2l < 15) {

          phi <- rwrpcauchy(1, ifelse(mpar$pars$psi< 0.5, 0.28, 0.35) * pi, 0.8)
          ds[i, 1] <- xy[i-1, 1] + s[i] * sin(phi)
          ds[i, 2] <- xy[i-1, 2] + s[i] * cos(phi)
          
        } else {
          
        ds[i, ] <- switch(move,
                          brw = {
                            
                            brw(n=1, 
                                      data, 
                                      xy = xy[i-1,], 
                                      coa = mpar$pars$coa, 
                                      dir = dir,
                                      buffer = mpar$pars$buffer, 
                                      rho = mpar$pars$rho[1],
                                      a = mpar$pars$a,
                                      b = s[i], 
                                      taxis = mpar$taxis,
                                      u = u[i],
                                      v = v[i], 
                                      shelf = mpar$shelf,
                                      beta = mpar$pars$beta)
                          },
                          rw = {

                            rw(n=1,
                                data,
                                xy = xy[i-1,],
                                buffer = mpar$pars$buffer,
                                rho = mpar$pars$rho[2],
                                a = mpar$pars$a,
                                b = s[i],
                                taxis = mpar$taxis,
                                u = u[i],
                                v = v[i],
                                shelf = mpar$shelf,
                                beta = mpar$pars$beta)
                          },
                          drift = {
                            ds[i, ] <- rbind(xy[i-1, 1:2])
                          })
        }
      
      xy[i, 1:2] <- cbind(ds[i, 1] + u[i], 
                          ds[i, 2] + v[i])

      
      if((extract(data$land, rbind(xy[i, ])) == 0 | is.na(extract(data$land, rbind(xy[i, ]))))  & any(!is.na(xy[i,]))) {
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
        dx = ds[, 1] - lag(xy[, 1]),
        dy = ds[, 2] - lag(xy[, 2]),
        u = u,
        v = v,
        ts = ts,
        w = w,
        fl = fl,
        s = s
      )
    } else if(!mpar$growth) {
      X <-
        data.frame(
          x = xy[, 1],
          y = xy[, 2],
          dx = ds[, 1] - lag(xy[, 1]),
          dy = ds[, 2] - lag(xy[, 2]),
          u = u,
          v = v, 
          ts = ifelse(mpar$temp, ts, NA)
        ) 
    } 
    if(sum(is.na(X$ts)) == nrow(X)) {
      X <- X %>% select(-ts)
    } else {
      X$ts[nrow(X)] <- ifelse(X$ts[nrow(X)] == 0, NA, X$ts[nrow(X)])
    }
    
    sim <- X %>% as_tibble() 
    
    ## remove records after sim is stopped for being stuck on land, etc...
    if(mpar$land | mpar$boundary) {
      if(mpar$growth) {
        sim <- sim %>%
          filter(!is.na(x) & !is.na(y) & w != 0 & fl != 0 & s != 0)
      } else {
        sim <- sim %>%
          filter(!is.na(x) & !is.na(y))
      }
    }
    
    nsim <- nrow(sim)

    sim <- sim %>%
      mutate(id = id) %>%
      mutate(date = seq(mpar$pars$start.dt, by = 3600, length.out = nsim)) %>%
      select(id, date, everything())
    
    param <- mpar  
    out <- list(sim = sim, params = param)  
    class(out) <- "simsmolt"
    
    return(out)
    
  }
