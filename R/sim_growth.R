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
#' @importFrom raster extract xyFromCell
#' @importFrom CircStats rwrpcauchy
#' @importFrom dplyr %>% mutate
#' @importFrom tibble as_tibble
#' @importFrom stats runif rbinom
#' @importFrom lubridate week yday
#' @importFrom stringr str_split
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
      start.dt = ISOdatetime(2018,07,09,00,00,00, tz="UTC"),
      start = c(990, 1245),    #c(750, 200)
      coa = NULL, #c(400, 2390),      #c(610, 394),
      mdir = -20/180*pi,
      rho = 0.8,
      ntries = 1,
      temp = TRUE,
      ts.q = 0.75,
      advect = TRUE,
      shelf = TRUE,
      growth = TRUE,
      buffer = 5,
      b = 1.6, ## assumed sustained travel speed of body-lengths / s
      a = 0, ## scale parameter of Weibull dist for move steps; bigger = less variable step lengths; 0 = no variability
      w0 = 185 ## starting mass g
    )
    pnms <- names(mpar.full)
    
    if (is.null(data))
      stop("Can't find output from sim_setup()\n")
    if (class(data$land)[1] != "RasterLayer") stop("distance2land must be a RasterLayer")

    if(mpar$advect & !all(c("u","v") %in% names(data))) {
      cat("Turning off current-advected movement as ocean data do not contain currents\n")
      mpar$advect <- FALSE
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
    
    ## keep track of whether simulation hits a boundary or land (set to FALSE at start)
    mpar$boundary <- FALSE
    mpar$land <- FALSE
      
    
    ## define location matrix & initialise start position
    ## ds - active swimming displacements
    ## dl - displacements to deflect away from land
    xy <- matrix(NA, N, 2)
    xy[1,] <- cbind(mpar$start)       #cbind(sloc[1], sloc[2])
    ds <- matrix(NA, N, 2)
    ds[1,] <- c(NA, NA)
    
    ## define other vectors
    u <- v <- vector("numeric", N)
    
    if(mpar$growth) {
      s <- ts <- w <- fl <- vector("numeric", N)
     # move_dir <- rep(NA, N)
      w[1] <- mpar$w0
      fl[1] <- (w[1] / 8987.9) ^ (1 / 2.9639)
      s[1] <- fl[1] * mpar$b * 3.6 ## initial swim speed (fl * b body-lengths / s) - in km/h
    } else {
      fl <- (mpar$w0 / 8987.9) ^ (1 / 2.9639)
      s <- rep(fl * mpar$b * 3.6, N)
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
                 ts[i-1] <- extract(data$ts[[(yday(mpar$start.dt + i * 3600) - d1)]], rbind(xy[i - 1, ])) - 273
                 if(is.na(ts[i-1])) {
                   ## calc mean Temp within 2 km buffer of location @ time i-1
                   ts[i-1] <- extract(data$ts[[(yday(mpar$start.dt + i * 3600) - d1)]], rbind(xy[i - 1, ]), buffer = 2, df = TRUE)[,2] %>%
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
        browser()
      }
        
      ## convert W to fL - based on Byron et al 2014 (fig A1 in S1)
      fl[i] <- (w[i] / 8987.9) ^ (1 / 2.9639)
      ## smolts can't shrink their length... (this would unrealistically affect swim speed & energetics), so
      ## if change in weight implies reduction in forklength, then stick with last forklength  - mimics loss of condition
      if(fl[i] < fl[i-1]) fl[i] <- fl[i-1] 
      
      ## determine size-based average move step for current time step
      ## assume avg swim speed of b bL/s
      s[i] <- fl[i] * mpar$b * 3.6 ## forklength * b m/s converted to km/h
      } 
      
      ### Current Advection
      if (mpar$advect) {
        ## determine envt'l forcing
        ## determine advection due to current, convert from m/s to km/h
        if(data$ocean == "doy") {
          u[i] <- extract(data$u[[(yday(mpar$start.dt + i * 3600) - d1)]], rbind(xy[i - 1, ])) * 3.6 
          v[i] <- extract(data$v[[(yday(mpar$start.dt + i * 3600) - d1)]], rbind(xy[i - 1, ])) * 3.6
          } else if(data$ocean == "cl") {
          u[i] <- extract(data$u, rbind(xy[i - 1, ])) * 3.6 
          v[i] <- extract(data$v, rbind(xy[i - 1, ])) * 3.6
          }
        
        u[i] <- ifelse(is.na(u[i]), 0, u[i])
        v[i] <- ifelse(is.na(v[i]), 0, v[i])
        
      } else if(!mpar$advect) {
        u[i] <- v[i] <- 0
      }

      ### Temperature-dependent movement
      if (mpar$temp) {
        
        ## reverse migration bias from -20 deg (~magnetic N) to 160 deg (S), if smolt in < 5 C water for 12 h
        if(i > 12) { 
          
          ## movement direction influenced by ts spatial gradient of current ts 
          ##  implies 0 or -ve growth, for current mass(w[i]) and speed (s[i])
          g.rng <- growth(w[i], seq(6, 20, l = 100), s[i])
          ts.mig <- seq(6, 20, l=100)[which(g.rng >= w[i])] %>% min()
          ts.rng <- seq(6, 20, l = 100)[which(g.rng >= quantile(g.rng, mpar$ts.q))]  %>% range()
          
          dir <- ifelse(all(ts[i - 1:12] <= ts.mig), (mpar$mdir + pi) %% (pi), mpar$mdir)
          move <- ifelse(ts[i-1] >= ts.rng[1] & ts[i-1] <= ts.rng[2], "rw", mpar$move)
            
        } else {
          dir <- mpar$mdir
          move <- mpar$move
        }
      
      #   ds[i, ] <- temp_brw(
      #       n = 1,
      #       i = i,
      #       mpar = mpar,
      #       d1 = d1,
      #       data = data,
      #       xy = xy[i - 1,],
      #       ts = ts[i-1],
      #       ts.rng = ts.rng,
      #       b = s[i]
      #     )
      # } else if(!mpar$temp) {
      } else {
        dir <- mpar$mdir
        move <- mpar$move
      }
        
        ## Temperature-independent movement
        ds[i, ] <- switch(move,
                          brw = {
                            
                            biased_rw(n=1, 
                                      data, 
                                      xy = xy[i-1,], 
                                      coa = NULL, 
                                      dir = dir,
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
                                rho = mpar$rho,
                                a = mpar$a,
                                b = s[i])
                          },
                          drift = {
                            ds[i, ] <- rbind(xy[i-1, 1:2])
                          })
 #     }
        
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
          u = u,
          v = v, 
          ts = ts
        )
    } 
    X$ts[nrow(X)] <- X$ts[nrow(X) - 1]
    
    sim <- X %>% as_tibble() 
    ## remove records after sim is stopped for being stuck on land, etc...
    sim <- sim %>%
      filter(!is.na(x) & !is.na(y) & w != 0 & fl != 0 & s != 0)
    nsim <- nrow(sim)
    
    sim <- sim %>%
      mutate(id = id) %>%
      mutate(date = seq(mpar$start.dt, by = 3600, length.out = nsim)) %>%
      select(id, date, everything())
    
    param <- mpar  
    out <- list(sim = sim, params = param)  
    class(out) <- "simsmolt"
    
    return(out)
    
  }
