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
#' @importFrom raster extract
#' @importFrom CircStats rwrpcauchy
#' @importFrom dplyr %>%
#' @importFrom tibble as_tibble
#' @importFrom stats runif rbinom
#' @export
#' 
sim_new <-
  function(id=1, 
           N = 2760,
           data = NULL,
           mpar = list(),
           pb = TRUE
  ) {
    ## default move parameters
    
    mpar.full <- list(
      start = NULL,
      rho = 0.8,
      ntries = 1,
      taxis = "no",
      a = 100, ## scale parameter of Weibull dist for move steps; bigger = less variable step lengths
      g = 0.01375,
      sm = 47
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
    
    
    
    if(is.null(mpar$start)) {
      sloc <- c(560, 82.5) #c(750, 200)
    } else {
      sloc <- mpar$start
    }
      
    
    ## define location/survivourship, sw matrix & start position
    X <- matrix(NA, N, 13)
    X[1,] <- cbind(sloc[1], sloc[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
    theta_s <- c()
    theta_s[1] <- 0 / 180 * pi
    
    
    ## recursion
    for (i in 2:N) {
      if(i==2 && pb)  tpb <- txtProgressBar(min = 2, max = N, style = 3)
      
      ## define starting mass (mpar$sm) as 47
      ## use 47.35 g as avg mass of smolts at tagging - based on 17 cm forklength and L - W power relationship of fL = (m/8987.9)^(1/2.9639)
      ## from Byron et al. 2014 Fish Oceanogr
      
      ## apply constant (for now) daily growth at rate mpar$g & convert to forklength
      X[i-1, 11] <- ((mpar$sm + i * mpar$g) / 8987.9)^(1/2.9639)
      
      ## determine hourly average move step (Weibull 'b' parameter)
      ## assume avg swim speed of 1.6 bL/s
      X[i-1, 12] <- X[i-1, 11] * 1.6 * 3.6 ## forklength * 1.6 m/s converted to km/h
      
      ## heading to move in
      #theta_s[i] <- rwrpcauchy(1, theta_s[i-1], mpar$rho) %% (2*pi)
      
      ## apply envt'l forcing
      ## determine advection due to current
      u <- extract(data$u, rbind(X[i -  1, 1:2])) * 3.6
      v <- extract(data$v, rbind(X[i -  1, 1:2])) * 3.6
      u <- ifelse(is.na(u), 0, u)
      v <- ifelse(is.na(v), 0, v)
      
      ## employ rheotaxis or not by adjusting move direction (default is a bias to North)
      ## the strength of adherence to taxis is dictated by rho: rho = 1 - perfect taxis, rho = 0 - random movement, no taxis
      ## magnitude of current relative to magnitude of swimming, normalized to 1
      X[i-1, 9] <- rel_mag <- ifelse(sqrt(u^2+v^2) / X[i-1, 12] > 1, 1, sqrt(u^2+v^2) / X[i-1, 12])
      
      if(mpar$taxis == "pos") {
        ## swim against current
        X[i-1, 10] <- mu_c <- (atan2(u, v) + pi)
        ## weighted direction
        mu <- (mu_c * rel_mag) %% (2*pi)
        
      } else if(mpar$taxis == "neg") {
        ## swim with current
        X[i-1, 10] <- mu_c <- atan2(u, v)
        mu <- (mu_c * rel_mag) %% (2*pi)
        
      } else if(mpar$taxis == "no") {
        ## biased movement to North with advection
        X[i-1, 10] <- mu_c <- NA
        mu <- 0
      }
        
      
      
      ## move vectors
      ## draw n proposal steps for active swimming
      d_s <-
        rw(
          n = mpar$ntries,
         # mu = theta_s[i],
          mu = mu, 
          rho = mpar$rho,
          a = mpar$a,
          b = X[i-1, 12]
        ) 
      
      
      X[i, 1:2] <-
        cbind(X[i - 1, 1] + d_s[, 1] + u, 
              X[i - 1, 2] + d_s[, 2] + v)
      
      ## keep track of displacement & advection 
      X[i-1, 3:4] <- cbind(d_s[,1], d_s[,2])
      X[i-1, 5:6] <- cbind(u,v)
      X[i-1, 7:8] <- cbind(d_s[, 3], d_s[,4]) # record stepwise displacements & headings
      ## survival is guaranteed
      X[i, 13] <- 1 
      
      if(pb){
        setTxtProgressBar(tpb, i)
        if(i==N | X[i, 13] == 0) close(tpb)
      }
      
      if(X[i, 13] == 0) break # stop if dead
    }
    ## truncate matrix if death occurs before i = N
    if(length(na.omit(X[, 13])) < N) {
      end <- which(X[, 13] < 1)
      X <- X[1:end, ]
      N <- nrow(X)
    }
    X <- data.frame(X)
    names(X) <- c("x", "y", "dx", "dy", "u", "v", "d", "phi", "rel_mag", "mu_c", "fL", "b", "s")
    
    sim <- X %>% as_tibble() %>%
      mutate(id = rep(id, N)) %>%
      mutate(date = seq(ISOdatetime(2018,07,09,00,00,00), by = 3600, length.out = N)) %>%
      mutate(d_r = sqrt((dx+u)^2 + (dy+v)^2), phi_r = atan2((dx+u), (dy+v)) %% (2*pi)) %>%
      select(id, date, everything())
    
    param <- mpar  
    out <- list(sim = sim, params = param)  
    class(out) <- "simsmolt"
    
    return(out)
    
  }
