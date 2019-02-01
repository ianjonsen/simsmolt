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
           coa = c(1000, 1000),
           y.chgpts = c(835, 1100),
           data = NULL,
           mpar = list()
           ) {
    ## default move parameters
    mpar.full <- list(
      a = 2,
      b = 0.864,
      rho = 0.8,
      ntries = 100,
      surv = 0.9936,
      taxis = "no"
    )
    pnms <- names(mpar.full)
    
    if (is.null(data))
      stop("Can't find output from sim_setup()\n")
    if (class(data$bathy)[1] != "RasterLayer") stop("bathymetry must be a RasterLayer")
    if (!is.null(data$uv)) {
      if (class(data$uv)[1] != "RasterBrick" || (!names(data$uv) %in% c("u","v"))) 
        stop("current data must be supplied as a RasterBrick with u and v layers")
    }
    if (is.null(y.chgpts)) y.chgpts <- c(5000, 5000)
    
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
      mpar <- mpar.full
    } 
    
    ## define location/survivourship, s matrix & start position
    X <- matrix(NA, N, 3)
    s <- matrix(NA, N, 2)
    X[1,] <- cbind(tag[1], tag[2], 1)
    delta <- theta <- ps <- c()
    thetaL <- runif(1, -20, 10)
    rhoL <- runif(1, 0.6, 0.8)
  
    ## recursion
    for (i in 2:N) {
      if(X[i - 1, 2] < y.chgpts[1]) {
        delta[i] <- sqrt((coa[1] - X[i - 1, 1]) ^ 2 + (coa[2] - X[i - 1, 2]) ^ 2)
        theta[i] <- atan2(coa[1] - X[i - 1, 1], coa[2] - X[i - 1, 2])
      } else if(X[i - 1, 2] >= y.chgpts[1] & X[i - 1, 2] < y.chgpts[2]) {
        ## if in Lab Sea then head NW
        delta[i] <- NA
        theta[i] <- thetaL / 180 * pi
        rho <- rhoL
      } else if(X[ i - 1, 2] >= y.chgpts[2]) {
        delta[i] <- NA
        theta[i] <- (thetaL - 20) / 180 * pi
        rho <- rhoL
      }
    
      switch(mpar$taxis, no = {
        ## draw n proposal steps
        d <-
          rw(
            n = mpar$ntries,
            mu = theta[i],
            rho = mpar$rho,
            a = mpar$a,
            b = mpar$b
          )
        
        ## add advection due to current, if raster is supplied
        if (is.null(data$uv)) {
          tmp <- cbind(X[i - 1, 1] + d[, 1], X[i - 1, 2] + d[, 2])
        } else {
          uv <-
            extract(data$uv, rbind(X[i - 1, 1:2])) * 3600 / 1000 # to convert from m/s to km/h
          tmp <-
            cbind(X[i - 1, 1] + d[, 1] + uv[1], X[i - 1, 2] + d[, 2] + uv[2])
        }
      },
      rp = {
        if (is.null(data$uv)) stop("a current raster must be supplied to simulate rheotaxis")
        uv <- extract(data$uv, rbind(X[i - 1, 1:2]), method = "bilinear") * 3600 / 1000
        mu <- (atan2(uv[1], uv[2]) + pi) %% (2*pi) ## move dir is opposite current dir
        d <- rw(n = mpar$ntries, mu = mu, rho = 0.95, a = mpar$a, b = mpar$b)
        tmp <- cbind(X[i - 1, 1] + d[, 1] , X[i - 1, 2] + d[, 2] ) # FIXME: remove uv advection as a test
      },
      rn = {
        if (is.null(data$uv)) stop("a current raster must be supplied to simulate rheotaxis")
        uv <- extract(data$uv, rbind(X[i - 1, 1:2]), method = "bilinear") * 3600 / 1000
        mu <- atan2(uv[1], uv[2]) %% (2*pi) ## move dir is with current dir
        d <- rw(n = mpar$ntries, mu = mu, rho = 0.95, a = mpar$a, b = mpar$b)
        tmp <- cbind(X[i - 1, 1] + d[, 1] + uv[1], X[i - 1, 2] + d[, 2] + uv[2])
      })
      
      ## check if proposed update is on land or in deep ocean
      p.step <- extract(data$bathy, tmp, method = "simple")
      idx <- which(p.step == max(p.step))[1]
      ps[i] <- p.step[p.step == max(p.step)][1]
      X[i, 1:2] <- tmp[idx,]
      X[i, 3] <- rbinom(1, 1, mpar$surv^(1/24)) # rescale daily survival to hourly prob. 
      if(X[i, 3] == 0) break
    }

    end <- ifelse(length(which(X[, 3] == 0)) == 0, dim(X)[2], which(X[, 3] == 0))
    X <- X[1:end, ]
    X <- data.frame(X)
    names(X) <- c("x", "y", "s")

    sim <- data.frame(X, ps=ps, delta=delta, theta=theta) %>%
      as_tibble()
    
    
    data$tag <- tag
    data$coa <- coa
    data$chgpts <- y.chgpts
    data$taxis <- mpar$taxis
  
    out <- list(sim = sim, data = data)  
    class(out) <- "simsmolt"
    
    return(out)
    
}
