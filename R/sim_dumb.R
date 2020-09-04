#' @title simulate dumb movements with or without current advection, using inputs from sim_setup
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
sim_dumb <-
  function(id=1, 
           N = 2760,
           data = NULL,
           mpar = list(),
           pb = TRUE
  ) {
    ## default move parameters
    
    mpar.full <- list(
      start = "def", # or "sobi" - postsmolts start in SoBI or E of northern tip of NF
      a = 2,
      b = 0.864,
      rho = 0.8,
      ntries = 1,
      buffer = c(10, 10),
      mindist = 8,
      maxdist = 850 
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
                   def = c(sample(seq(675, 725, l=50), 1), 82.5)
    )
    
    ## define location/survivourship, sw matrix & start position
    X <- matrix(NA, N, 3)
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
      
      ## apply envt'l forcing
      ## add advection due to current
      u <- extract(data$u, rbind(X[i -  1, 1:2])) * 3.6
      v <- extract(data$v, rbind(X[i -  1, 1:2])) * 3.6
      u <- ifelse(is.na(u), 0, u)
      v <- ifelse(is.na(v), 0, v)
      tmp <-
        cbind(X[i - 1, 1] + d_s[, 1] + u, 
              X[i - 1, 2] + d_s[, 2] + v)
      
      idx <- sample(idx, 1)
      X[i, 1:2] <- tmp[idx, ]
      
      ## survival is guarranteed
      X[i, 3] <- 1 
      
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
