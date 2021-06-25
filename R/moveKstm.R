#' @title random walk movement kernel for St Marys River kelts
#'
#' @description utility function not to be called by user
#'
#' @importFrom CircStats rwrpcauchy
#' @importFrom stats rweibull
#' @importFrom raster extract xyFromCell
#' @export
#'
moveKstm <- function(data, xy = NULL, mpar, i, step, ts, w) {

  phi <- NULL
  ## regular movement
  ## calculate distance to land
  d2l <- extract(data$land, rbind(xy))

  if (d2l > mpar$pars$buffer) {
    switch(mpar$scenario,
           rs = {
             if(i < mpar$pars$N/2) {
               ## state 1: migration toward reconditioning area
               ## if current xy not inside land buffer then employ biased migration
               delta <- c(mpar$pars$coa[1,1] - xy[1], mpar$pars$coa[1,2] - xy[2])

               mu <- atan2(delta[1], delta[2])
               rho <- tanh(mpar$pars$r * sqrt(sum(delta^2)))
               phi <- rwrpcauchy(1, mu, rho)

             } else if (i >= mpar$pars$N/2) {
               ## state 2: migration back to spawning river
               mu <- atan2(mpar$pars$coa[2,1] - xy[1], mpar$pars$coa[2,2] - xy[2])
               phi <- rwrpcauchy(1, mu, mpar$pars$rho)
             }
           },
           as = {
             if(i < round(mpar$pars$N * mpar$pars$pN)) {
               ## state 1: migration toward E NL
               if(xy[1] < mpar$pars$NFline.x) {
                 delta <- c(mpar$pars$coa[1,1] - xy[1], mpar$pars$coa[1,2] - xy[2])

                 mu <- atan2(delta[1], delta[2])
                 rho <- tanh(mpar$pars$r * sqrt(sum(delta^2)))
                 phi <- rwrpcauchy(1, mu, rho)

               } else if(xy[1] >=  mpar$pars$NFline.x) {
                ## state 2a: turn toward W Greenland
                delta <- c(mpar$pars$coa[2,1] - xy[1], mpar$pars$coa[2,2] - xy[2])

                mu <- atan2(delta[1], delta[2])
                rho <- tanh(mpar$pars$r * sqrt(sum(delta^2)))
                phi <- rwrpcauchy(1, mu, rho)

               }
               if(xy[2] > 775) {
                  ## state 2: migration toward W Greenland, unconstrained in x axis
                 delta <- c(mpar$pars$coa[2,1] - xy[1], mpar$pars$coa[2,2] - xy[2])

                 mu <- atan2(delta[1], delta[2])
                 rho <- tanh(mpar$pars$r * sqrt(sum(delta^2)))
                 phi <- rwrpcauchy(1, mu, rho)
               }
             } else if (i >= round(mpar$pars$N * mpar$pars$pN)) {
               ## state 3: migration back to E NL
               if(xy[2] > mpar$pars$NFline.y) {
                 mu <- atan2(mpar$pars$coa[3,1] - xy[1], mpar$pars$coa[3,2] - xy[2])
                 phi <- rwrpcauchy(1, mu, mpar$pars$rho)

               } else if(xy[2] <= mpar$pars$NFline.y |
                         xy[1] < mpar$pars$coa[3,1]) {
                 ## state 4: migration back to spawning river
                 mu <- atan2(mpar$pars$coa[4,1] - xy[1], mpar$pars$coa[4,2] - xy[2])
                 phi <- rwrpcauchy(1, mu, mpar$pars$rho)
               }
             }
           })

    ## Temperature-dependent direction reversal (instantaneous)
#    g.rng <- growth(w, seq(1, 25, l = 100), step)
#    tsm <- seq(1, 25, l = 100)[which(g.rng >= w)] * 0.5 #0.75
    if(ts <= 5) {
      cells <- extract(data$ts[[(yday(mpar$pars$start.dt + i * 3600))]],
                       y = cbind(xy[1], xy[2]),
                       buffer = 5,
                       df = TRUE,
                       cellnumbers = TRUE)
      ## take 1st location at max ts
      cell.tmax <- cells[cells[, 3] == max(cells[, 3], na.rm = TRUE), "cells"]
      xys <- xyFromCell(data$ts[[(yday(mpar$pars$start.dt + i * 3600))]], cell.tmax) %>% as.vector()
      mu <- atan2(xys[1] - xy[1], xys[2] - xy[2])
      phi <- rwrpcauchy(1, mu, mpar$pars$rho)

    }

  } else {
    ## direct kelt to move parallel to land
    if(xy[2] < 1200) {
      phi <- as.numeric(extract(data$land_dir, rbind(xy)) + 0.5 * pi)
      if (d2l <= 2) {
        phi <-
          (phi + 0.5 * pi) ## move in opposite direction of land if within 2km
      }
    } else {
      ## if kelt near Greenland then move parallel to land in either direction
      sg <- sample(c(-1,1), size = 1)
      phi <- as.numeric(extract(data$land_dir, rbind(xy)) + sg * 0.5 * pi)
      if (d2l <= 5) {
        phi <-
          (phi + sg * 0.5 * pi) ## move in opposite direction of land if within 5km
      }
    }
  }

  new.xy <- c(xy[1] + sin(phi) * step, xy[2] + cos(phi) * step)

  if(mpar$shelf & d2l <= 20) {
    pv <- c(extract(data$shelf[[1]], rbind(new.xy))[1],
            extract(data$shelf[[2]], rbind(new.xy))[1])
    new.xy <- new.xy + pv * mpar$pars$beta
  }

  ## check if new xy close/on land
  new.d2l <- extract(data$land, rbind(new.xy))

  ## if new location on land (0) then adjust so it's in water
  if(!is.na(new.d2l) & new.d2l == 0) {
    ## find all nearby cells within 5 km & select the one farthest from land
    cells <- extract(data$land, rbind(new.xy), buffer = 5, cellnumbers = TRUE, df = TRUE)
    cell.max <- cells[cells[, 3] == max(cells[, 3], na.rm = TRUE), 2][1]
    new.xy <- xyFromCell(data$land, cell.max) %>% as.vector()

  } else if(is.na(new.d2l)) {
    new.xy <- c(NA,NA)
  }

  cbind(new.xy[1], new.xy[2], phi %% (2*pi))

}
