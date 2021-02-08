#' @title temperature-biased random walk function
#' 
#' @description utility function not to be called by user
#' 
#' @importFrom CircStats rwrpcauchy
#' @importFrom stats rweibull
#' @importFrom raster extract xyFromCell
#' @export
#' 
temp_brw <- function(n = 1, i = NULL, mpar = NULL, d1 = NULL, data, xy = NULL, buffer = NULL, ts = NULL, ts.rng = NULL, dir = NULL, a, b) {
  
  if (a > 0)
    st <- rweibull(n, a, b)
  else if(a == 0 & b > 0) {
    st <- b
  } else {
    st <- 0
  }

  d2l <- extract(data$land, rbind(xy))
  ## Temperature-independent if in `SoBi box`
  if(all(xy[1] >= 950,
         xy[1] <= 1065,
         xy[2] >= 1200,
         xy[2] <= 1305
  ) & d2l <= buffer & d2l > 2) {
    mu <- (extract(data$land_dir, rbind(xy))  + 0.5 * pi) %% (pi)
    rho <- mpar$rho
  } else if(d2l <= 2) {
      mu <- mu + 0.5 * pi %% (pi) ## move in opposite direction of land if within 2km
      rho <- 0.95
  } else {
    ## Temperature-dependent if outside `SoBi box`
    ## move direction becomes more concentrated around mu as ts declines below lower limit for growth
    rho.l <-
      exp(ts.rng[1] - 1.25 * ts) / (1 + exp(ts.rng[1] - 1.25 * ts))
    rho.u <-
      exp(-ts.rng[2] + 1.25 * ts) / (1 + exp(-ts.rng[2] + 1.25 * ts))
    
    if (rho.l > rho.u) {
      ## move toward warmer water
      rho <- rho.l
      switch(data$ocean,
             cl = {
               ## choose cell with highest ts, within 1 move step (km)
               cells <-
                 extract(
                   data$ts,
                   rbind(xy),
                   buffer = st,
                   cellnumbers = TRUE,
                   df = TRUE
                 )
               cell.max <-
                 cells[cells[, 3] %in% max(cells[, 3], na.rm = TRUE), 2][1]
               cell.xy <- xyFromCell(data$ts, cell.max)
             },
             doy = {
               cells <- extract(
                 data$ts[[yday(mpar$start.dt + i * 3600) - d1]],
                 rbind(xy),
                 buffer = st,
                 cellnumbers = TRUE,
                 df = TRUE
               )
               cell.max <-
                 cells[cells[, 3] %in% max(cells[, 3], na.rm = TRUE), 2][1]
               cell.xy <-
                 xyFromCell(data$ts[[yday(mpar$start.dt + i * 3600) - d1]], cell.max)
             })
      
    } else if (rho.l < rho.u) {
      ## move toward colder water
      rho <- rho.u
      switch(data$ocean,
             cl = {
               ## choose cell with highest ts, within 1 move step (km)
               cells <-
                 extract(
                   data$ts,
                   rbind(xy),
                   buffer = st,
                   cellnumbers = TRUE,
                   df = TRUE
                 )
               cell.max <-
                 cells[cells[, 3] %in% min(cells[, 3], na.rm = TRUE), 2][1]
               cell.xy <- xyFromCell(data$ts, cell.max)
             },
             doy = {
               cells <- extract(
                 data$ts[[yday(mpar$start.dt + i * 3600) - d1]],
                 rbind(xy),
                 buffer = st,
                 cellnumbers = TRUE,
                 df = TRUE
               )
               cell.max <-
                 cells[cells[, 3] %in% min(cells[, 3], na.rm = TRUE), 2][1]
               cell.xy <-
                 xyFromCell(data$ts[[yday(mpar$start.dt + i * 3600) - d1]], cell.max)
             })
    }
    mu <- atan2(cell.xy[1] - xy[1], cell.xy[2] - xy[2])
  }
  
  phi <- rwrpcauchy(n, mu, rho)
  new.xy <- c(xy[1] + sin(phi) * st, xy[2] + cos(phi) * st)
  new.d2l <- extract(data$land, rbind(new.xy))
  
  ## if new location on land (0) then adjust so it's in water
  if(new.d2l == 0) {
    print("\n adjusting location off land")
    ## find all nearby cells within 3 km & select the one farthest from land
    cells <- extract(data$land, rbind(new.xy), buffer = 3, cellnumbers = TRUE, df = TRUE)
    cell.max <- cells[cells[, 3] == max(cells[, 3], na.rm = TRUE), 2][1]
    new.xy <- xyFromCell(data$land, cell.max)
  }

  return(new.xy)
}