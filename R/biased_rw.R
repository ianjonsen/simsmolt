#' @title biased random walk function
#' 
#' @description utility function not to be called by user
#' 
#' @importFrom CircStats rwrpcauchy
#' @importFrom stats rweibull
#' @importFrom raster extract xyFromCell
#' @export
#' 
biased_rw <- function(n = 1, data, xy = NULL, coa = NULL, dir = NULL, buffer = NULL, rho, a, b) {
  
  if(is.null(coa) & is.null(dir)) stop("Cannot implement a biased random walk without a centre of attraction or direction")
  if(!is.null(coa) & !is.null(dir)) stop("Only one of a centre of attraction or direction may be specified, not both")
 
  if (a > 0)
    st <- rweibull(n, a, b)
  else if(a == 0 & b > 0) {
    st <- b
  } else {
    st <- 0
  }
  
  d2l <- extract(data$land, rbind(xy))

  if(d2l > buffer | xy[1] < 300) {
    if(!is.null(coa) & is.null(dir)) {
    mu <- atan2(coa[1] - xy[1], coa[2] - xy[2]) 
    } else {
      mu <- dir
    }
    
  } else if (xy[1] >= 300 & !all(xy[1] >= 950, 
                                 xy[1] <= 1065, 
                                 xy[2] >= 1200, 
                                 xy[2] <= 1305) & 
             d2l <= buffer) {
    
    ## direct smolt to move eastward & parallel to shore (avoid land) only after passing through most of GoM
    mu <- (extract(data$land_dir, rbind(xy)) + 0.5 * pi) %% (2*pi)
    if(d2l <= 2) {
      mu <- mu + 0.5 * pi %% (2*pi) ## move in opposite direction of land if within 2km
      rho <- 0.95
      }
    } else if(xy[1] >= 300 & all(xy[1] >= 950, 
                                  xy[1] <= 1065, 
                                  xy[2] >= 1200, 
                                  xy[2] <= 1305) & 
              d2l <= buffer) {
      mu <- (extract(data$land_dir, rbind(xy))  + 0.5 * pi) %% (pi)
      if(d2l <= 2) {
        mu <- mu + 0.5 * pi %% (pi) ## move in opposite direction of land if within 2km
        rho <- 0.95
      }
    }
  
  phi <- rwrpcauchy(n, mu, rho)  
  
  new.xy <- c(xy[1] + sin(phi) * st, xy[2] + cos(phi) * st)
  new.d2l <- extract(data$land, rbind(new.xy))
  
  ## if new location on land (0) then adjust so it's in water
  if(new.d2l == 0) {
    ## find all nearby cells within 3 km & select the one farthest from land
    cells <- extract(data$land, rbind(new.xy), buffer = 3, cellnumbers = TRUE, df = TRUE)
    cell.max <- cells[cells[, 3] == max(cells[, 3], na.rm = TRUE), 2][1]
    new.xy <- xyFromCell(data$land, cell.max)
  }
  
  return(new.xy)
}