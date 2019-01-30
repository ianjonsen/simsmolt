##' @importFrom dplyr %>%
##' @method summary simsmolt
##' @export
summary.simsmolt <- function(s, ...) {
  if (length(list(...)) > 0) {
    warning("additional arguments ignored")
  }
  
  ## num h smolt is in Lab Sea array
  h.in.grid <-
    point.in.SpatialPolygons(s$sim$x, s$sim$y, s$data$labsea_poly) %>%
    sum()
  tr <- s$trans %>% group_by(array) %>%
    summarise(n())
  dt <- s$detect %>% group_by(array) %>%
    summarise(n())

  mrec <- apply(s$data$sobi, 2, mean)
  time2sobi <- which(s$sim$y > mrec[2])[1] / 24 ## time in days
  
  structure(list(
    time2sobi = time2sobi,
    h.in.grid = h.in.grid,
    trans = tr,
    detections = dt
  ),
  class = "simsmolt")
  
}

##' @method print summary.simsmolt
##' @export
print.summary.simsmolt <- function(x, digits = 3,
                               ...)
{
  
  cat(sprintf("%3.0f d to reach SoBI\n", x$time2sobi/24))
  cat(sprintf("%3.0d h in Lab Sea grid\n", x$h.in.grid))
  cat(sprintf("%3.0i total transmissions\n", x$trans))
  cat(sprintf("%3.0i total detections\n", x$detections))
  
  invisible(x)
}## print.summary.simsmolt
