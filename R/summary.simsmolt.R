##' @importFrom dplyr %>% summarise group_by n
##' @importFrom prevR point.in.SpatialPolygons
##' @method summary simsmolt
##' @export
summary.simsmolt <- function(x, ...) {
  if (length(list(...)) > 0) {
    warning("additional arguments ignored")
  }
  
  ## num h smolt is in Lab Sea array
  h.in.grid <-
    point.in.SpatialPolygons(x$sim$x, x$sim$y, x$data$labsea_poly) %>%
    sum()
  tr <- x$trans %>% group_by(array) %>%
    summarise(tr=n())
  dt <- x$detect %>% group_by(array) %>%
    summarise(dt=n())

  mrec <- apply(x$data$sobi, 2, mean) / 1000
  time2sobi <- which(x$sim$y > mrec[2])[1] / 24 ## time in days
  
  structure(list(
    time2sobi = time2sobi,
    h.in.grid = h.in.grid,
    trans = tr,
    dets = dt
  ),
  class = "summary.simsmolt")
  
}

##' @method print summary.simsmolt
##' @export
print.summary.simsmolt <- function(x, digits = 3,
                               ...)
{
  browser()
  cat(sprintf("%3.0f d to reach SoBI\n", x$time2sobi/24))
  cat(sprintf("%3.0d h in Lab Sea grid\n", x$h.in.grid))
  cat(sprintf("%3.0i total transmissions\n", x$trans$tr))
  cat(sprintf("%3.0i total detections\n", x$dets$dt))
  
  invisible(x)
}## print.summary.simsmolt
