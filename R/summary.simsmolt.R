##' @importFrom dplyr %>%
##' @method summary simsmolt
##' @export
summary.simsmolt <- function(s, ...) {
  if (length(list(...)) > 0) {
    warning("additional arguments ignored")
  }
  
  ## count number of h smolt is in Lab Sea array
  h.in.grid <-
    point.in.SpatialPolygons(s$sim$x, s$sim$y, s$data$labsea_poly) %>%
    sum()
  
  dt <- s$detect %>% group_by(array) %>%
    summarise(n())

  mrec <- apply(s$data$sobi, 2, mean)
  time2sobi <- which(s$sim$y > mrec[2])[1] / 24 ## time in days
  
  
}