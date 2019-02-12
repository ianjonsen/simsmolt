##' @importFrom dplyr %>% summarise group_by n
##' @method summary simsmolt
##' @export
summary.simsmolt <- function(x, ...) {
  if (length(list(...)) > 0) {
    warning("additional arguments ignored")
  }

  if(class(x)[2] == "rowwise_df" || class(x)[2] == "grouped_df") {
    ## summarise multiple replicates
    if(names(x)[2] != "rep") stop("expecting simulation output objects to be named 'rep'")
    all.tr <- lapply(x$rep, function(.) .$trans) %>% do.call(rbind, .)
    all.dt <- lapply(x$rep, function(.) .$detect) %>% do.call(rbind, .)
    mNdt <- sapply(x$rep, function(.) (max(.$sim$y) - min(.$sim$y)) / (nrow(.$sim) / 24))
    
    browser() 
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
  cat(sprintf("%3.0f d to reach SoBI\n", x$time2sobi/24))
  cat(sprintf("%3.0d h in Lab Sea grid\n", x$h.in.grid))
  cat(sprintf("%3.0i transmissions on SoBI line, %3.0i transimissions on LabSea array\n", x$trans$tr[2], x$trans$tr[1]))
  cat(sprintf("%3.0i detections on SoBI line, %3.0i detections on LabSea array\n", x$dets$dt[2], x$dets$dt[1]))
  
  invisible(x)
}## print.summary.simsmolt
