#' @title simulate acoustic transmissions & detection, using sim_move & presim output
#' 
#' @description simulates transmissions & detections along simulated track segments within a defined range of acoustic array(s)
#' 
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#' 
#' @param delay - min & max time intervals (s) between transmissions
#' @param burst - duration of each transmission (s)
#' @param psim - a list of required data from \code{presim}
#' @export
#' 
sim_detect <-
  function(d, delay = c(20,60), burst = 5.0, psim=NULL){
    
    source("fn/trap.r") ## simulate tag TRansmissions Along Path - adapted from glatos
    
    ## simulate tag transmissions along track but only within +/-10 km of avg receiver location
    ##  otherwise trap() output is far too big to generate along full track
    ##    - convert locs from km to m grid; vel in m/s
    ## mean loc for SoBI line
    tr_sobi <- tr_labsea <- NULL
    mrec <- apply(psim$sobi / 1000, 2, mean)
    dfn <-
      function(r, xy)
        sqrt((r["x"] - xy[, "x"]) ^ 2 + (r["y"] - xy[, "y"]) ^ 2)
    in.rng <- which(dfn(mrec, d[, c("x", "y")]) <= 10)
    if (length(in.rng) > 0)
      tr_sobi <- trap(d[in.rng, c("x", "y")] * 1000, delayRng = c(20, 60)) ## delay from Chaput et al 2018
   
    ## create SpatialPolygon around LabSea receivers with a 2km buffer
    rec <- sp::Polygon(psim$labsea[chull(psim$labsea),] / 1000)
    rec.box <-
      sp::SpatialPolygons(
        list(Polygons(list(rec), ID = 1)),
        integer(1),
        proj4string = sp::CRS(psim$prj)
      )
    rb.buff <- rec.box %>% raster::buffer(., 1)
    in.grid <- prevR::point.in.SpatialPolygons(d$x, d$y, rb.buff)
    if (sum(in.grid) > 0)
      ## FIXME: will return approx error if only 1 location in grid
      tr_labsea <- trap(d[in.grid, c("x", "y")] * 1000, delayRng = delay, burstDur = burst)
    
    ## define logistic detection range (m) function
    ## parameterised from analysis of SoBI sentinel tag detections
    ## in July 2009 & July 2010 (see ~/Dropbox/collab/otn/fred/r/fn/sentinel.r)
    pdrf <- function(dm, b = c(3.017, -0.0139)) {
      plogis(b[1] + b[2] * dm)
    }
    
    ## simulate detections given receiver locations & simulated transmission along track
    ##  FIXME: NEED TO ADAPT GLATOS VERSION SO PDRF CAN BE MODIFIED IN FN CALL
    dt_sobi <- dtL <- ndtL <- ndtS <- h.in.grid <- NULL
    if (!is.null(tr_sobi)) {
      dt_sobi <-
        glatos::detect_transmissions(trnsLoc = tr_sobi,
                                     recLoc = psim$sobi,
                                     detRngFun = pdrf)
      ndt_sobi <- ifelse(nrow(dt_sobi) == 0, 0, nrow(dt_sobi))
    }
    if (!is.null(tr_labsea)) {
      ## count number of h smolt is in receiver grid
      h.in.grid <-
        prevR::point.in.SpatialPolygons(out$x, out$y, rec.box) %>% sum()
      dt_labsea <-
        glatos::detect_transmissions(trnsLoc = tr_labsea,
                                     recLoc = psim$labsea,
                                     detRngFun = pdrf)
      ndt_labsea <- ifelse(nrow(dt_labsea) == 0, 0, nrow(dt_labsea))
    }
    time2sobi <- which(out$y > mrec[2])[1]
    browser()
  }