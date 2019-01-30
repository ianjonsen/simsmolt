#' @title simulate acoustic transmissions & detection, using sim_move & presim output
#' 
#' @description simulates transmissions & detections along simulated track segments within a defined range of acoustic array(s)
#' 
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#' 
#' @param delay - min & max time intervals (s) between transmissions
#' @param burst - duration of each transmission (s)
#' @param psim - a list of required data from \code{presim}
#' @importFrom sp Polygon Polygons SpatialPolygons CRS
#' @importFrom raster buffer
#' @importFrom prevR point.in.SpatialPolygons
#' @importFrom dplyr %>% bind_rows mutate arrange desc
#' @importFrom glatos detect_transmissions
#' @importFrom stats plogis
#' @export
#' 
sim_detect <-
  function(s, delay = c(20,60), burst = 5.0, psim=NULL){
    
    ## simulate tag transmissions along track but only within +/-10 km of avg receiver location
    ##  otherwise trap() output is far too big to generate along full track
    ##    - convert locs from km to m grid; vel in m/s
    ## mean loc for SoBI line
    tr_sobi <- tr_labsea <- NULL
    mrec <- apply(s$data$sobi / 1000, 2, mean)
    dfn <-
      function(r, xy)
        sqrt((r["x"] - xy[, "x"]) ^ 2 + (r["y"] - xy[, "y"]) ^ 2)
    in.rng <- which(dfn(mrec, s$sim[, c("x", "y")]) <= 10)
    if (length(in.rng) > 0)
      tr_sobi <- sim_transmit(s$sim[in.rng, c("x", "y")] * 1000, delayRng = c(20, 60)) ## delay from Chaput et al 2018
   
    ## create 1 km buffer around Lab Sea array
    rb.buff <- s$data$labsea_poly %>% buffer(., 1)
    
    ## determine which smolt track locations are in Lab Sea array
    in.grid <- point.in.SpatialPolygons(s$sim$x, s$sim$y, rb.buff)
    if (sum(in.grid) > 0)
      ## FIXME: will return approx error if only 1 location in grid
      tr_labsea <- sim_transmit(s$sim[in.grid, c("x", "y")] * 1000, delayRng = delay, burstDur = burst)
    
    ## define logistic detection range (m) function
    ## parameterised from analysis of SoBI sentinel tag detections
    ## in July 2009 & July 2010 (see ~/Dropbox/collab/otn/fred/r/fn/sentinel.r)
    pdrf <- function(dm, b = c(3.017, -0.0139)) {
      plogis(b[1] + b[2] * dm)
    }
    
    ## simulate detections given receiver locations & simulated transmission along track
    ##  FIXME: NEED TO ADAPT GLATOS VERSION SO PDRF CAN BE MODIFIED IN FN CALL
    dt_sobi <- dt_labsea <- NULL
    if (!is.null(tr_sobi)) {
      dt_sobi <-
        detect_transmissions(trnsLoc = tr_sobi,
                                     recLoc = s$data$sobi,
                                     detRngFun = pdrf)
    }
    if (!is.null(tr_labsea)) {
      dt_labsea <-
        detect_transmissions(trnsLoc = tr_labsea,
                                     recLoc = s$data$labsea,
                                     detRngFun = pdrf)
    }
    
    if (exists("dt_sobi") & exists("dt_labsea")) {
      s$trans <- bind_rows(tr_sobi, tr_labsea) %>%
        mutate(array = rep(c("sobi", "labsea"), c(nrow(tr_sobi), nrow(tr_labsea)))) %>%
        arrange(desc(array), etime)
        s$detect <- bind_rows(dt_sobi, dt_labsea) %>%
          mutate(array = rep(c("sobi", "labsea"), c(nrow(dt_sobi), nrow(dt_labsea)))) %>%
          arrange(desc(array), etime, recv_id, trns_id)
    } else if (exists("dt_sobi") & !exists("dt_labsea")) {
      s$trans <- tr_sobi %>%
        mutate(array = rep(c("sobi"), c(nrow(tr_sobi)))) %>%
        arrange(desc(array), etime)
      s$detect <- dt_sobi %>%
        mutate(array = rep(c("sobi"), c(nrow(dt_sobi)))) %>%
        arrange(desc(array), etime, recv_id, trns_id)
    } else if (!exists("dt_sobi") & exists("dt_labsea")) {
      s$trans <- tr_labse %>%
        mutate(array = rep(c("labsea"), c(nrow(tr_labsea)))) %>%
  arrange(desc(array), etime)
      s$detect <- dt_labsea %>%
        mutate(array = rep(c("labsea"), c(nrow(dt_labsea)))) %>%
        arrange(desc(array), etime, recv_id, trns_id)
    }
    s
  }