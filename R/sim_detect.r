#' @title simulate acoustic transmissions & detection, using sim_move & presim output
#' 
#' @description simulates transmissions & detections along simulated track segments within a defined range of acoustic array(s)
#' 
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#' 
#' @param s - a simsmolt class list containing output from sim_setup and sim_move
#' @param delay - min & max time intervals (s) between transmissions
#' @param burst - duration of each transmission (s)
#' @param b - logistic regression parameters for detection range function (pdrf) - parameters from glm - estimated SoBI sentinel data (old)
#' @importFrom sp Polygon Polygons SpatialPolygons CRS
#' @importFrom raster buffer
#' @importFrom prevR point.in.SpatialPolygons
#' @importFrom dplyr %>% bind_rows mutate arrange desc
#' @importFrom stats plogis
#' @export
#' 
sim_detect <-
  function(s, data, delay = c(20,60), burst = 5.0){
    
    ## simulate tag transmissions along track but only within +/-10 km of avg receiver location
    ##  otherwise trap() output is far too big to generate along full track
    ##    - convert locs from km to m grid; vel in m/s

    recs <- data$recs 
    trans <- tmp.tr <- dt <- tmp.dt <- NULL
    yrec <- recs$y %>% unique()
    
    in.rng <- lapply(1:length(yrec), function(i) {
      which(abs(yrec[i] - s$sim[, "y"]) <= 1.5)
    })

    ## drop rec lines that smolt did not cross
    in.rng <- in.rng[which(sapply(in.rng, length) > 0 )]
    
#    if (any(sapply(in.rng, length) > 0)) {
      trans <- lapply(1:length(in.rng), function(i){
        path <- s$sim[in.rng[[i]], c("id","x","y")]
        path[, c("x","y")] <- path[, c("x","y")] * 1000
        sim_transmit(path, delayRng = delay, burstDur = burst) %>%
        mutate(line = rep(paste0("l", i), nrow(.)))
      }) %>% 
      do.call(rbind, .) 
      
   
    ## define logistic detection range (m) function
    ## parameterised from analysis of SoBI sentinel tag detections
    ## in July 2009 & July 2010 (see ~/Dropbox/collab/otn/fred/r/fn/sentinel.r)
    pdrf <- function(dm, b = c(3.017, -0.0139)) {
      plogis(b[1] + b[2] * dm)
    }
 
    ## simulate detections given receiver locations & simulated transmission along track
    ## FIXME: need to adapt glatos::detect_transmissions so logistic parameters for pdrf can be
    ## FIXME: passed in via sim_detect arguments
#    if(nrow(trans) > 0) {
      recs <- recs %>%
        mutate(x = x * 1000, y = y * 1000)
      detect <- trans %>% 
        group_by(line) %>%
        do(glatos_detect_transmissions(trnsLoc = ., recLoc = recs[, c("x","y","z")], detRngFun = pdrf))
      
      s$trans <- trans %>%
        select(id, line, et, x, y) %>%
        arrange(line, et)
      
      s$detect <- detect %>%
        arrange(line, etime, recv_id, trns_id)
#    }
 
    return(s)
  }