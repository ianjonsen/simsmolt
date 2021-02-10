##' Map simulated track
##'
##' @title plot
##' @param s a fitted object of class simsmolt
##' @param data data object created by sim_setup
##' @param raster select a raster (from data) to use as background (bathy, uv, or NULL)
##' @param fast generate plot quickly (poorer resolution; default = TRUE)
##' @param rec should receiver locations be displayed (logical)
##' @param track should smolt track(s) be displayed (logical)
##' @param alpha translucence for smolt track(s)
##' @param lwd width of smolt track(s)
##' @param col colour for receiver locations
##' @param size of smolt track end point(s)
##' @param xlim plot x limits
##' @param ylim plot y limits
##'
##' @importFrom ggplot2 ggplot coord_fixed geom_raster aes theme_minimal coord_sf
##' @importFrom ggplot2 scale_color_brewer scale_fill_gradientn geom_contour 
##' @importFrom ggplot2 geom_polygon geom_path theme_dark fortify geom_point
##' @importFrom ggplot2 aes_string theme_classic theme element_rect ylab xlab
##' @importFrom raster rasterToPoints crop nlayers
##' @importFrom rasterVis gplot
##' @importFrom dplyr "%>%"
##' @importFrom sp spTransform
##' @method plot simsmolt
##' @export

plot.simsmolt <- function(s, data, xlim = NULL, ylim = NULL, 
                          raster = "ts", bathy = FALSE, fast = TRUE, layer = NULL, rec = FALSE, track = TRUE, 
                          alpha = 0.9, lwd = 0.2, size = 0.2, col = "red", pal = "Cividis", maxp = 5e4,
                          crs = "+proj=laea +lat_0=41 +lon_0=-71 +units=km +ellps=WGS84",
                          ...) {

  if(bathy) {
    bathy.c <- rasterToPoints(data$bathy) %>% data.frame()
    names(bathy.c) <- c("x", "y", "z")
  }
    
    if(!fast) {
    if (!is.null(raster)) {
      switch(raster,
             bathy = {
               ras <- bathy.c
             },
             uv = {
               if (nlayers(data$u) > 1) {
                 u <- data$u[[floor(nlayers(data$u) / 2)]]
                 v <- data$v[[floor(nlayers(data$v) / 2)]]
                 ras <- sqrt(u ^ 2 + v ^ 2) * 3.6
               } else {
                 ras <- sqrt(data$u ^ 2 + data$v ^ 2) * 3.6
               }
               
             },
             ts = {
               if (nlayers(data$ts) > 1) {
                 if (is.null(layer))
                   ras <- data$ts[[nlayers(data$ts)]] - 273
                 else
                   ras <- data$ts[[layer]] - 273
                 
               } else {
                 ras <- data$ts - 273
               }
             })
      if (inherits(ras, "RasterLayer")) {
        ras <- rasterToPoints(ras) %>% data.frame()
        names(ras) <- c("x", "y", "z")
      }
    }
  } else {
    if(nlayers(data$ts) > 1) {
      if (is.null(layer))
        ras <- data$ts[[nlayers(data$ts)]] - 273
      else
        ras <- data$ts[[layer]] - 273
    } else {
      ras <- data$ts - 273
    }
    ras.cont <- rasterToPoints(ras) %>% data.frame()
    names(ras.cont) <- c("x", "y", "z")
  }

  
#  data(countriesLow, package = "rworldmap")
#  coast <- spTransform(countriesLow, data$prj) %>%
#    crop(., data$bathy) %>%
#    fortify(.) %>%
#    rename(x = long, y = lat)
  
  if (is.null(xlim))
    xlim <- c(0, 2307)
  if (is.null(ylim))
    ylim <- c(0, 3198)

  if (!is.na(class(s)[2]) && (class(s)[2] == "tbl_df")) {
    compl <- sapply(s$rep, function(.) !inherits(., "try-error"))
    cat(sprintf("dropping %i failed runs\n\n", sum(!compl)))
    s <- s[compl, ]
    
    detect <-
      lapply(s$rep, function(.)
        .$detect) %>% do.call(rbind, .)
   
    sim <- lapply(s$rep, function(.) .$sim) %>% do.call(rbind, .)
    sim.last <- lapply(s$rep, function(.) .$sim[nrow(.$sim), ]) %>% do.call(rbind, .)

    } else if(is.na(class(s)[2])){
  
  sim <- s$sim
  sim.last <- s$sim[nrow(s$sim), ] 
  detect <- s$detect

    }
  
  ## generate plot
    m <- gplot(ras, maxpixels = maxp) +
      geom_raster(aes(fill = value)) +
      geom_contour(
        data = ras.cont,
        aes(x, y, z = z),
        breaks = 6, 
        col = "steelblue1",
        lwd = 0.4
      ) +
      scale_fill_gradientn(colours = hcl.colors(n=100, pal)) +
      theme_dark()
    
    if(bathy) {
      m <- m + 
        geom_contour(
          data = bathy.c,
          aes(x, y, z = z),
          breaks = -900,
          col = "white",
          alpha = 0.5,
          lwd = 0.4
        )
    }
  
  #m <- m + geom_polygon(data = coast, aes_string(x="x", y="y", group="group"), fill = "black")
  
  if(rec) {
  m <-
    m + geom_point(
      data = data$recLocs_asf,
      aes(x, y),
      colour = col,
      size = 0.4
    ) 
  }
  
  if(track) {
   m <- m + geom_path(data = sim,
              aes(x, y, group=id),
              colour = "pink",
              alpha = alpha, size = lwd) +
    
    geom_point(data = sim.last, 
               aes(x, y),
               colour = "dodgerblue",
               size = size)
   }
  
  if (rec && !is.null(detect) && nrow(detect) > 0) {
    m <-
      m + geom_point(data = detect, 
                     aes(recv_x, recv_y), 
                     colour = "red",
                     size = 0.6)
    
  }
  if(!is.null(crs)) {
    m <- m + coord_sf(crs = crs)
  }
  m <- m + ylab("") + xlab("")
  
  return(m) 
  
  
}
