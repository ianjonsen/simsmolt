##' Map simulated track
##'
##' @title plot
##' @param x a fitted object of class simsmolt
##' @param data data object created by sim_setup
##' @param xlim plot x limits
##' @param ylim plot y limits
##' @param raster select a raster (from data) to use as background (bathy, uv, or NULL)
##' @param res downsampling factor to product a plot faster (default = 5, full-resolution = 0)
##' @param rec should receiver locations be displayed (logical)
##' @param track should smolt track(s) be displayed (logical)
##' @param last display position at end of simulation only (logical; if T then track forced to be F)
##' @param alpha translucence for smolt track(s)
##' @param lwd width of smolt track(s)
##' @param col colour for receiver locations
##' @param size of smolt track end point(s)

##'
##' @importFrom ggplot2 ggplot coord_fixed aes coord_sf geom_sf
##' @importFrom ggplot2 scale_fill_gradientn stat_smooth geom_line
##' @importFrom ggplot2 geom_path theme_minimal geom_point
##' @importFrom ggplot2 aes_string theme ylab xlab guides guide_legend
##' @importFrom raster extent crop nlayers
##' @importFrom dplyr "%>%"
##' @importFrom patchwork wrap_plots
##' @importFrom sf st_transform st_as_sf st_crop st_bbox
##' @importFrom stars geom_stars st_as_stars st_contour
##' @method plot simsmolt
##' @export

plot.simsmolt <- function(x, data = NULL, xlim = NULL, ylim = NULL, 
                          raster = "ts", res = 5, esrf = TRUE, layer = NULL, rec = FALSE, track = TRUE, last = FALSE,
                          alpha = 0.5, lwd = 0.2, size = 0.2, col = c("salmon", "red", "dodgerblue"), pal = "Cividis", maxp = 5e4,
                          landfill = grey(0.6),
                          crs = "+proj=laea +lat_0=41 +lon_0=-71 +units=km +datum=WGS84",
                          ...) {
  
  ## process simulated tracks
  if (!is.na(class(x)[2]) && (class(x)[2] == "tbl_df")) {
    
    ## handle multiple tracks
    compl <- sapply(x$rep, function(.) !inherits(., "try-error"))
    if(sum(!compl) > 0) cat(sprintf("dropping %i failed runs", sum(!compl)))
    x <- x[compl, ]
    
    detect <-
      lapply(x$rep, function(.)
        .$detect) %>% do.call(rbind, .)
    
    sim <- lapply(x$rep, function(.) .$sim) %>% do.call(rbind, .)
    sim.last <- lapply(x$rep, function(.) .$sim[nrow(.$sim), ]) %>% do.call(rbind, .)
    
  } else if(is.na(class(x)[2])){
    
    ## handle single track
    sim <- x$sim
    sim.last <- x$sim[nrow(x$sim), ] 
    detect <- x$detect
    
  }
  
  if(!is.null(data)) {
    bathy <- stars::st_as_stars(data$bathy)
    bathy.c <- stars::st_contour(bathy, contour_lines = TRUE, breaks = c(-1000, -999.99)) # to get continental shelf without many holes
  }
    
    if (!is.null(raster)) {
      switch(raster,
             uv = {
               if (nlayers(data$u) > 1) {
                 u <- data$u[[ifelse(!is.null(layer), layer, nlayers(data$u))]]
                 v <- data$v[[ifelse(!is.null(layer), layer, nlayers(data$u))]]
                 ras <- sqrt(u ^ 2 + v ^ 2) * 3.6
               } else {
                 ras <- sqrt(data$u ^ 2 + data$v ^ 2) * 3.6
               }
               
             },
             ts = {
               if (nlayers(data$ts) > 1) {
                 ras <- data$ts[[ifelse(!is.null(layer), layer, nlayers(data$u))]] - 273
               } else {
                 ras <- data$ts - 273
               }
             })
      if (inherits(ras, "RasterLayer")) {
        ras <- stars::st_as_stars(ras)
      }
    }
  
  if (is.null(xlim))
    xlim <- c(extent(data$bathy)[1], extent(data$bathy)[2])
  if (is.null(ylim))
    ylim <- c(extent(data$bathy)[3], extent(data$bathy)[4])
  
  ## get coastline
#  coast <- sf::st_as_sf(rworldmap::countriesLow) %>%
#    st_transform(crs = crs) #%>%
    #st_crop(., c("xmin" = xlim[1], "xmax" = xlim[2], "ymin" = ylim[1], "ymax" = ylim[2]))
  
  ## generate plot
    m <- ggplot() +
      geom_stars(data = ras, downsample = res) + 
      geom_sf(
        data = bathy.c,
        col = "white",
        alpha = 0.25,
        lwd = 0.4
      ) +
      scale_fill_gradientn(colours = hcl.colors(n=100, pal), na.value = "transparent") +
      theme_minimal() + 
      guides(fill = guide_legend(title = ifelse(raster == "ts", "T ºC", "||uv|| km/h")))
  
      if(esrf) {
        m <- m + geom_sf(data = data$esrf, col = NA, fill = "orange", alpha = 0.25)
      }
    
  ## add land  
#    m <- m + geom_sf(data = coast, 
#            fill = landfill,
#            lwd = 0
#            )
  
  if(track & !last) {
   m <- m + geom_path(data = sim,
              aes(x, y, group=id),
              colour = col[1],
              alpha = alpha, 
              size = lwd) +
    geom_point(data = sim.last, 
               aes(x, y),
               colour = col[2],
               size = size,
               alpha = 1)
  } else if(last) {
     m <- m + geom_point(data = sim.last, 
                         aes(x, y),
                         colour = col[2],
                         size = size + 0.5,
                         alpha = 1)
   }
  
  if(rec) {
      m <-
        m + geom_point(
          data = data$recLocs_asf,
          aes(x, y),
          colour = col[3],
          size = 0.4
        ) 
    }
    
    
  if (rec && !is.null(detect) && nrow(detect) > 0) {
    m <-
      m + geom_point(data = detect, 
                     aes(recv_x, recv_y), 
                     colour = col[3],
                     size = 0.6)
    
  }
  if(!is.null(crs)) {
    m <- m + 
      coord_sf(crs = crs) +
      xlim(xlim[1], xlim[2]) +
      ylim(ylim[1], ylim[2])
  }
    
  
  m <- m + theme(
    axis.title = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(1, 0.8)
  )
  
  sim24 <- lapply(x$rep, function(.) {
    .$sim[seq(1, nrow(.$sim), by = 24), ]
    }) %>%
    do.call(rbind, .)
  
  ### Generate time-series plots
  ## Mass gain/loss daily time-series
  wp <- ggplot(sim24, aes(date, w, group = id)) +
    geom_line(lwd = 0.25,
              col = grey(0.6),
              alpha = 0.25) +
    theme_minimal() +
    labs(title = "Mass (g)") +
    theme(axis.title = element_blank())
  
  ## Fork-length daily time-series
  fp <- ggplot(sim24, aes(date, fl, group = id)) +
    geom_line(lwd = 0.25,
              col = grey(0.6),
              alpha = 0.25) +
    theme_minimal() +
    labs(title = "Fork-length (m)") +
    theme(axis.title = element_blank())
  
  ## Temperature daily time-series
  tp <- ggplot(sim24, aes(date, ts, group = id)) +
    geom_line(lwd = 0.25,
              col = grey(0.6),
              alpha = 0.25) +
    theme_minimal() +
    labs(title = "Temperature (ºC)") +
    theme(axis.title = element_blank())
  
  design <- "AAAAAA
             AAAAAA
             AAAAAA
             AAAAAA
             BBCCDD"
  wrap_plots(m, wp, fp, tp, design = design)
  
}
