##' Map simulated track
##'
##' @title plot
##' @param x a fitted object of class simsmolt
##' @param data data object created by sim_setup
##' @param xlim plot x limits
##' @param ylim plot y limits
##' @param raster select a raster (from data) to use as background (bathy, uv, or NULL)
##' @param res downsampling factor to product a plot faster (default = 5)
##' @param rec should receiver locations be displayed (logical)
##' @param track should smolt track(s) be displayed (logical)
##' @param alpha translucence for smolt track(s)
##' @param lwd width of smolt track(s)
##' @param col colour for receiver locations
##' @param size of smolt track end point(s)

##'
##' @importFrom ggplot2 ggplot coord_fixed aes coord_sf geom_sf
##' @importFrom ggplot2 scale_fill_gradientn 
##' @importFrom ggplot2 geom_polygon geom_path theme_minimal geom_point
##' @importFrom ggplot2 aes_string theme element_rect ylab xlab guides guide_legend
##' @importFrom raster extent crop nlayers
##' @importFrom dplyr "%>%"
##' @importFrom sf st_transform st_as_sf
##' @importFrom stars geom_stars st_as_stars st_contour
##' @method plot simsmolt
##' @export

plot.simsmolt <- function(x, data = NULL, xlim = NULL, ylim = NULL, 
                          raster = "ts", res = 5, NLpoly = TRUE, layer = NULL, rec = FALSE, track = TRUE, 
                          alpha = 0.9, lwd = 0.2, size = 0.2, col = "red", pal = "Cividis", maxp = 5e4,
                          crs = "+proj=laea +lat_0=41 +lon_0=-71 +units=km +ellps=WGS84",
                          ...) {
  
  ## process simulated tracks
  if (!is.na(class(x)[2]) && (class(x)[2] == "tbl_df")) {
    
    ## handle multiple tracks
    compl <- sapply(x$rep, function(.) !inherits(., "try-error"))
    cat(sprintf("dropping %i failed runs", sum(!compl)))
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
  
#  data(countriesLow, package = "rworldmap")
#  coast <- spTransform(countriesLow, data$prj) %>%
#    crop(., data$bathy) %>%
#    fortify(.) %>%
#    rename(x = long, y = lat)
  
  if (is.null(xlim))
    xlim <- c(extent(data$bathy)[1], extent(data$bathy)[2])
  if (is.null(ylim))
    ylim <- c(extent(data$bathy)[3], extent(data$bathy)[4])
  

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
      guides(fill = guide_legend(title = "T ºC"))
  
      if(NLpoly) {
        NLp <- readRDS("../simdata/ESRF_regions/NLpoly.RDS")
        m <- m + geom_sf(data = NLp, col = NA, fill = "orange", alpha = 0.25)
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
              colour = "salmon",
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
