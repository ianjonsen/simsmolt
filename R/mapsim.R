##' Produce a simple map of simulated track(s) with or without receivers
##'
##' @title mapsim
##' @param x a fitted object of class simsmolt
##' @param data data object created by sim_setup
##' @param xlim plot x limits
##' @param ylim plot y limits
##' @param res downsampling factor to product a plot faster (default = 5, full-resolution = 0)
##' @param rec should receiver locations be displayed (logical)
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
##' @importFrom patchwork wrap_plots area
##' @importFrom sf st_transform st_as_sf st_crop st_bbox
##' @importFrom stars geom_stars st_as_stars st_contour
##' @importFrom wesanderson wes_palette
##' @export

mapsim <- function(x, data = NULL, xlim = NULL, ylim = NULL, 
                          res = 5, esrf = FALSE, rec = TRUE, 
                          track = TRUE, last = FALSE,
                          alpha = 0.5, lwd = 0.25, size = 0.2, col = "red", pal = "Blues 3", 
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
  
  wespal <- wes_palette("Darjeeling2", type = "discrete")
  
  if (is.null(xlim))
    xlim <- c(extent(data$bathy)[1], extent(data$bathy)[2])
  if (is.null(ylim))
    ylim <- c(extent(data$bathy)[3], extent(data$bathy)[4])
  
  if(!is.null(data)) {
    bathy <- stars::st_as_stars(data$bathy)
#    bathy.c <- stars::st_contour(bathy, contour_lines = TRUE, breaks = c(-1000, -999.99)) # to get continental shelf without many holes
  }
  
  ## generate plot
  m <- ggplot() +
    geom_stars(data = bathy, downsample = res) + 
    # geom_sf(
    #   data = bathy.c,
    #   col = "white",
    #   alpha = 0.75,
    #   lwd = 0.4
    # ) +
    scale_fill_gradientn(colours = hcl.colors(n=100, pal), na.value = grey(0.8)) +
    theme_minimal()
  
  if(esrf) {
    m <- m + geom_sf(data = data$esrf, col = NA, fill = "white", alpha = 0.25)
  }
  
  if(rec) {
    m <-
      m + geom_point(
        data = data$recLocs,
        aes(x, y),
        colour = col,
        size = 1
      ) 
  }
  
  if(track) {
    m <- m + geom_path(data = sim,
                       aes(x, y, group=id),
                       alpha = alpha, 
                       size = lwd,
                     col = wespal[2])
    }
  if(last){
    m <- m + geom_point(data = sim.last, 
                 aes(x, y),
                 colour = wespal[1],
                 size = size,
                 alpha = 1)
  }
  
  if(!is.null(crs)) {
    m <- m + 
      coord_sf(crs = crs) +
      xlim(xlim[1], xlim[2]) +
      ylim(ylim[1], ylim[2])
  }
  
  
  m <- m + theme(
    axis.title = element_blank(),
    legend.position = "none"
  )
  
  suppressWarnings(m)
}
    
  