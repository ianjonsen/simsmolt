##' Map simulated track
##'
##' @title plot
##' @param s a fitted object of class simsmolt
##' @param xlim plot x limits
##' @param ylim plot y limits
##'
##' @importFrom ggplot2 ggplot coord_fixed geom_raster aes theme_minimal geom_point scale_color_brewer scale_fill_gradient2
##' @importFrom raster rasterToPoints
##' @method plot simsmolt
##' @export

plot.simsmolt <- function(s, xlim = NULL, ylim = NULL, ...) {
  bathy <- s$data$bathy
  bathy[bathy > 0] <- NA
  bathy <- rasterToPoints(bathy) %>% data.frame()
  land <- rasterToPoints(s$data$land) %>% data.frame()
  
  names(land)[3] <- "d"
  
  if (is.null(xlim))
    xlim <- extendrange(s$sim$x, f = 1)
  if (is.null(ylim))
    ylim <- extendrange(s$sim$y, f = 0.25)
  
  coa <- data.frame(x = s$data$coa[1], y = s$data$coa[2])
  
  m <- ggplot() +
    coord_fixed(
      ratio = 1.84,
      xlim = xlim,
      ylim = ylim,
      expand = TRUE
    ) +
    #  coord_quickmap(xlim = xrng, ylim = yrng) +
    geom_raster(data = land, aes(x, y, fill = d)) +
    #  scale_fill_gradient2(low = "#053061", mid = "#43a2ca", high = "#e0f3db", guide = "none") +
    scale_fill_viridis_c(direction = -1) +
    theme_minimal() +
    geom_point(data = coa,
               aes(x, y),
               colour = "red",
               size = 1) +
    geom_contour(
      data = bathy,
      aes(x, y, z = z),
      breaks = seq(-700, -100, by = 100),
      col = "white",
      lwd = 0.25
    )
  if (!is.null(s$trans) && nrow(s$trans) > 0) {
    m <-
      m + geom_point(
        data = s$trans,
        aes(x / 1000, y / 1000),
        col = "lightblue",
        alpha = 0.35,
        size = 0.75
      )
  }
  m <-
    m + geom_point(
      data = s$data$recs,
      aes(x, y),
      colour = "blue",
      size = 0.4
    ) +
    geom_point(data = s$sim,
               aes(x, y),
               colour = "firebrick",
               size = 0.01)
  
  if (!is.null(s$detect) && nrow(s$detect) > 0) {
    m <-
      m + geom_point(data = s$detect, aes(recv_x / 1000, recv_y / 1000, colour = line)) +
      scale_color_brewer(type = "qual", palette = 1)
  }
  
  return(m)
  
}
