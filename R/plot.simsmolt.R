##' Map simulated track
##'
##' @title plot
##' @param s a fitted object of class simsmolt
##' @param data data object created by sim_setup
##' @param xlim plot x limits
##' @param ylim plot y limits
##' @param ca plot centres of attraction
##'
##' @importFrom ggplot2 ggplot coord_fixed geom_raster aes theme_minimal geom_point fortify
##' @importFrom ggplot2 scale_color_brewer scale_fill_viridis_c geom_contour geom_path
##' @importFrom raster rasterToPoints
##' @importFrom sp spTransform
##' @method plot simsmolt
##' @export

plot.simsmolt <- function(s, data, xlim = NULL, ylim = NULL, ca = FALSE, raster = TRUE, ...) {
  
  bathy.grd <- data$bathy
  bathy <- rasterToPoints(bathy.grd) %>% data.frame()
  
  data(countriesLow, package = "rworldmap")
  coast <- spTransform(countriesLow, data$prj) %>%
    crop(., bathy.grd) %>%
    fortify(.) %>%
    rename(x = long, y = lat)
  
  recs <- data$recs
  
  if (is.null(xlim))
    xlim <- c(0,1250)
  if (is.null(ylim))
    ylim <- c(0,1750)
  
  if (!is.na(class(s)[2]) & (class(s)[2] == "rowwise_df" | class(s)[2] == "grouped_df")) {
    compl <- sapply(s$rep, function(.) !inherits(., "try-error"))
    cat(sprintf("dropping %i failed runs\n\n", sum(!compl)))
    s <- s[compl, ]
    
    if (ca) {
      coa <-
        data.frame(
          x = lapply(s$rep, function(.)
            .$params$coa[,1]) %>% do.call(c, .),
          y = lapply(s$rep, function(.)
            .$params$coa[,2]) %>% do.call(c, .)
        )
    }
    
    detect <-
      lapply(s$rep, function(.)
        .$detect) %>% do.call(rbind, .)
   
    sim <-
      lapply(1:nrow(s), function(i)
        data.frame(id = i, s$rep[[i]]$sim)) %>% do.call(rbind, .)
    sim.last <- lapply(1:nrow(s), function(i)
      data.frame(id = i, s$rep[[i]]$sim[nrow(s$rep[[i]]$sim), ])) %>%do.call(rbind, .)
    
    } else if(is.na(class(s)[2])){
  
  sim <- s$sim
  sim.last <- s$sim[nrow(s$sim), ] 
  detect <- s$detect
  
  if(ca) coa <- data.frame(x = s$params$coa[,1], y = s$params$coa[,2])
    }
  
  ## generate plot
  m <- ggplot() +
    coord_fixed(
      ratio = 1.84,
      xlim = xlim,
      ylim = ylim,
      expand = TRUE
    )
  
  if(raster) {
    m <- m + geom_raster(data = bathy, aes(x, y, fill = z)) +
      scale_fill_viridis_c(direction = 1) +
      theme_dark()
  }

  m <- m + geom_polygon(data = coast, aes_string(x="x", y="y", group="group"), fill = "black")
  
  if (ca) {
    m <- m + geom_point(data = coa,
                        aes(x, y),
                        colour = "green",
                        size = 2)
  }
  m <- m + geom_contour(
    data = bathy,
    aes(x, y, z = z),
    breaks = -900,
    col = "white",
    lwd = 0.25
  )
  
  m <-
    m + geom_point(
      data = recs,
      aes(x, y),
      colour = "blue",
      size = 0.4
    ) +
    geom_path(data = sim,
              aes(x, y, group=id),
              colour = "salmon",
              size = 0.1, alpha = 0.25) +
    
    geom_point(data = subset(sim.last, s == 1), 
               aes(x, y),
               colour = "dodgerblue",
               size = 0.01) +
    
    geom_point(data = subset(sim.last, s == 0), 
               aes(x, y),
               colour = "black",
               size = 0.01)
  
  if (!is.null(detect) && nrow(detect) > 0) {
    m <-
      m + geom_point(data = detect, 
                     aes(recv_x, recv_y), 
                     colour = "red",
                     size = 0.5)
    
  }
  
  return(m) 
  
  
}
