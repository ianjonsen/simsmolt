##' Map simulated track
##'
##' @title plot
##' @param s a fitted object of class simsmolt
##' @param data data object created by sim_setup
##' @param raster select a raster (from data) to use as background (bathy, u, v, or NULL)
##' @param rec should receiver locations be displayed (logical)
##' @param track should smolt track(s) be displayed (logical)
##' @param alpha translucence for smolt track(s)
##' @param lwd width of smolt track(s)
##' @param col colour for receiver locations
##' @param size of smolt track end point(s)
##' @param xlim plot x limits
##' @param ylim plot y limits
##' @param ca plot centres of attraction
##'
##' @importFrom ggplot2 ggplot coord_fixed geom_raster aes theme_minimal
##' @importFrom ggplot2 scale_color_brewer scale_fill_viridis_c geom_contour 
##' @importFrom ggplot2 geom_polygon geom_path theme_dark fortify geom_point
##' @importFrom ggplot2 aes_string theme_classic theme element_rect ylab xlab
##' @importFrom raster rasterToPoints crop calc
##' @importFrom dplyr "%>%"
##' @importFrom sp spTransform
##' @method plot simsmolt
##' @export

plot.simsmolt <- function(s, data, xlim = NULL, ylim = NULL, ca = FALSE, 
                          raster = NULL, rec = TRUE, track = TRUE, 
                          alpha = 0.25, lwd = 0.2, size = 0.2, col = "blue", option = "D", ...) {
  
  bathy.c <- rasterToPoints(data$bathy) %>% data.frame()
  names(bathy.c) <- c("x","y","z")
  
  if(!is.null(raster)) {
    switch(raster, 
         bathy = {
           ras <- bathy.c
         },
         u = {
           
          # if(m==0) ras <- calc(data$u, mean)
          # else ras <- data$u[[m]]
           ras <- calc(data$u, mean)
           ras <- rasterToPoints(ras) %>% data.frame()
           names(ras) <- c("x","y","z")
         },
         v = {
           #if(m==0) ras <- calc(data$v, mean)
           #else ras <- data$v[[m]]
           ras <- calc(data$v, mean)
           ras <- rasterToPoints(ras) %>% data.frame()
           names(ras) <- c("x","y","z")
         })
}
  
  browser() 
  
  data(countriesLow, package = "rworldmap")
  coast <- spTransform(countriesLow, data$prj) %>%
    crop(., data$bathy) %>%
    fortify(.) %>%
    rename(x = long, y = lat)
  
  if (is.null(xlim))
    xlim <- c(0,1250)
  if (is.null(ylim))
    ylim <- c(0,1750)
  
  if (!is.na(class(s)[2]) && (class(s)[2] == "rowwise_df" || class(s)[2] == "grouped_df")) {
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
      data.frame(id = i, s$rep[[i]]$sim[nrow(s$rep[[i]]$sim), ])) %>% do.call(rbind, .)
    
    } else if(is.na(class(s)[2])){
  
  sim <- s$sim
  sim.last <- s$sim[nrow(s$sim), ] 
  detect <- s$detect
  
  if(ca) coa <- data.frame(x = s$params$coa[,1], y = s$params$coa[,2])
    }

  
  ## generate plot
  m <- ggplot() +
    coord_fixed(
      ratio = diff(ylim) / diff(xlim),
      xlim = xlim,
      ylim = ylim,
      expand = FALSE
    )
  
  if(!is.null(raster)) {
    m <- m + geom_raster(data = ras, aes(x, y, fill = z)) +
      scale_fill_viridis_c(direction = 1, guide = "none", option=option) +
      theme_dark() +
      geom_contour(
        data = bathy.c,
        aes(x, y, z = z),
        breaks = -900,
        col = "white",
        lwd = 0.25
      )
    
  } else {
    m <- m + theme_classic() + 
      theme(panel.background = element_rect(grey(0.3))) + 
      geom_contour(
        data = bathy.c,
        aes(x, y, z = z),
        breaks = c(-500, -700, -900, -1100, -1500, -2000, -2500, -5000, -10000),
        col = grey(0.8),
        lwd = 0.05
      ) +
      geom_contour(
        data = bathy.c,
        aes(x, y, z = z),
        breaks = -900,
        col = "white",
        lwd = 0.1
      )
  }

  m <- m + geom_polygon(data = coast, aes_string(x="x", y="y", group="group"), fill = "black")
  
  if (ca) {
    m <- m + geom_point(data = coa,
                        aes(x, y),
                        colour = "green",
                        size = 2)
  }
  
  
  if(rec) {
  m <-
    m + geom_point(
      data = data$recLocs,
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
    
    geom_point(data = subset(sim.last, s == 1), 
               aes(x, y),
               colour = "dodgerblue",
               size = size) +
    
    geom_point(data = subset(sim.last, s == 0), 
               aes(x, y),
               colour = "black",
               size = size)
   }
  
  if (rec && !is.null(detect) && nrow(detect) > 0) {
    m <-
      m + geom_point(data = detect, 
                     aes(recv_x, recv_y), 
                     colour = "red",
                     size = 0.6)
    
  }
  
  m <- m + ylab("") + xlab("")
  
  return(m) 
  
  
}
