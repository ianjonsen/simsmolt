##' Map simulated track
##'
##' @title plot
##' @param s a fitted object of class simsmolt
##' @param xlim plot x limits
##' @param ylim plot y limits
##'
##' @importFrom ggplot2 ggplot coord_fixed geom_raster aes theme_minimal geom_point 
##' @importFrom ggplot2 scale_color_brewer scale_fill_viridis_c geom_contour geom_path
##' @importFrom raster rasterToPoints
##' @method plot simsmolt
##' @export

plot.simsmolt <- function(s, data, xlim = NULL, ylim = NULL, ca = FALSE, ...) {
  
  if (!is.na(class(s)[2]) & (class(s)[2] == "rowwise_df" | class(s)[2] == "grouped_df")) {
    compl <- sapply(s$rep, function(.) !inherits(., "try-error"))
    cat(sprintf("dropping %i failed runs\n\n", sum(!compl)))
    s <- s[compl, ]
    
    bathy <- data$bathy
    bathy[bathy > 0] <- NA
    bathy <- rasterToPoints(bathy) %>% data.frame()
    land <- data$land
    land <- rasterToPoints(land) %>% data.frame()
    names(land)[3] <- "d"
    
    if (is.null(xlim))
      xlim <- extendrange(sapply(s$rep, function(.)
        .$sim$x), f = 1)
    if (is.null(ylim))
      ylim <- extendrange(sapply(s$rep, function(.)
        .$sim$y), f = 0.25)
    
    if (ca) {
      coa <-
        data.frame(
          x = sapply(s$rep, function(.)
            .$params$coa[1]),
          y = sapply(s$rep, function(.)
            .$params$coa[2])
        )
      coa$y <- ifelse(coa$y > ylim[2], ylim[2] - 0.5, coa$y)
    }
    
    trans <-
      lapply(s$rep, function(.)
        .$trans) %>% do.call(rbind, .)
    detect <-
      lapply(s$rep, function(.)
        .$detect) %>% do.call(rbind, .)
    recs <- data$recs
    sim <-
      lapply(1:nrow(s), function(i)
        data.frame(id = i, s$rep[[i]]$sim)) %>% do.call(rbind, .)
    
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
      theme_minimal()
    if (ca) {
      m <- m + geom_point(data = coa,
                          aes(x, y),
                          colour = "red",
                          size = 1)
    }
    m <- m + geom_contour(
      data = bathy,
      aes(x, y, z = z),
      breaks = seq(-400, -700, by = -100),
      col = "white",
      lwd = 0.2
    )
    if (!is.null(trans) && nrow(trans) > 0) {
      m <-
        m + geom_point(
          data = trans,
          aes(x / 1000, y / 1000),
          col = "lightblue",
          alpha = 0.35,
          size = 0.75
        )
    }
    m <-
      m + geom_point(
        data = recs,
        aes(x, y),
        colour = "blue",
        size = 0.4
      ) +
      geom_path(data = sim,
                 aes(x, y, group=id),
                 colour = "firebrick",
                 size = 0.1, alpha = 0.25)
    
    if (!is.null(detect) && nrow(detect) > 0) {
      m <-
        m + geom_point(data = detect, aes(recv_x / 1000, recv_y / 1000, colour = line)) +
        scale_color_brewer(type = "qual", palette = 1)
    }
    
    return(m) 
    
    } else if(is.na(class(s)[2])){
      
  bathy <- data$bathy
  bathy[bathy > 0] <- NA
  bathy <- rasterToPoints(bathy) %>% data.frame()
  land <- rasterToPoints(data$land) %>% data.frame()
  
  names(land)[3] <- "d"
  
  if (is.null(xlim))
    xlim <- extendrange(s$sim$x, f = 1)
  if (is.null(ylim))
    ylim <- extendrange(s$sim$y, f = 0.25)
  if(ca) {
  coa <- data.frame(x = s$params$coa[1], y = s$params$coa[2])
  }
  
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
    theme_minimal() 
  
  if(ca) {
    m <- m + geom_point(data = coa,
               aes(x, y),
               colour = "red",
               size = 1) 
  }
  
    m <- m + geom_contour(
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
      data = data$recs,
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
}
