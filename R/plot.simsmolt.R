##' Map simulated track
##'
##' @title plot
##' @param s a fitted object of class simsmolt
##' 
##' @importFrom ggplot2 ggplot coord_cartesian geom_raster aes theme_minimal geom_point
##' @importFrom raster rasterToPoints
##' @method plot simsmolt
##' @export

plot.simsmolt <- function(s, ...) {

bathy <- rasterToPoints(s$data$bathy) %>% data.frame()
names(bathy)[3] <- "z"
  
m <- ggplot() + 
  coord_cartesian(xlim = c(s$data$tag[1], s$data$coa[1]), ylim = c(s$data$tag[2], s$data$coa[2])) + 
  geom_raster(data = bathy, aes(x, y, fill = z)) + 
  theme_minimal() + 
  geom_point(data = s$sim, aes(x, y), colour = "yellow", size = 0.025) +
  geom_point(data = s$data$sobi, aes(x/1000, y/1000), colour = "red", size = 0.01) +
  geom_point(data = s$data$labsea, aes(x/1000, y/1000), colour = "red", size = 0.01)
   if(!is.null(s$detect)) {
     m <- m + geom_point(data = s$detect, aes(recv_x/1000, recv_y/1000, colour = array), alpha = 0.45) 
   }

print(m)

}
