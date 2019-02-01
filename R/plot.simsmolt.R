##' Map simulated track
##'
##' @title plot
##' @param s a fitted object of class simsmolt
##' 
##' @importFrom ggplot2 ggplot coord_fixed geom_raster aes theme_minimal geom_point scale_color_brewer scale_fill_gradient2
##' @importFrom raster rasterToPoints
##' @method plot simsmolt
##' @export

plot.simsmolt <- function(s, ...) {

bathy <- rasterToPoints(s$data$bathy) %>% data.frame()
names(bathy)[3] <- "z"
  
m <- ggplot() + 
  coord_fixed(ratio = 1.84, xlim = c(100, 1250), ylim = c(700, 2000)) + 
  geom_raster(data = bathy, aes(x, y, fill = z)) + 
  scale_fill_gradient2(low = "#053061", mid = "#43a2ca", high = "#e0f3db", guide = "none") +
  theme_minimal() + 
  geom_point(data = s$sim, aes(x, y), colour = "firebrick", size = 0.025) +
  geom_point(data = s$data$sobi, aes(x/1000, y/1000), colour = "#2166ac", size = 0.01) +
  geom_point(data = s$data$labsea, aes(x/1000, y/1000), colour = "#2166ac", size = 0.01)
   if(!is.null(s$detect)) {
     m <- m + geom_point(data = s$detect, aes(recv_x/1000, recv_y/1000, colour = array)) +
       scale_color_brewer(type = "qual", palette = 1)
   }

return(m)

}
