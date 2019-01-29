##' Visualise fixed and random relationships
##'
##' @title plot
##' @param m a fitted object of class simsmolt
##' @param ani logical (default FALSE) to indicate whether plot should be rendered as an animation
##' @param nf number of frames to render
##' @param fps frame rate of animation in frames / s
##' 
##' @importFrom ggplot2 ggplot coord_cartesian geom_raster aes theme_minimal geom_point
##' @importFrom gganimate transition_time shadow_wake ease_aes animate
##' @importFrom tibble as_tibble
##' @importFrom raster rasterToPoints
##' @method plot simsmolt
##' @export

plot.simsmolt <- function(sim, ani = FALSE, nf = 100, fps = 10, ...) {

bathy <- rasterToPoints(sim$data$bathy) %>% as_tibble()
  
m <- ggplot(out) + 
  coord_cartesian(xlim = c(tag[1], coa[1]), ylim = c(tag[2], coa[2] * 1.1)) + 
  geom_raster(data = bathy, aes(x, y, fill = layer)) + 
  theme_minimal() + 
  geom_point(aes(x, y), colour = "yellow", size = 0.025) +
  geom_point(data = out$data$sobi, aes(x, y), colour = "red", size = 0.01) +
  geom_point(data = out$data$labsea, aes(x/1000, y/1000), colour = "red", size = 0.01) 

if(ani) {
  m <- m +
    transition_time(etime) +
    shadow_wake(wake_length=0.02, size = 0.05, alpha = 0.1) +
    ease_aes("linear")
  
  animate(m, nf, fps)
}

#if(!is.null(trS)) {
#  m <- m + geom_point(data = trS, aes(x/1000, y/1000), colour = "blue", size = 0.01) + 
#    geom_point(data = dtS, aes(recv_x/1000, recv_y/1000), colour = "green", alpha = 0.45) 
#}
#if(!is.null(trL)) {
#  m <- m + geom_point(data = trL, aes(x/1000, y/1000), colour = "blue", size = 0.01) + 
#    geom_point(data = dtL, aes(recv_x/1000, recv_y/1000), colour = "green", alpha = 0.45)
#}

}
