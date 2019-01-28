require(tidyverse)

out$etime <- 1:nrow(out)
out$id <- rep(1,nrow(out))
m <- ggplot(out) + 
  coord_cartesian(xlim = c(tag[1], coa[1]), ylim = c(tag[2], coa[2] * 1.1)) + 
  geom_raster(data = bathy.df, aes(x, y, fill = layer)) + 
  theme_minimal() + 
  geom_point(data = sobi, aes(x, y), colour = "red", size = 0.01) +
  geom_point(data = Lrecs, aes(x/1000, y/1000), colour = "red", size = 0.01)
if(!is.null(trS)) {
  m <- m + 
    geom_point(data = trS, aes(x/1000, y/1000), colour = "blue", size = 0.01) + 
    geom_point(data = dtS, aes(recv_x/1000, recv_y/1000), colour = "green", alpha = 0.45) 
}
if(!is.null(trL)) {
  m <- m + 
    geom_point(data = trL, aes(x/1000, y/1000), colour = "blue", size = 0.01) + 
    geom_point(data = dtL, aes(recv_x/1000, recv_y/1000), colour = "green", alpha = 0.45)
}
  m <- m + 
    geom_point(aes(x, y), colour = "yellow", size = 0.025) +
  gganimate::transition_time(etime) +
  gganimate::shadow_wake(wake_length=0.02, size = 0.05, alpha = 0.1) +
  gganimate::ease_aes("linear")
  
gganimate::animate(m, 900, 10, duration = 10)

  
  
  

