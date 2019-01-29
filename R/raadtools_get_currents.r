## script to grab AVISO geostrophic currents
##  using raadtools via AMP VM in Hobart

## define annual date range - based on (rough) avg deployment 
##    dates & detection dates on SoBI line
##    mid-Jun to mid-Nov
dts <- lapply(0:18, function(i){
  seq(ISOdate(2000,6,15) + 86400 * (365*i), 
      by="week", 
      length.out=23)
  }) %>% 
  do.call(c, .) %>%
  format(., "%Y-%m-%d")

## get current data - u & v separately
utmp <- readcurr(date=dts, "weekly", 
                 uonly=TRUE, 
                 xylim = extent(c(290-360,320-360,45,70)))
vtmp <- readcurr(date=dts, "weekly", 
                 vonly=TRUE, 
                 xylim = extent(c(290-360,320-360,45,70)))

## create climatologies
u <- mean(utmp)
v <- mean(vtmp)

## put them together & export
uv <- brick(u,v)
names(uv) <- c("u","v")
writeRaster(uv, filename="uv")

