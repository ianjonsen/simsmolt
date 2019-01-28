## create landmask & current baselayers for simulations
require(tidyverse)
require(raster)

## baselayer
tmp <- raster("dat/bathy.nc")
prj <- "+proj=laea +datum=WGS84 +lat_0=45.00833 +lon_0=-66.99167 +ellps=WGS84 +units=km"
bathy.xy <- projectRaster(tmp, crs=prj)

## for higher-res option - slows down ggplot but smooths coastline a bit
#ex1 <- extent(bathy.xy)
#tmp2 <- raster(xmn=ex1[1], ymn=ex1[3], xmx=ex1[2], ymx=ex1[4], resolution = c(0.5,0.5), crs = prj)
#hires <- resample(tmp1, tmp2)

## apply logistic function to depth
prob.xy <- calc(bathy.xy, function(x) plogis(-4 - 0.08*x - 0.0001*x^2))
prob.xy[prob.xy > 0.2] <- 1
#prob.xy[prob.xy < 0.001] <- 0

## convert probability surface to df (for ggplot)
prob.df <- rasterToPoints(prob.xy) %>%
  tbl_df()

## write locally
writeRaster(prob.xy, filename = "dat/prob_xy.grd", overwrite = TRUE)
write_delim(prob.df, path = "dat/prob_df.csv", delim = ",")

#bathy.df <- rasterToPoints(bathy.xy) %>% 
#  tbl_df() %>% 
#  mutate(z = ifelse(z < -700, NA, z)) %>%
#  mutate(z = ifelse(z > 0, NA, z)) %>%
#  mutate(z = ifelse(z < 0, 1, z)) %>%
#  mutate(z = ifelse(is.na(z), 0, z))

## convert bathy.df back to a raster
#bathy.xy <- rasterFromXYZ(bathy.df, 
#                          res = c(0.99, 1.82), 
                          #res = c(0.5, 0.5),
#                          crs = "+proj=laea +datum=WGS84 +lat_0=57.5 +lon_0=55 +units=km +ellps=WGS84 +towgs84=0,0,0"
#)
#writeRaster(bathy.xy, filename="dat/bathy_xy.grd")
#write_delim(bathy.df, path="dat/bathy_df.csv", delim=",")

## read AVISO current climatology (mid-June to mid-Nov; 2000 to 2018), 
##  project, resample & crop to bathy.rst extent & resolution
uv <- brick("dat/uv.grd") %>%
  projectRaster(., crs = prj) %>%
  resample(., prob.xy) %>%
  crop(., extent(prob.xy)) 

## replace NA values with mean velocity; to avoid extract() issues in simulation
##  issues caused by replacing true NA's on land are avoided via the land mask
uv[["u"]][which(is.na(values(uv[["u"]])))] <- cellStats(uv[["u"]], mean)
uv[["v"]][which(is.na(values(uv[["v"]])))] <- cellStats(uv[["v"]], mean)

##  mask uv: 1) by bathy.rst NA values; and 2) bt bathy.rst 0 values
uv <- mask(uv, prob.xy, maskvalue=NA) %>%
  mask(., prob.xy, maskvalue=0)
writeRaster(uv, filename="dat/uv_xy.grd", overwrite = TRUE)

## convert uv raster to df
u.df <- rasterToPoints(uv[["u"]]) %>%
  tbl_df()
v.df <- rasterToPoints(uv[["v"]]) %>%
  tbl_df()
write_delim(u.df, path="dat/u_df.csv", delim=",")
write_delim(v.df, path="dat/v_df.csv", delim=",")