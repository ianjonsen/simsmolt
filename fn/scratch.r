require(tidyverse)

## rerddap test
soda_info <- rerddap::info("erdSoda331oceanmday")
test <- rerddap::griddap(soda_info, latitude=c(45.00833,70), longitude=c(290,320), time=c("1980-01-16","1980-01-16"), depth=c(5,5), fields="all")

ggplot(test$data, aes(lon, lat, fill=u)) + geom_raster(interpolate=FALSE) +theme_bw() + scale_fill_viridis_c() + coord_fixed(1.3, xlim = c(290,320),  ylim = c(45, 70))

prj <- "+proj=laea +datum=WGS84 +lat_0=45.00833 +lon_0=293.0083 +ellps=WGS84 +units=km"
u.rst <- test$data %>% 
  select(lon, lat, u) %>%
  raster::rasterFromXYZ(., crs = "+proj=longlat +datum=WGS84") %>%
  raster::projectRaster(., crs=prj)

u.rst <- raster::rasterFromXYZ(test$data, crs = "+proj=longlat +datum=WGS84")






tmp <- raster::raster("bathy.nc")
prj <- "+proj=laea +datum=WGS84 +lat_0=45.00833 +lon_0=-66.99167 +ellps=WGS84 +units=km"
tmp1 <- raster::projectRaster(tmp, crs=prj)

bathy <- raster::rasterToPoints(tmp1) %>% 
  tbl_df() %>% 
  mutate(z = ifelse(z < -700, NA, z)) %>%
  mutate(z = ifelse(z > 0, NA, z)) %>%
  mutate(z = ifelse(z < 0, 1, z)) %>%
  mutate(z = ifelse(is.na(z), 0, z))

ggplot() +
  coord_quickmap(
    expand = FALSE) +
  geom_raster(
    data = bathy,
    aes(x, y, fill = z)
  ) + 
  theme_minimal() 
#  marmap::scale_fill_etopo()

## convert bathy df back to a raster
bathy.rst <- raster::rasterFromXYZ(bathy, 
                                   res = c(0.99, 1.82), 
                                   crs = "+proj=laea +datum=WGS84 +lat_0=57.5 +lon_0=55 +units=m +ellps=WGS84 +towgs84=0,0,0"
                                   )

foo <- crwr(bathy.rst, initPos=c(700,745), nsteps = 200)

ggplot() + coord_quickmap() + geom_raster(data=bathy, aes(x,y,fill=z))+theme_minimal() + geom_path(data=foo, aes(x,y), colour="yellow")


## downscale resolution to make polygon conversion more manageable
#bathy <- raster::aggregate(bathy.rst, factor=5, na.rm = FALSE)


## convert raster to polygon
#bathy.poly <- raster::rasterToPolygons(tmp, 
#                                       n=4, 
#                                       dissolve = TRUE, 
#                                       na.rm = FALSE
#                                       )

#bathy.poly@data$id = rownames(bathy.poly@data)
#bathy.poly.points = broom::tidy(bathy.poly, region="id")
#bathy.poly.df = plyr::join(bathy.poly.points, bathy.poly@data, by="id")
#bathy2 <- subset(bathy.poly.df, id ==2)
#saveRDS(bathy, file="bathy_poly_df.RDS")

## visualise
#m <- ggplot(bathy.poly.df) + 
#  aes(long, lat, group = group, fill=id) + 
#  geom_polygon()
