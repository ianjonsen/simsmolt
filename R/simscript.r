require(tidyverse)
require(raster)
source("fn/simulate2.r")

bathy.xy <- raster("dat/prob_xy.grd")
uv <- brick("dat/uv_xy.grd")

## for ggplot
bathy.df <- read_csv("dat/prob_df.csv")
u.df <- read_csv("dat/u_df.csv")
v.df <- read_csv("dat/v_df.csv")

## simple topo/bathymetric cost function
## plot(seq(-20,400,l=40), plogis(-5 + 0.2*seq(-20,400,l=40) - 0.001*seq(-20, 400, l=40)^2))
## could be used to turn landmask into a probabilistic map
## rw chooses N steps (displacement & direction) and select most probable outcome by sampling
## prob map and choosing max value
## keeps movements off land

## get ASF receiver locations
asf <- read_csv("dat/stations.csv") %>% 
  dplyr::select(-notes) %>% 
  dplyr::filter(collectioncode=="ASF", 
                grepl("Acoustic", station_type),
                stationstatus == "active") %>% 
  rename(lon=longitude, lat=latitude) %>%
  dplyr::select(-FID, -collectioncode, -station_type, -stationclass, -the_geom)

## project from longlat to laea
prj_ll <- "+proj=longlat +ellps=WGS84"
prj_laea <- "+proj=laea +datum=WGS84 +lat_0=45.00833 +lon_0=-66.99167 +ellps=WGS84 +units=km"
loc <- data.frame(x=asf$lon, y=asf$lat)
coordinates(loc) <- c("x","y")
proj4string(loc) <- CRS(prj_ll)
locp <- sp::spTransform(loc, CRS(prj_laea)) %>% data.frame()
asf <- bind_cols(asf, locp) %>%
  dplyr::filter(x >= 0)
## get SoBI receiver locations
sobi <- asf %>% filter(grepl("Strait", locality)) 
## set receiver locations in m
Srecs <- sobi %>% dplyr::select(x, y) * 1000

## simulate tracks
#set.seed(1)
tag <- c(150,240) #330
coa <- c(900,950)

out <- simulate2(
  N = 1800,
  tag = tag,
  coa = coa,
  a = 3,
  b = 0.864,
  rho = 0.8,
  land = bathy.xy,
  curr = uv,
  ntries = 1000
)


  source("fn/trap.r") ## simulate tag TRansmissions Along Path - adapted from glatos
  
  ## simulate tag transmissions along track but only within +/-10 km of avg receiver location
  ##  otherwise trap() output is far too big to generate along full track
  ##    - convert locs from km to m grid; vel in m/s
  ## mean loc for SoBI line
  trS <- trL <- NULL
  mSrec <- apply(Srecs / 1000, 2, mean)
  dfn <-
    function(r, xy)
      sqrt((r["x"] - xy[, "x"]) ^ 2 + (r["y"] - xy[, "y"]) ^ 2)
  in.rng <- which(dfn(mSrec, out[, c("x", "y")]) <= 10)
  if (length(in.rng) > 0)
    trS <- trap(out[in.rng, c("x", "y")] * 1000, delayRng = c(20, 60)) ## delay from Chaput et al 2018
  
  ## `place` receivers in Labrador Sea & generate Polygon object for track sampling
  Lrecs <-
    expand.grid(x = seq(780, 830, l = 26), y = seq(850, 856, l = 4)) * 1000
  rec <- sp::Polygon(Lrecs[chull(Lrecs),] / 1000)
  ## create SpatialPolygon with a 2km buffer
  rec.box <-
    sp::SpatialPolygons(
      list(Polygons(list(rec), ID = 1)),
      integer(1),
      proj4string = sp::CRS(
        "+proj=laea +datum=WGS84 +lat_0=57.5 +lon_0=55 +units=km +ellps=WGS84 +towgs84=0,0,0"
      )
    )
  rb.buff <- rec.box %>% buffer(., 1)
  in.grid <- prevR::point.in.SpatialPolygons(out$x, out$y, rb.buff)
  if (sum(in.grid) > 0)
    ## FIXME: will return approx error if only 1 location in grid
    trL <- trap(out[in.grid, c("x", "y")] * 1000, delayRng = c(20, 60))
  
  ## define logistic detection range (m) function
  pdrf <- function(d, b = c(0.5,-1 / 120)) {
    1 / (1 + exp(-(b[1] + b[2] * d)))
  }
  
  ## simulate detections given receiver locations & simulated transmission along track
  ##  PROB NEED TO ADAPT GLATOS VERSION SO PDRF CAN BE MODIFIED IN FN CALL
  dtS <- dtL <- ndtL <- ndtS <- h.in.grid <- NULL
  if (!is.null(trS)) {
    dtS <-
      glatos::detect_transmissions(trnsLoc = trS,
                                   recLoc = Srecs,
                                   detRngFun = pdrf)
    ndtS <- ifelse(nrow(dtS) == 0, 0, nrow(dtS))
  }
  if (!is.null(trL)) {
    ## count number of h smolt is in receiver grid
    h.in.grid <-
      prevR::point.in.SpatialPolygons(out$x, out$y, rec.box) %>% sum()
    dtL <-
      glatos::detect_transmissions(trnsLoc = trL,
                                   recLoc = Lrecs,
                                   detRngFun = pdrf)
    ndtL <- ifelse(nrow(dtL) == 0, 0, nrow(dtL))
  }
  time2sobi <- which(out$y > mSrec[2])[1]
  
  if (out[nrow(out), "y"] > ((mSrec[2]) + 0.5)) {
    cat(sprintf("%3.0f d to reach SoBI\n", time2sobi/24))
    cat(sprintf("%3.0d SoBI detections\n", ndtS))
    cat(sprintf("%3.0d h in Lab Sea grid\n", h.in.grid))
    cat(sprintf("%3.0i Lab Sea detections\n", ndtL))
  } else {
    cat("track did not reach SoBI\n")
  }


m <- ggplot(out) + 
  coord_cartesian(xlim = c(tag[1], coa[1]), ylim = c(tag[2], coa[2] * 1.1)) + 
  geom_raster(data = bathy.df, aes(x, y, fill = layer)) + 
  theme_minimal() + 
#  geom_point(aes(x=coa[1],y=coa[2]), colour='red', size = 2) +
  geom_point(aes(x, y), colour = "yellow", size = 0.025) +
  geom_point(data = sobi, aes(x, y), colour = "red", size = 0.01) +
  geom_point(data = Lrecs, aes(x/1000, y/1000), colour = "red", size = 0.01)
  if(!is.null(trS)) {
    m <- m + geom_point(data = trS, aes(x/1000, y/1000), colour = "blue", size = 0.01) + 
      geom_point(data = dtS, aes(recv_x/1000, recv_y/1000), colour = "green", alpha = 0.45) 
    }
  if(!is.null(trL)) {
    m <- m + geom_point(data = trL, aes(x/1000, y/1000), colour = "blue", size = 0.01) + 
      geom_point(data = dtL, aes(recv_x/1000, recv_y/1000), colour = "green", alpha = 0.45)
  }

