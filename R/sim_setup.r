#' @title Pre-simulation setup
#' 
#' @description load required rasters, receiver locations
#' 
#' 
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#' 
#' @param bathy - required bathymetry layer to define simulation domain
#' @param tag - start location(s) of simulated animals
#' @param coa - optional Centre-Of-Attraction location(s) to provide movement bias(es)
#' @param uv - optional current layers, supplied as u and v components, velocity must be in m/s
#' @param sst - optional sea surface temperature layer(s)
#' @param rec - optional acoustic receiver locations
#' @importFrom raster raster brick
#' @importFrom sp coordinates<- proj4string<- CRS spTransform
#' @importFrom dplyr select filter rename bind_cols %>% 
#' @importFrom readr read_csv
#' @export
#' 
sim_setup <-
  function(bathy = file.path("..","simdata","prob_xy.grd"),
           uv = NULL,
           sst = NULL,
           rec = TRUE
           ) {
    require(tidyverse)
    
    if(is.null(bathy)) stop("path to bathymetry layer must be supplied\n")
    
    ## FIXME: this needs to be generalized - provide spatial extent for query to download ETOPO2 data?
    ## FIXME:   or rely on user supplying their own bathymetry data
    
    ## load required raster layers
    bathy.xy <- raster(bathy)
    if(!is.null(uv)) uv <- brick(uv)
    
    ## load receiver location data
    if(rec) {
      ## FIXME: this needs to be generalized - do receiver data munging prior to using this function...
      ## FIXME:  could generalize by accessing OTN server to pull requested receiver data from anywhere...
      ## FIXME:  prep code would prob require consistent receiver location / history format on OTN server
      
      ## get ASF receiver locations
      asf <- read_csv(file.path("..","simdata","stations.csv")) %>% 
        select(-notes) %>% 
        filter(collectioncode=="ASF", 
                      grepl("Acoustic", station_type),
                      stationstatus == "active") %>% 
        rename(lon=longitude, lat=latitude) %>%
        select(-FID, -collectioncode, -station_type, -stationclass, -the_geom)
      
      ## project from longlat to laea
      prj_ll <- "+proj=longlat +ellps=WGS84"
      prj_laea <- "+proj=laea +datum=WGS84 +lat_0=45.00833 +lon_0=-66.99167 +ellps=WGS84 +units=km"
      loc <- data.frame(x=asf$lon, y=asf$lat)
      coordinates(loc) <- c("x","y")
      proj4string(loc) <- CRS(prj_ll)
      locp <- spTransform(loc, CRS(prj_laea)) %>% data.frame()
      asf <- bind_cols(asf, locp) %>%
        filter(x >= 0)
      ## get SoBI receiver locations
      sobi <- asf %>% filter(grepl("Strait", locality)) 
      ## set receiver locations in m
      Srecs <- sobi %>% select(x, y) * 1000
      
      ## `place` receivers in Labrador Sea & generate Polygon object for track sampling
      Lrecs <-
        expand.grid(x = seq(780, 830, l = 26), y = seq(850, 856, l = 4)) * 1000
      
      ## create SpatialPolygon around LabSea receivers (for summary & detection purposes)
      rec <- Polygon(Lrecs[chull(Lrecs),] / 1000)
      rec.box <-
        SpatialPolygons(
          list(Polygons(list(rec), ID = 1)),
          integer(1),
          proj4string = CRS(prj_laea)
        )
      
      
    }
    if(!is.null(uv)) list(bathy = bathy.xy, uv = uv, sobi = Srecs, labsea = Lrecs, prj = prj_laea)
    else list(bathy = bathy.xy, sobi = Srecs, labsea = Lrecs, labsea_poly = rec.box, prj = prj_laea)
  }
