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
#' @export
#' 
presim <-
  function(bathy = "dat/prob_xy.grd",
           uv = NULL,
           sst = NULL,
           rec = TRUE
           ) {
    require(tidyverse)
    
    if(is.null(bathy)) stop("path to bathymetry layer must be supplied\n")
    
    ## load required raster layers
    bathy.xy <- raster::raster(bathy)
    if(!is.null(uv)) uv <- raster::brick(uv)
    
    ## load receiver location data
    if(rec) {
      ## FIXME: this needs to be generalized - do receiver data munging prior to using this function...
      ## FIXME:  could generalize by accessing OTN server to pull requested receiver data from anywhere...
      ## FIXME:  prep code would prob require consistent receiver location / history format on OTN server
      
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
      sp::coordinates(loc) <- c("x","y")
      sp::proj4string(loc) <- sp::CRS(prj_ll)
      locp <- sp::spTransform(loc, sp::CRS(prj_laea)) %>% data.frame()
      asf <- bind_cols(asf, locp) %>%
        dplyr::filter(x >= 0)
      ## get SoBI receiver locations
      sobi <- asf %>% filter(grepl("Strait", locality)) 
      ## set receiver locations in m
      Srecs <- sobi %>% dplyr::select(x, y) * 1000
    }
    if(!is.null(uv)) list(bathy = bathy.xy, uv = uv, sobi = Srecs, prj = prj_laea)
    else list(bathy = bathy.xy, sobi = Srecs, prj = prj_laea)
  }
