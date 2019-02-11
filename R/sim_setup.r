#' @title Pre-simulation setup
#'
#' @description load required rasters, receiver locations
#'
#'
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#'
#' @param bathy - required bathymetry file
#' @param land - required distance to land to define simulation domain
#' @param land.dir - required direction to nearest land
#' @param b700.dist - required distance to 700 m isobath
#' @param b700.dir -  required direction to 700 m isobath
#' @param tag - start location(s) of simulated animals
#' @param coa - optional Centre-Of-Attraction location(s) to provide movement bias(es)
#' @param uv - optional current layers, supplied as u and v components, velocity must be in m/s
#' @param sst - optional sea surface temperature layer(s)
#' @param rec - optional acoustic receiver locations
#' @importFrom raster raster brick projectRaster
#' @importFrom sp coordinates<- proj4string<- CRS spTransform
#' @importFrom dplyr select filter rename bind_cols %>%
#' @importFrom readr read_csv
#' @export
#'
sim_setup <-
  function(bathy = file.path("..", "simdata", "bathy.nc"),
           land = file.path("..", "simdata", "d2land_xy.grd"),
           land.dir = file.path("..", "simdata", "land_dir2.grd"),
           b700.dist = file.path("..", "simdata", "b700_dist.grd"),
           b700.dir = file.path("..", "simdata", "b700_dir.grd"),
           uv = NULL,
           sst = NULL,
           rec = TRUE) {
    
    if (is.null(land) |
        is.null(land.dir) | is.null(b700.dist) | is.null(b700.dir))
      stop("path to required rasters must be supplied, see '?sim_setup'\n")
    
    ## FIXME: this needs to be generalized - provide spatial extent for query to download ETOPO2 data?
    ## FIXME:   or rely on user supplying their own bathymetry data
    
    prj_laea <-
      "+proj=laea +datum=WGS84 +lat_0=45.00833 +lon_0=-66.99167 +ellps=WGS84 +units=km"
    ## load required raster layers
    bathy <- raster(bathy) %>%
      projectRaster(., crs=prj_laea)
    
      
    land <- raster(land)
    land[land < 1] <- NA
    
    land.dir <- raster(land.dir)
    b700.dist <- raster(b700.dist)
    b700.dir <- raster(b700.dir)
    
    if (!is.null(uv))
      uv <- brick(uv)
    
    ## load receiver location data
    if (rec) {
      ## FIXME: this needs to be generalized - do receiver data munging prior to using this function...
      ## FIXME:  could generalize by accessing OTN server to pull requested receiver data from anywhere...
      ## FIXME:  prep code would prob require consistent receiver location / history format on OTN server
      
      ## get ASF receiver locations
      # asf <- read_csv(file.path("..", "simdata", "stations.csv")) %>%
      #   select(-notes) %>%
      #   filter(
      #     collectioncode == "ASF",
      #     grepl("Acoustic", station_type),
      #     stationstatus == "active"
      #   ) %>%
      #   rename(lon = longitude, lat = latitude) %>%
      #   select(-FID,
      #          -collectioncode,
      #          -station_type,
      #          -stationclass,
      #          -the_geom)
      # 
      # ## project from longlat to laea
      # prj_ll <- "+proj=longlat +ellps=WGS84"
      # prj_laea <-
      #   "+proj=laea +datum=WGS84 +lat_0=45.00833 +lon_0=-66.99167 +ellps=WGS84 +units=km"
      # loc <- data.frame(x = asf$lon, y = asf$lat)
      # coordinates(loc) <- c("x", "y")
      # proj4string(loc) <- CRS(prj_ll)
      # locp <- spTransform(loc, CRS(prj_laea)) %>% data.frame()
      # asf <- bind_cols(asf, locp) %>%
      #   filter(x >= 0)
      # ## get SoBI receiver locations
      # sobi <- asf %>% filter(grepl("Strait", locality))
      # ## set receiver locations in m
      # Srecs <- sobi %>% select(x, y) * 1000
      
      ## Nain, NL - 56.542222, -61.692778
      ## `place` receivers in Labrador Sea & generate Polygon object for track sampling
 
        ## 4 lines from just N of SoBI to Nain, NL
        ## 10 km spacing
        nain <- cbind(-61.692778,56.542222) %>% 
          rgdal::project(., proj=prj_laea) %>%
          as.data.frame()
        names(nain) <- c("x","y")
        recLines <- list(data.frame(line=rep("Lab_4", 22), x=seq(365,575, by = 10), y=rep(1292, 22)))
        recLines[[2]] <- data.frame(line=rep("Lab_3", 36), x=seq(520, 870, by = 10), y=rep(1292-154, 36))
        recLines[[3]] <- data.frame(line=rep("Lab_2", 27), x=seq(740, 1000, by = 10), y=rep(1292-154*2, 27))
        recLines[[4]] <- data.frame(line=rep("Lab_1", 40), x=seq(765, 1160, by = 10), y=rep(1292-154*3, 40))
        recLines <- do.call(rbind, recLines)
      
#      Lrecs <-
#        expand.grid(x = seq(780, 830, l = 26), y = seq(850, 856, l = 4)) * 1000
      
      ## create SpatialPolygon around LabSea receivers (for summary & detection purposes)
#      rec <- Polygon(Lrecs[chull(Lrecs), ] / 1000)
#      rec.box <-
#        SpatialPolygons(list(Polygons(list(rec), ID = 1)),
#                        integer(1),
#                        proj4string = CRS(prj_laea))
      
      
    }
    if (!is.null(uv))
      list(
        bathy = bathy,
        land = land,
        land.dir = land.dir,
        b700.dist = b700.dist,
        b700.dir = b700.dir,
        uv = uv,
        recs = recLines,
        prj = prj_laea
      )
    else
      list(
        bathy = bathy,
        land = land,
        land.dir = land.dir,
        b700.dist = b700.dist,
        b700.dir = b700.dir,
        recs = recLines,
#        labsea_poly = rec.box,
        prj = prj_laea
      )
  }
