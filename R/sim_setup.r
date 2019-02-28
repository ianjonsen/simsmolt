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
#' @param b900.dist - required distance to 900 m isobath
#' @param b900.dir -  required direction to 900 m isobath
#' @param tag - start location(s) of simulated animals
#' @param coa - optional Centre-Of-Attraction location(s) to provide movement bias(es)
#' @param u - optional zonal current layers, velocity must be in m/s
#' @param v - optional meridional current layers, velocity must be in m/s
#' @param sst - optional sea surface temperature layer(s)
#' @param rec - optional acoustic receiver locations
#' @importFrom raster raster brick projectRaster extract
#' @importFrom sp coordinates<- proj4string<- CRS spTransform SpatialPointsDataFrame
#' @importFrom dplyr select filter rename bind_cols %>% tibble
#' @importFrom readr read_csv
#' @export
#'
sim_setup <-
  function(bathy = file.path("..", "simdata", "bathy_xy.grd"),
           land = file.path("..", "simdata", "d2landc_xy.grd"),
           land.dir = file.path("..", "simdata", "land_dirc.grd"),
           b900.dist = file.path("..", "simdata", "b900_dist.grd"),
           b900.dir = file.path("..", "simdata", "b900_dir.grd"),
           u = file.path("..", "simdata", "ru.grd"),
           v = file.path("..", "simdata", "rv.grd"),
           sst = NULL,
           rspace = 10) {
    
    if (is.null(land) |
        is.null(land.dir) | is.null(b900.dist) | is.null(b900.dir))
      stop("path to required rasters must be supplied, see '?sim_setup'\n")
    
    ## FIXME: this needs to be generalized - provide spatial extent for query to download ETOPO2 data?
    ## FIXME:   or rely on user supplying their own bathymetry data
    
    prj_laea <-
      "+proj=laea +datum=WGS84 +lat_0=51.00833 +lon_0=-64.74167 +ellps=WGS84 +units=km"
    ## load required raster layers
    bathy <- raster(bathy)
    
    land <- raster(land)
    #land[land < 1] <- NA
    
    land.dir <- raster(land.dir)
    b900.dist <- raster(b900.dist)
    b900.dir <- raster(b900.dir)
    
    if (all(!is.null(u), !is.null(v)))
      u <- brick(u)
      v <- brick(v)
    
    ## load receiver location data

      ## FIXME: this needs to be generalized - do receiver data munging prior to using this function...
      ## FIXME:  could generalize by accessing OTN server to pull requested receiver data from anywhere...
      ## FIXME:  prep code would prob require consistent receiver location / history format on OTN server
      
 
        ## 4 lines from just N of SoBI to Nain, NL
        ## 10 km spacing
        nain <- cbind(-61.692778,56.542222) %>% 
          rgdal::project(., proj=prj_laea) %>%
          as.data.frame()
        names(nain) <- c("x","y")
        x1 <- seq(615, 1015, by = rspace)
        x2 <- seq(595, 1000, by = rspace)
        x3 <- seq(440, 870, by = rspace)
        x4 <- seq(220, 575, by = rspace)
        recLines <- list(data.frame(line=rep("l1", length(x1)), x=x1, y=rep(607-154*3, length(x1))))
        recLines[[2]] <- data.frame(line=rep("l2", length(x2)), x=x2, y=rep(607-154*2, length(x2)))
        recLines[[3]] <- data.frame(line=rep("l3", length(x3)), x=x3, y=rep(607-154, length(x3)))
        recLines[[4]] <- data.frame(line=rep("l4", length(x4)), x=x4, y=rep(607, length(x4)))
        
        recLines <- do.call(rbind, recLines) 
        recLines$line <- as.character(recLines$line)
        recLines$z <- extract(bathy, recLines[, c("x","y")])
        recLines <- recLines %>% 
          filter(!is.na(x), !is.na(y), !is.na(z)) %>%
          filter(z < -10, z > -600) # drop receivers at >= -10 m depth & < -600 m depth
        
        ## adjust depth, assuming receivers placed 100 m off seafloor
        recLines <- recLines %>%
          mutate(z = ifelse(z < -120, z + 100, z)) %>%
          mutate(id = rownames(.))
        
        
#      Lrecs <-
#        expand.grid(x = seq(780, 830, l = 26), y = seq(850, 856, l = 4)) * 1000
      
      ## create SpatialPolygon around LabSea receivers (for summary & detection purposes)
#      rec <- Polygon(Lrecs[chull(Lrecs), ] / 1000)
#      rec.box <-
#        SpatialPolygons(list(Polygons(list(rec), ID = 1)),
#                        integer(1),
#                        proj4string = CRS(prj_laea))
      
      

    if (all(!is.null(u), !is.null(v)))
      list(
        bathy = bathy,
        land = land,
        land.dir = land.dir,
        b900.dist = b900.dist,
        b900.dir = b900.dir,
        u = u,
        v = v,
        recs = recLines,
        prj = prj_laea
      )
    else
      list(
        bathy = bathy,
        land = land,
        land.dir = land.dir,
        b900.dist = b900.dist,
        b900.dir = b900.dir,
        recs = recLines,
        prj = prj_laea
      )
  }
