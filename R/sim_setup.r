#' @title Pre-simulation setup
#'
#' @description load required rasters, receiver locations
#'
#'
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#'
#' @param config - path to config.R script containing file.paths for required & optional environmental layers (see Details)
#' @param rec - optional acoustic receiver locations (default is no receivers)
#' @param rspace - nominal spacing (km) between receivers
#' @param rnum - the number of receivers to be used on grid arrays (to approx. match to num used on lines)
#' @param ocean - temporal form of environmental data layers: "cl" = climatological; "doy" = daily means; FALSE = no env't data
#' @param doy.rng - rnage of days of year to be loaded from disk as a rasterStack (default = 1:365)
#' @param uv - are current velocity data to be loaded from disk
#' @importFrom raster raster stack brick projectRaster extract
#' @importFrom sp coordinates<- proj4string<- CRS spTransform SpatialPointsDataFrame spsample
#' @importFrom sf st_as_sf st_sample st_coordinates st_distance
#' @importFrom dplyr select filter rename bind_cols %>% tibble distinct
#' @export
#'
sim_setup <-
  function(config = file.path("/", "Users", "jonsen", "OneDrive - Macquarie University", "collab", "otn", "simdata", "config.R"), 
           rec = "none", rspace = NULL, rnum = NULL,
           ocean = "doy", doy.rng = NULL, uv = FALSE, crs = NULL, esrf = TRUE) {
    
    ## FIXME: this needs to be generalized - provide spatial extent for query to download ETOPO2 data?
    ## FIXME:   or rely on user supplying their own bathymetry data
       
    ## load receiver location data

      ## FIXME: this needs to be generalized - do receiver data munging prior to using this function...
      ## FIXME:  could generalize by accessing OTN server to pull requested receiver data from anywhere...
      ## FIXME:  prep code would prob require consistent receiver location / history format on OTN server

    source(config)  
    if(is.null(crs)) prj <- "+proj=laea +lat_0=41 +lon_0=-71 +units=km +datum=WGS84"
    
      if(rec == "lines") {
        ## 4 lines from just N of SoBI to Nain, NL
        ## 10 km spacing
       # nain <- cbind(-61.692778,56.542222) %>% 
          #rgdal::project(., proj=prj_laea) %>%
          #as.data.frame()
        #names(nain) <- c("x","y")
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
        recLocs <- recLines %>%
          mutate(z = ifelse(z < -120, z + 100, z)) %>%
          mutate(id = rownames(.))

        recPoly_sf <- NULL
        
      } else if (rec == "rnd") {
        ## 1 vlarge grid w random placement
        poly <- Polygon(data.frame(x = c(615,595,440,320,220, 450,600,790,900,1015, 615), 
                           y = c(145,299,453,520,607, 607,520,453,299,145, 145)))
        recPoly_sf <- SpatialPolygons(list(Polygons(list(poly), ID = 1)), 
                                   integer(1), proj4string = CRS(prj)) %>%
          st_as_sf()
        
        ## use same random sample each time (for now)
        set.seed(10)
        recLocs <- st_sample(recPoly_sf, size = 350) %>%
          st_coordinates() %>%
          as_tibble() %>%
          rename(x = X, y = Y)
        ## equivalent # of recs to rpsace = 5
        recLocs <- recLocs[1:350, ]
        recLocs <- recLocs %>%
          mutate(z = extract(bathy, recLocs[, c("x","y")])) %>%
          filter(z > -600, z < -10)
        if (rspace == 10) {
          n <- nrow(recLocs)
          recLocs <- recLocs[sample(1:n, round(n * 0.5)), ]
        }
        ## adjust depth, assuming receivers placed 100 m off seafloor
        recLocs <- recLocs %>%
          mutate(z = ifelse(z < -120, z + 100, z)) %>%
          mutate(id = rownames(.))
        
        
        dist <- recLocs %>%
          st_as_sf(., coords = c("x","y"), crs = prj) %>%
          st_distance()
        diag(dist) <- NA
        min.dist <- apply(dist, 2, min, na.rm = TRUE)
        
      } else if (rec == "grid") {
        ## 1 vlarge grid w random placement

        poly <- Polygon(data.frame(x = c(615,595,440,320,220, 450,600,790,900,1015, 615), 
                                   y = c(145,299,453,520,607, 607,520,453,299,145, 145)))
        recPoly <- SpatialPolygons(list(Polygons(list(poly), ID = 1)), 
                                      integer(1), proj4string = CRS(prj)) 
        
        
        ## use same random sample each time b/c we want grid to always be in same place
        set.seed(pi)  #rnum = 165, 135
        grid <- spsample(recPoly, n=rnum, type="regular", proj4string = prj) %>% st_as_sf()
        
        recPoly_sf <- st_as_sf(recPoly)
        recLocs <- grid %>% st_coordinates() %>%
          as_tibble() %>%
          rename(x = X, y = Y)
        ## aadd & filter on bathymetry
        recLocs <- recLocs %>%
          mutate(z = extract(raster(bathy), recLocs[, c("x","y")])) %>%
          filter(z > -600, z < -10)
        ## adjust depth, assuming receivers placed 100 m off seafloor
        recLocs <- recLocs %>%
          mutate(z = ifelse(z < -120, z + 100, z)) %>%
          mutate(id = rownames(.))
      
      } else if (rec == "asf") {
        # stn <- read.csv(file.path(recs, "stations.csv")) %>%
        #   filter(stationstatus == "active" & stationclass == "deployed") %>%
        #   rename(lat = latitude, lon = longitude) %>%
        #   filter(lat >= 41, lat <= 68, lon >= -71, lon <= -43) %>%
        #   select(-notes, -the_geom)
        # stn <- stn[grep("Acoustic", stn$station_type), ]
        # stn <- stn %>%
        #   sf::st_as_sf(coords = c("lon","lat"), crs = 4326) %>%
        #   sf::st_transform(., crs = prj)
        # recLocs <- sf::st_coordinates(stn) %>% as.data.frame()
        # names(recLocs) <- c("x","y")
        # sf::st_geometry(stn) <- NULL
        # stn <- cbind(stn, recLocs)
        
        ## grab ASF - PHS/SOBI receiver details
        asf <- read.csv(file.path(recs, "ASF_2017-2020_Tx_SOBIandPHS.csv")) 
        names(asf) <- tolower(names(asf))
        asf_stn <- asf %>%
          select(date, receiver, year, otn_array, station_name, locality, region, lat, long) %>%
          rename(lon = long) %>%
          distinct(receiver, year, .keep_all = TRUE) %>%
          sf::st_as_sf(coords = c("lon","lat"), crs = 4326) %>%
          sf::st_transform(., crs = prj)
        recLocs_asf <- sf::st_coordinates(asf_stn) %>% as.data.frame()
        names(recLocs_asf) <- c("x","y")
        sf::st_geometry(asf_stn) <- NULL
        asf_stn <- cbind(asf_stn, recLocs_asf)
        
      } else if(rec %in% c("esrf_g", "esrf_l")) {

        switch(rec,
               esrf_g = {
                 recs <- read.csv(file.path(recs, "esrf_grid.csv"))
               },
               esrf_l = {
                 recs <- read.csv(file.path(recs, "esrf_lines.csv"))
               })

        recs <- recs %>%
          rename(lon=LONG, lat=LAT) %>%
          sf::st_as_sf(., coords = c("lon","lat"), crs = 4326) %>% 
          sf::st_transform(., crs = prj) %>% 
          sf::st_coordinates() %>% 
          as_tibble() %>% 
          rename(x=X, y=Y) 
        recs <- recs %>%
          mutate(z = extract(raster(bathy), recs[, c("x","y")])) %>%
          mutate(id = rownames(.))
        if(rec == "esrf_grid") {
          poly <- sp::Polygon(recs[c(chull(recs[, c("x","y")]),54), c("x","y")])
          recPoly <- SpatialPolygons(list(Polygons(list(poly), ID = 1)), 
                                   integer(1), proj4string = CRS(prj))
          recPoly_sf <- st_as_sf(recPoly)
        }
    }
    
    out <- list(
      bathy = raster(bathy),
      land = raster(d2land),
      land_dir = raster(land_dir)
    )

    if (!is.null(ocean)) {
      switch(ocean, 
             cl = {
               if(uv) {
                out[["u"]] <- raster(file.path(riops, "riops_365_u.grd"))
                out[["v"]] <- raster(file.path(riops, "riops_365_v.grd"))
               }
               out[["ts"]] <- raster(file.path(riops, "riops_365_t.grd"))
             },
          doy = {
            if(is.null(doy.rng)) doy.rng <- 1:365
            if(uv) {
              out[["u"]] <- stack(file.path(riops, "riops_doy_u.grd"))[[doy.rng]] 
              out[["v"]] <- stack(file.path(riops, "riops_doy_v.grd"))[[doy.rng]] 
            }
            out[["ts"]] <- stack(file.path(riops, "riops_doy_t.grd"))[[doy.rng]]
            
          })
    }
    
    if(!rec %in% c("none", "asf", "esrf_g", "esrf_l")) {
      out[["recLocs"]] <- recLocs
      out[["recPoly"]] <- recPoly_sf
      out[["rec"]] <- rec
    } else if (rec == "asf") {
      #out[["recLocs"]] <- stn
      out[["recLocs"]] <- asf_stn
    } else if (rec %in% c("esrf_g","esrf_l")) {
      out[["recLocs"]] <- recs
      out[["recPoly"]] <- ifelse(rec == "esrf_grid", recPoly_sf, NULL)
    } 
    
    out[["sobi.box"]] <- c(980,1030,1230,1275)
    if(esrf) out[["esrf"]] <- readRDS(file.path(polygons, "NLpoly.RDS"))
    
    out[["ocean"]] <- ocean
    out[["prj"]] <- prj
    
    return(out)
  }
