#' @title Pre-simulation setup
#'
#' @description load required rasters, receiver locations
#'
#'
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#'
#' @param config - path to config.R script containing file.paths for required & optional environmental layers (see Details)
#' @importFrom raster raster stack brick projectRaster extract
#' @importFrom sp coordinates<- proj4string<- CRS spTransform SpatialPointsDataFrame spsample
#' @importFrom sf st_as_sf st_sample st_coordinates st_distance
#' @importFrom dplyr select filter rename bind_cols %>% tibble distinct
#' @export
#'
sim_setup <-
  function(config = file.path("/", "Users", "jonsen", "OneDrive - Macquarie University", "collab", "otn", "esrf", "fn", "config.R")) {

    ## FIXME: this needs to be generalized - provide spatial extent for query to download ETOPO2 data?
    ## FIXME:   or rely on user supplying their own bathymetry data

    ## load receiver location data

      ## FIXME: this needs to be generalized - do receiver data munging prior to using this function...
      ## FIXME:  could generalize by accessing OTN server to pull requested receiver data from anywhere...
      ## FIXME:  prep code would prob require consistent receiver location / history format on OTN server

    source(config)
    if(is.null(prj)) prj <- "+proj=laea +lat_0=41 +lon_0=-71 +units=km +datum=WGS84"

    out <- list(
      bathy = raster(bathy),
      land = raster(d2land),
      land_dir = raster(land_dir)
    )

    if (!is.null(d2shelf))
      out[["shelf"]] <- stack(d2shelf)

    out[["u"]] <- stack(file.path(riops, "riops_doy_u.grd"))
    out[["v"]] <- stack(file.path(riops, "riops_doy_v.grd"))
    out[["ts"]] <- stack(file.path(riops, "riops_doy_t.grd"))

    out[["recLocs"]] <- esrf_rec
    out[["recPoly"]] <- recPoly_sf

    out[["sobi.box"]] <- c(980,1030,1230,1275)
    out[["esrfPoly"]] <- readRDS(file.path(polygon.file, "NLpoly.RDS"))

    out[["prj"]] <- prj

    return(out)
  }
