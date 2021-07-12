
## simulate migration from Campbellton River, NS
require(tidyverse, quietly = TRUE)
require(simsmolt)

d <- sim_setup()
spars <- suppressMessages(readRDS("~/OneDrive - Macquarie University/collab/otn/esrf/data/kelt_sim_pars.RDS") %>%
                            filter(river == "StMarys"))
start.xy <- spars %>%
  sf::st_as_sf(., coords = c("x","y"), crs = "+proj=laea +lon_0=-71 +lat_0=41 +units=km +datum=WGS84 +no_defs") %>%
  sf::st_transform(., crs = d$prj) %>%
  sf::st_coordinates(.) %>%
  as.data.frame() %>% rename(x=X,y=Y)

## simulate V13 tags from St Marys River, NS
## expected tag lifespan 444 d - 25 d before ocean migration starts (guess)
## Ocean migration start date is first day at start location (spars x,y) >= 7 C
## add 20 d to give 21 d min smoltification
start.dt <- lubridate::as_datetime(which(raster::extract(d$ts, cbind(start.xy$x+20, start.xy$y), method = "bilinear") - 273
                                         >= 7)[1] * 86400,
                                   origin = "2021-01-01 00:00:00",
                                   tz = "UTC")
## calculate number of days b/w tagging & start of ocean migration
om.delay <-
  round(as.numeric(start.dt - lubridate::as_datetime(spars$dt)), 0)

n <- 5
coa1 <- data.frame(x = runif(n,7865,7970), y = runif(n,2030,2190))
coa2 <- data.frame(x = runif(n,6605,7085), y = runif(n,3425,3775))        ## for "as" scenarios
coa3 <- data.frame(x = runif(n,7820,7910), y = runif(n,1860,1925))

#coa1 <- data.frame(x = rnorm(1,900,20), y = rnorm(1,375,20))      ## for "rs" scenarios

sim <- lapply(1:n, function(i) {
  sim_kelt(
      id = i,
      data = d,
      river = "stm",
      mpar = with(subset(spars, river == "StMarys"),
                  sim_par(growth=FALSE,
                          scenario = "as",
                          shelf = FALSE,
                          N = (100 - om.delay) * 24,
                          pN = 1, #0.75
                          start = c(start.xy$x+20,start.xy$y),
                          start.dt = start.dt,
                          coa = cbind(c(coa1$x[i], coa2$x[i], start.xy$x+20),
                                      c(coa1$y[i], coa2$y[i], start.xy$y)), ## for "as" cbind(c(coa1$x[i],x-10), c(coa1$y[i],y-40)),
                          r = 0.001,#0.0005, #0.001, #0.00125,
                          w0 = s.w0 + om.delay, # guess 1 g/d growth after tagging
                          uvm = 1,
                          b = 1.6,
                          rho = 0.6,
                          surv = 1, #0.9968, # increase by 0.0032 1 - (1 - 0.9936)/2 for kelts
                          Dreten = 0 # assume no tag expulsion
                          )), # = p(0.5) @ 500 m  + < 0.01 @ 1000 m   V13

      pb = TRUE) %>%
        sim_detect(., data = d, delay = c(70,150), noise = 0.5)
})

run <- tibble(id = 1:n, rep=sim)
class(run) <- append(class(run), "simsmolt", 0)

m <- mapsim(run, data = d, lwd = 0.75)  #, xlim=c(500,1500), ylim=c(300,1500)

print(m)
