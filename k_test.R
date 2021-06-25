
## simulate migration from Campbellton River, NS
require(tidyverse, quietly = TRUE)

spars <- suppressMessages(readRDS("~/OneDrive - Macquarie University/collab/otn/esrf/data/kelt_sim_pars.RDS") %>%
                            filter(river == "StMarys"))
d <- suppressMessages(sim_setup(config = "~/OneDrive - Macquarie University/collab/otn/esrf/fn/config.R"))


## simulate V13 tags from LeHave River, NS
## expected tag lifespan 444 d - 25 d before ocean migration starts (guess)
## Ocean migration start date is first day at start location (spars x,y) >= 7 C
## add 20 d to give 21 d min smoltification
start.dt <- lubridate::as_datetime(which(extract(d$ts, cbind(spars$x+10, spars$y-5), method = "bilinear") - 273
                                         >= 7)[1] * 86400,
                                   origin = "2021-01-01 00:00:00",
                                   tz = "UTC")
## calculate number of days b/w tagging & start of ocean migration
om.delay <-
  round(as.numeric(start.dt - lubridate::as_datetime(spars$dt)), 0)

n <- 5
coa1 <- data.frame(x = rep(1550,n), y = runif(n,600,775))
coa2 <- data.frame(x = runif(n,1200,1800), y = rep(2450,n))        ## for "as" scenarios
coa3 <- data.frame(x = runif(n,1350,1450), y = 550)

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
                          start = cbind(x+10,y-5),
                          start.dt = start.dt,
                          coa = cbind(c(coa1$x[i], coa2$x[i], x+10),
                                      c(coa1$y[i], coa2$y[i], y-5)), ## for "as" cbind(c(coa1$x[i],x-10), c(coa1$y[i],y-40)),
                          r = 0.00125, #0.001, #0.00125,
                          w0 = s.w0 + om.delay, # guess 1 g/d growth after tagging
                          uvm = 1,
                          b = 1.6,
                          rho = 0.6,
                          surv = 1, #0.9968, # increase by 0.0032 1 - (1 - 0.9936)/2 for kelts
                          Dreten = 0, # assume no tag expulsion
                          buffer = 5,
                          pdrf = c(5, -0.01))), # = p(0.5) @ 500 m  + < 0.01 @ 1000 m   V13

      pb = TRUE) %>%
        sim_detect(., data = d, delay = c(70,150), noise = 0.5)
})

run <- tibble(id = 1:n, rep=sim)
class(run) <- append(class(run), "simsmolt", 0)

m <- mapsim(run, data = d, lwd = 0.75)  #, xlim=c(500,1500), ylim=c(300,1500)

print(m)
