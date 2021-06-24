
## simulate migration from Campbellton River, NS
require(tidyverse, quietly = TRUE)

spars <- suppressMessages(readRDS("~/OneDrive - Macquarie University/collab/otn/esrf/data/kelt_sim_pars.RDS") %>%
                            filter(river == "StMarys"))
d <- suppressMessages(sim_setup(config = "~/OneDrive - Macquarie University/collab/otn/esrf/fn/config.R"))


## simulate V13 tags from LeHave River, NS
## expected tag lifespan 444 d - 25 d before ocean migration starts (guess)
## Ocean migration start date is first day at start location (spars x,y) >= 7 C
## add 20 d to give 21 d min smoltification
start.dt <- lubridate::as_datetime(which(extract(d$ts, cbind(spars$x, spars$y), method = "bilinear") - 273
                                         >= 8)[1] * 86400,
                                   origin = "2021-01-01 00:00:00",
                                   tz = "UTC")
## calculate number of days b/w tagging & start of ocean migration
om.delay <-
  round(as.numeric(start.dt - lubridate::as_datetime(spars$dt)), 0)

#coa1 <- data.frame(x = runif(1,1200,1800), y = rep(2450,1))        ## for "as" scenarios
coa1 <- data.frame(x = rnorm(1,700,20), y = rnorm(1,450,20))     ## for "rs" scenarios

sim <- lapply(1:1, function(i) {
  sim_kelt(
      id = i,
      data = d,
      river = "stm",
      mpar = with(subset(spars, river == "StMarys"),
                  sim_par(growth=FALSE,
                          scenario = "rs",
                          shelf = TRUE, ## but only for land avoidance (set in moveKcam.R) & only for "as"
                          N = (104 - om.delay) * 24,
                          start = cbind(x,y),
                          start.dt = start.dt,
                          mdir = NULL,
                          coa = cbind(c(coa1$x[i], x),c(coa1$y[i],y)), ## for "as" cbind(c(coa1$x[i],x-10), c(coa1$y[i],y-40)),
                          r = 0.00125, #0.001, #0.00125,
                          w0 = s.w0 + om.delay, # guess 1 g/d growth after tagging
                          uvm = 1,
                          b = 1.6,
                          rho = 0.4,
                          surv = 1, #0.9968, # increase by 0.0032 1 - (1 - 0.9936)/2 for kelts
                          Dreten = 0, # assume no tag expulsion
                          buffer = 5,
                          pdrf = c(5, -0.01), # = p(0.5) @ 500 m  + < 0.01 @ 1000 m   V13
                          beta = c(-2,-2))),
      pb = TRUE) %>%
        sim_detect(., data = d, delay = c(70,150), noise = 0.5)
})

run <- tibble(id = 1:1, rep=sim)
class(run) <- append(class(run), "simsmolt", 0)

print(mapsim(run, data = d, lwd = 0.5) + geom_point(data = data.frame(x=coa1$x,y=coa1$y), aes(x,y), shape=17, col="red"))
