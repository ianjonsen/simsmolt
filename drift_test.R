
## simulate migration from Campbellton River, NS
require(tidyverse, quietly = TRUE)

d <- suppressMessages(sim_setup(config = "~/OneDrive - Macquarie University/collab/otn/esrf/fn/config.R"))


## simulate V13 tags from LeHave River, NS
## expected tag lifespan 444 d - 25 d before ocean migration starts (guess)
## Ocean migration start date is first day at start location (spars x,y) >= 7 C
## add 20 d to give 21 d min smoltification
start.dt <- lubridate::ymd("2021-07-07", tz = "UTC")

n <- 10
#tmp <- approx(x=c(925, 1200), y=c(1500, 1600), n = n)
#st.loc <- data.frame(x = tmp$x, y = tmp$y)

st.loc <- data.frame(x = rep(1450,n), y = rep(825, n))
st.dt <- seq(start.dt, by=86400, length = n)

sim <- lapply(1:n, function(i) {
  sim_kelt(
    id = i,
    data = d,
    river = "drift",
    mpar = sim_par(
      growth = FALSE,
      scenario = NULL,
      shelf = FALSE,
      N = 150 * 24,
      pN = 1,
      start = c(st.loc$x[i], st.loc$y[i]),
      start.dt = st.dt[i],
      coa = NULL,
      surv = 1,
      Dreten = 0,
      buffer = 5),
    pb = TRUE)
})

run <- tibble(id = 1:n, rep=sim)
class(run) <- append(class(run), "simsmolt", 0)

m <- mapsim(run, data = d, lwd = 0.75)  #, xlim=c(500,1500), ylim=c(300,1500)

print(m)
