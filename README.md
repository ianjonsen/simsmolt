# simsmolt
Simulate salmon smolt migration and potential for acoustic detection

### Roadmap: 
while focused on a specific project, this code will be generalised so detections of acoustically-tagged animals can be simulated under a variety of ecological scenarios

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.org/ianjonsen/simsmolt.svg?branch=master)](https://travis-ci.org/ianjonsen/simsmolt)

## tidy simulation
### simulate a single smolt's migration for 75 d (1800 h)
`out <- sim_setup() %>% sim_move(N=1800, data=.) %>% sim_detect()`  
`summary(out)`  
`plot(out)`

### simulate multiple, independent smolts
`d <- sim_setup()`  
`out <- data.frame(id=1:5) %>% group_by(id) %>% do(rep = sim_move(N=1200, data = d, mpar=list(surv=1)) %>% sim_detect(.))`

### simulate multiple smolts across multiple processor cores
`d <- sim_setup()`  
`cls <- multidplyr::create_cluster(11)`  
`out <- data.frame(id = 1:100) %>%`  
  `multidplyr::partition(id, cluster = cls %>%`  
  `do(rep = sim_move(N=1200, data = d, mpar=list(surv=1, rho=0.8), pb=FALSE) %>% sim_detect(.)) %>%`  
  `collect()`
