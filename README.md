# simsmolt
Simulate salmon smolt migration and potential for acoustic detection

### Roadmap: 
while focused on a specific project, this code will be generalised so detections of acoustically-tagged animals can be simulated under a variety of ecological scenarios

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.org/ianjonsen/simsmolt.svg?branch=master)](https://travis-ci.org/ianjonsen/simsmolt)

## tidy simulation
### simulate a single smolt's migration for 75 d (1800 h)
`d <- sim_setup()`  
`out <- sim_move(N=1800, data=d) %>% sim_detect(., data = d)`  
`summary(out)`  
`plot(out)`

### simulate multiple, independent smolts
`d <- sim_setup()`  
`out <- data.frame(id=1:5) %>%`  
    `group_by(id) %>%`  
    `do(rep = try(sim_move(N=1200, data = d, mpar=list(surv=1)) %>% sim_detect(., data  = d)))`  
`class(out) <- append(class(out), "simsmolt", 0)`  
`plot(out)`

### simulate multiple smolts across multiple processor cores
`d <- sim_setup()`  
`cls <- multidplyr::create_cluster(11)`  
`multidplyr::cluster_library(cls, c("simsmolt", "dplyr"))`  
`multidplyr::cluster_copy(cls, d)`  
`out <- data.frame(id = 1:100) %>%`  
    `multidplyr::partition(id, cluster = cls) %>%`  
    `do(rep = try(sim_move(N=1200, data = d, mpar=list(rho=0.8), pb=FALSE) %>% sim_detect(., data = d))) %>%`  
    `collect()`  
`class(out) <- append(class(out), "simsmolt", 0)`  
`plot(out)`  
`summary(out)`
