# simsmolt
Simulate salmon smolt migration and potential for acoustic detection

### Roadmap: 
while focused on a specific project, this code will be generalised so detections of acoustically-tagged animals can be simulated under a variety of ecological scenarios

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.org/ianjonsen/simsmolt.svg?branch=master)](https://travis-ci.org/ianjonsen/simsmolt)

## tidy simulation
`out <- sim_setup() %>% sim_move(data=.) %>% sim_detect()`  
`summary(out)`  
`plot(out)`
