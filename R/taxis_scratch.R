##taxis bits

## employ rheotaxis or not by adjusting move direction (default is a bias to North)
## the strength of adherence to taxis is dictated by rho: rho = 1 - perfect taxis, rho = 0 - random movement, no taxis
## magnitude of current relative to magnitude of swimming, normalized to 1
# switch(mpar$taxis,
#        pos = {
#          ## swim against current
#          mu_c[i] <- (atan2(u[i], v[i]) + pi) %% (2 * pi)
#        },
#        neg = {
#          ## swim with current
#          mu_c[i] <- atan2(u[i], v[i]) %% (2 * pi)
#        },
#        none = {
#          ## current direction does not influence active swimming direction
#          mu_c[i] <- NA
#        })

## FIXME: Rheotaxis not really working as intended - maybe interaction with biased movement toward coa is messing it up
## try coding so currents and coa effects on movement can be turned on/off independently to investigate this



## I DON'T THINK THIS IS USEFUL AS IT BLURS ACTIVE SWIM DIRECITON BIAS WITH CURRENT DIRECTION, WHICH IS ALREADY BUILT INTO u,v
## calculate weighted circular mean, using current & swimming magnitudes as weights
# uv <- sqrt(u[i]^2 + v[i]^2)
# uv.w <- uv / (uv + s[i])
# s.w <- s[i] / (uv + s[i])
# 
# if(!is.na(mu_c[i])) {
#   sinr <- (sin(mu_c[i]) * uv.w + sin(mu_s[i]) * s.w)
#   cosr <- (cos(mu_c[i]) * uv.w + cos(mu_s[i]) * s.w)
#   mu[i] <- atan2(sinr, cosr)
# } else {
#   mu[i] <- mu_s[i]
# }
