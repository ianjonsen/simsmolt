#' @title growth (energetics) function
#' 
#' @description utility function not to be called by user
#' @param w mass at previous time step in g
#' @param ts water temperature at previous time step in C
#' @param s active swim speed at previous time step in km/h
#' @importFrom CircStats rwrpcauchy
#' @importFrom stats rweibull
#' @export
#' 
growth <- function(w, ts, s) {
  
  ## fixed values for Bioenergetics component (from Byron et al. 2014 S1, Table A1)
  # consumption
  CTO <- 17
  CQ <- 4
  CK1 <- 0.25
  CTL <- 24
  CTM <- 18
  CK4 <- 0.75
  CA <- 0.12
  CB <- -0.275
  p <- 0.6 # proportion of consumption
  # respiration
  RA <- 3.021
  RQ <- 1.1017
  SDA <- 0.172
  RTO <- 1.6803
  RB3 <- -0.0068
  RK4 <- 0.0902
  BACT <- 0.9778
  # wastes
  FA <- 0.212
  FB <- -0.222
  FG <- 0.631
  UA <- 0.026
  UB <- 0.580
  UG <- -0.299
  # energy density
  alpha <- 6647
  beta <- 0.5249
  EDprey <- 4650
  
  G1 <- (CTO - CQ)^-1 * log(0.98 * (1 - CK1) / (CK1 * 0.02))
  G2 <- (CTL - CTM)^-1 * log(0.98 * (1 - CK4) / (CK4 * 0.02))
  
  L1 <- exp(G1 * (ts - CQ))
  L2 <- exp(G2 * (CTL - ts))
  
  
  ## consumption (Co)
  Cmax <- CA * w^CB
  Ka <- (CK1 * L1)/(1 + CK1 * (L1 - 1)) 
  Kb <- (CK4 * L2)/(1 + CK4 * (L2 - 1))
  fT <- Ka * Kb
  Co <- Cmax * p * fT
  
  
  ## metabolism - note: convert swim speed (s[i-1]) from km/h back to m/s
  SS <- 0.15 * s / 3.6 ## 15% of sustained swim speed (1.6 body-length / s) converted from km/h back to m/s
  R <- RA/24 * 10^-2 * w^(RB3 * ts + RK4 * SS) * RQ^ts * RTO^SS * BACT^(ts * SS)
  ## egestion
  E <- FA/24 * ts^FB * exp(FG * p) * Co ## E = F from Byron et al
  ## excretion
  U <- UA/24 * ts^UB * exp(UG * p) * (Co - E)
  
  
  ## change in energy (energy density of prey = 4650 J)
  CJ <- Co * EDprey
  RJ <- R * 3240 * 4.184
  SJ <- SDA * (Co - E) * EDprey
  FJ <- E * EDprey
  UJ <- U * EDprey
  dJ <- CJ - (RJ + SJ + FJ + UJ)
  
  
  ## energy density of smolt @ w_t-1 
  EDpred <- alpha + beta * w
  ## smolt weight gain/loss
  w_new <- dJ / EDpred * w + w
  
  
  ############## NOT USED
  ## define starting mass (mpar$sm) as 185 g (@ SoBI) - assuming avg forklength = 27 cm (based on Jon Carr's email of 07/10/2020 re: size rng of postsmolts at SoBI: 24 - 30 cm)
  ## use 47.35 g as avg mass of smolts at tagging - based on 17 cm forklength and L - W power relationship of fL = (m/8987.9)^(1/2.9639)
  ## from Byron et al. 2014 Fish Oceanogr
  
  #########################
  
  return(w_new)
}
               