##' \code{sim_par} defines the simulation parameters & control scenarios used by \code{simulate}.
##'
##' The movement process used predominantly in the simulation is
##' selected by the \code{move} argument.  Additional 
##' parameters include: \code{temp} is movement temperature-dependent; \code{advect} do ocean
##' currents influence movement; \code{growth} do smolts grow in temp-dependent fashion; 
##' \code{taxis} do smolts active swim with/against current; \code{start.date}; 
##' \code{start} (location); \code{coa} centre of attraction (can be NULL); 
##' \code{mdir} directional bias (a 2-element vector);
##' \code{rho} concentration of step directions (for wrapped-Cauchy, a 2-element vector); ...
##'
##' @title Control Values for \code{simulate}.
##' @param move "brw" a biased RW; "rw" simple RW; "drift" no active swimming
##' @param temp logical
##' @param advect logical
##' @param growth logical
##' @param shelf logical; should smolts be constrained to stay in water > - 1000 m depth (continental shelf)
##' @param taxis "p" (+ve rheotaxis), "n" (-ve rheotaxis), NA (no taxis; default)
##' @param migs migration scenario 1 or 2 (stop migration upon arrival to Grand Banks)
##' @param land keep track of sim rep hitting land (TRUE)
##' @param boundary keep track of sim rep hitting sim boundary (TRUE)
##' @param ... additional simulation control parameters
##' @return Returns a list with components
##'   \item{\code{move}}{the main move process}
##'   \item{\code{temp}}{temperature-dependent movements}
##'   \item{\code{advect}}{ocean-current-dependent movements}
##'   \item{\code{growth}}{temperature-dependent growth}
##'   \item{\code{shelf}}{movements constrained to stay on shelf}
##'   \item{\code{taxis}}{behavioural response to currents}
##'   \item{\code{migs}}{migration strategy}
##'   \item{\code{land}}{record if smolt gets "stuck" on land, thereby ending simulation}
##'   \item{\code{boundary}}{record if smolt hits simulation boundary, thereby ending simulation}
##'   \item{\code{pars}}{list of additional, required control parameters}
##' @export

sim_par <-
  function(move = c("brw","rw","drift"),
           temp = TRUE,
           advect = TRUE,
           growth = TRUE,
           shelf = TRUE,
           taxis = c(NA,"p","n"),
           scenario = 2,
           land = FALSE,
           boundary = FALSE,
           ...) {
    
    move <- match.arg(move)
    taxis <- match.arg(taxis)
    
    dots <- list(...)
    
    pars <- list(
      start.dt = ISOdatetime(2021,07,09,00,00,00, tz = "UTC"),
      start = c(995, 1240),
      coa = NULL,
      mdir = c(75,-45)/180*pi, # bias direction for N migration (S migration is mdir - pi)
      NFline = runif(1, 1350, 1450), # location on x-axis at which smolt turns N to Lab Shelf
      rho = c(0.6, 0.4), # directional persistence for brw [1] and rw [2]
      turn = 2.5, # rate at which smolts turn N after rounding NF
      ntries = 1,
      ts.q = 0.75,
      psi = 0.9,
      uvm = 1, # magnitude of current vectors: if uvm < 1 current strength is down-scaled
      buffer = 5,
      b = 2,
      a = 0,
      w0 = 185,
      surv = 0.9936, ## daily survival rate
      reten = 0.845^(1/60), ## daily V7/8 tag retention rate (in first ~ 60 d - Brundson et al 2019 ICES JMarSci 76:7)
      Dreten = 60, ## number of days within which tags can be expulsed
      pdrf = c(5, -0.02), # = p(0.5) @ 250 m  + < 0.01 @ 500 m   [c(4.865, -0.0139)  (~ consistent w HFX line V9 @ high power)]
      beta = c(-2, -2) # potential fn params to keep smolts on shelf
    )
    
    ## overide default control pars
    pars[names(dots)] <- dots
    
    list(move = move,
         temp = temp,
         advect = advect,
         growth = growth,
         shelf = shelf,
         taxis = taxis,
         scenario = scenario,
         land = land,
         boundary = boundary,
         pars = pars)
  }
