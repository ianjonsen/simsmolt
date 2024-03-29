% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_par.R
\name{sim_par}
\alias{sim_par}
\title{Control Values for \code{simulate}.}
\usage{
sim_par(
  move = c("brw", "rw", "drift"),
  temp = TRUE,
  advect = TRUE,
  noise = NULL,
  growth = TRUE,
  shelf = TRUE,
  scenario = "rs",
  land = FALSE,
  boundary = FALSE,
  ...
)
}
\arguments{
\item{move}{"brw" a biased RW; "rw" simple RW; "drift" no active swimming}

\item{temp}{logical}

\item{advect}{logical}

\item{noise}{numeric; add small random component to u,v vectors - only applies when mode = "drift"}

\item{growth}{logical}

\item{shelf}{logical; should smolts be constrained to stay in water > - 1000 m depth (continental shelf)}

\item{land}{keep track of sim rep hitting land (TRUE)}

\item{boundary}{keep track of sim rep hitting sim boundary (TRUE)}

\item{...}{additional simulation control parameters}

\item{migs}{migration scenario 1 or 2 (stop migration upon arrival to Grand Banks)}
}
\value{
Returns a list with components
  \item{\code{move}}{the main move process}
  \item{\code{temp}}{temperature-dependent movements}
  \item{\code{advect}}{ocean-current-dependent movements}
  \item{\code{noise}}{add random noise to current u,v values}
  \item{\code{growth}}{temperature-dependent growth}
  \item{\code{shelf}}{movements constrained to stay on shelf}
  \item{\code{migs}}{migration strategy}
  \item{\code{land}}{record if smolt gets "stuck" on land, thereby ending simulation}
  \item{\code{boundary}}{record if smolt hits simulation boundary, thereby ending simulation}
  \item{\code{pars}}{list of additional, required control parameters}
}
\description{
\code{sim_par} defines the simulation parameters & control scenarios used by \code{simulate}.
}
\details{
The movement process used predominantly in the simulation is
selected by the \code{move} argument.  Additional
parameters include: \code{temp} is movement temperature-dependent; \code{advect} do ocean
currents influence movement; \code{growth} do smolts grow in temp-dependent fashion;
\code{start.date};
\code{start} (location); \code{coa} centre of attraction (can be NULL);
\code{mdir} directional bias (a 2-element vector);
\code{rho} concentration of step directions (for wrapped-Cauchy, a 2-element vector); ...
}
