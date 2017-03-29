#' ---
#' title: "MACS205 : TD 1 et 2 Équations différentielles"
#' output: html_document
#' author: Régis Gourdel
#' ---
#'
#'
#' Pour exécuter ce script, ouvrez le avec RStudio et pressez Ctr+Shift+Entrée.
#' Pour l'exécuter et le convertir (compiler)
#' en notebook html, pressez Ctrl+shift+K (il faut avoir installé le package `rmarkdown`).
#' Pour obtenir un pdf plutôt qu'un html (déconseillé), remplacez
#' `output: html_document` par `output: pdf_document` dans l'en-tête ci-dessus.
#'
#' **TD 1**
#'
#' **Exercice 3**
#'
#' *d/*
LVmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Ingestion    <- rIng  * Prey * Predator
    GrowthPrey   <- rGrow * Prey * (1 - Prey/K)
    MortPredator <- rMort * Predator
    
    dPrey        <- GrowthPrey - Ingestion
    dPredator    <- Ingestion * assEff - MortPredator
    
    return(list(c(dPrey, dPredator)))
  })
}

pars  <- c(rIng   = 0.2,    # /day, rate of ingestion
           rGrow  = 1.0,    # /day, growth rate of prey
           rMort  = 0.2 ,   # /day, mortality rate of predator
           assEff = 0.5,    # -, assimilation efficiency
           K      = 10)     # mmol/m3, carrying capacity

yini  <- c(Prey = 1, Predator = 2)
times <- seq(0, 200, by = 1)
out   <- deSolve::ode(yini, times, LVmod, pars)
summary(out)
plot(out)