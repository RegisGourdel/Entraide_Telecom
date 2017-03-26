#' ---
#' title: "MACS205 :TD Équations différentielles"
#' output: html_document
#' author: Régis Gourdel
#' ---

library(deSolve)

#' **Question 1**
#' 
#' *a/*
N = function(N0,b,t) { N0 * exp(b*t) }

#' *b/*
N(1, 0.2, 1000)

#' **Question 4**
#' 
#' On définit les paramètres pour toute la suite du problème :
#' 
pars = c(rs = 0.15,
         Ms = 800,
         rc = 0.18,
         Mc = 1200,
         ds = 1,
         dc = 3)

aux_gompertz = function(r,M,N,n) { r * N * log(M/max(1, N + n)) }

gompertz_competition = function(t, N, params) {
  ds = aux_gompertz(params[1], params[2], N[1], N[2])
  dc = aux_gompertz(params[3], params[4], N[2], N[1])
  return( list(c(ds,dc)) )
}

Nini  = c(Ns = 1000, Nc = 10)
times = seq(0, 200, by = 1)

out   = ode(Nini, times, gompertz_competition, pars)
summary(out)
plot(out)


#' **Question 5**

euler = function(yini, times, func, pars) {
  n = length(times)
  y1 = c(yini[1])
  y2 = c(yini[2])
  for (i in 2:n) {
    d  = unlist(func(times[i-1], c(y1[i-1], y2[i-1]), pars))
    dt = times[i] - times[i-1]
    y1 = c(y1, y1[i-1] + dt*d[1])
    y2 = c(y2, y2[i-1] + dt*d[2])
  }
  return( rbind(y1,y2) )
}

outEuler = euler(c(1000,10), times, gompertz_competition, pars)
NsEuler = outEuler[1,]
NcEuler = outEuler[2,]

plot(times, NsEuler, type = 'l', col = "blue", ylab = "Population de cellules saines")
plot(times, NcEuler, type = 'l', col = "orange", ylab = "Population de cellules cancéreuse")

#' **Question 6**
#' 
#' Avec les paramètres donnés on constate que la population de cellules saines tend à disparaître alors que le nombre de cellules cancéreurses grandit jusqu'à une certaine limite.
#'
#'
#' **Question 7**
#' 

injection = function(t) {
    if (t == 30 || t == 31 || t == 55 || t == 56 || t == 80 || t == 81)
        return(1)
    else
        return(0)
}

gompertz_chimiotherapie = function(t, N, params) {
    temps = aux_gompertz(params[1], params[2], N[1], N[2]) - ds * N[1] * injection(t)
    tempc = aux_gompertz(params[3], params[4], N[2], N[1]) - dc * N[2] * injection(t)
    return( list(c(temps,tempc)) )
}

outChimio = euler(c(1000,10), times, gompertz_chimiotherapie, pars)
NsChimio = outChimio[1,]
NcChimio = outChimio[2,]

plot(times, NsChimio, type = 'l', col = "blue", ylab = "Population de cellules saines")
plot(times, NcChimio, type = 'l', col = "orange", ylab = "Population de cellules cancéreuse")