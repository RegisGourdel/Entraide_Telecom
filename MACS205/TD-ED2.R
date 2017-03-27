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
parms = c(rs = 0.15,
         Ms = 800,
         ds = 1,
         rc = 0.18,
         Mc = 1200,
         dc = 3)

Nini  = c(Ns = 1000, Nc = 10)

times = seq(0, 150, length.out = 1000)

#' Ainsi que la fonction qui donne la dérivée:

gompertz_competition = function(t,N,params){
    with(as.list(c(N,parms)), {
        dNs = rs * N[1] * log(Ms / max(1,N[1] + N[2]))
        dNc = rc * N[2] * log(Mc / max(1,N[1] + N[2]))
        return( list(c(dNs,dNc)) )
    })
}

#' et l'on fait appel à ``ode``:

out = ode(Nini, times, gompertz_competition, pars)
summary(out)
plot(out)


#' **Question 5**

euler = function(yini, times, func, pars) {
  n = length(times)
  y1 = c(yini[1])
  y2 = c(yini[2])
  for (i in 2:n) {
    d  = unlist( func(times[i-1], c(y1[i-1], y2[i-1]), pars) )
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
lines(times, NcEuler, type = 'l', col = "orange", ylab = "Population de cellules cancéreuse")

#' **Question 6**
#' 
#' Avec les paramètres donnés on constate que la population de cellules saines tend à disparaître alors que le nombre de cellules cancéreurses grandit jusqu'à une certaine limite.
#'
#'
#' **Question 7**
#' 

injection = function(t) {
    if ((t >= 30 && t < 31) || (t >= 55 && t < 56) || (t >= 80 && t < 81))
        return(1)
    else
        return(0)
}

gompertz_chimiotherapie = function(t, N, params) {
    with(as.list(c(N,parms)), {
        dNs = rs * N[1] * log(Ms / max(1,N[1] + N[2])) - ds * N[1] * injection(t)
        dNc = rc * N[2] * log(Mc / max(1,N[1] + N[2])) - dc * N[2] * injection(t)
        return( list(c(dNs,dNc)) )
    })
}

plotChimio = function(pas) {
    times = seq(0, 130, by = pas)
    outChimio = euler(Nini, times, gompertz_chimiotherapie, pars)
    NsChimio = outChimio[1,]
    NcChimio = outChimio[2,]
    c_range = range(0, NsChimio, NcChimio)
    plot(times, NsChimio, type = 'l', col = "blue", ylim = c_range)
    lines(times, NcChimio, type = 'l', col = "orange")
}

#' On teste alors la modélisation de la chimio en affinant le pas:
#' 
plotChimio(1)
plotChimio(0.5)
plotChimio(0.1)
plotChimio(0.01)
#plotChimio(0.001)

