#' ---
#' title: "MACS205 : TD Équations différentielles"
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
pars = c(rs = 0.15,
         Ms = 800,
         rc = 0.18,
         Mc = 1200)

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

euler = function(Fini, times, func, pars) {
  n = length(times)
  values = rep(Fini, n)
  for (i in times[2:n]) {
    values[i] = values[i-1] + (times[i] - times[i-1]) * unlist(func(times[i-1], values[i-1], pars))
  }
  return(values)
}

out2 = euler(Nini, times, gompertz_competition, pars)
summary(out2)
plot(out2)