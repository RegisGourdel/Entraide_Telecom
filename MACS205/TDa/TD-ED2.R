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

#' Ainsi que la fonction qui donne la dérivée:

gompertz_competition = function(t,N,params){
    with(as.list(c(N,parms)), {
        dNs = rs * N[1] * log(Ms / max(1,N[1] + N[2]))
        dNc = rc * N[2] * log(Mc / max(1,N[1] + N[2]))
        return( list(c(dNs,dNc)) )
    })
}

#' et l'on fait appel à ``ode``:

times = seq(0, 150, length.out = 1000)
out = ode(Nini, times, gompertz_competition, pars)
plot(out)


#' **Question 5**
#' 
#' On réalise la fonction ``euler`` pour résoudre le système différentiel :
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

#' On teste la fonction avec ``gompertz_competition`` :
outEuler = euler(c(1000,10), times, gompertz_competition, pars)
NsEuler = outEuler[1,]
NcEuler = outEuler[2,]

plot(times, NcEuler, type = 'l', col = "orange", ylab = "Population de cellules cancéreuse")
lines(times, NsEuler, type = 'l', col = "blue", ylab = "Population de cellules saines")

#' **Question 6**
#' 
#' Avec les paramètres donnés on constate que la population de cellules saines tend à disparaître alors que le nombre de cellules cancéreuses grandit jusqu'à une certaine limite.
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

#' On crée une fonction pour automatiser le tracer des courbes :
plotChimio = function(func,times) {
    outChimio = func(Nini, times, gompertz_chimiotherapie, pars)
    NsChimio = outChimio[1,]
    NcChimio = outChimio[2,]
    c_range = range(0, NsChimio, NcChimio)
    plot(times, NsChimio, type = 'l', col = "blue", ylim = c_range)
    lines(times, NcChimio, type = 'l', col = "orange")
}

#' On teste alors la modélisation de la chimio en affinant le pas:
#' 
#' La première, avec un pas de 1 dévie fortement vers $-\infty$ pour les cellules cancéreuses :
plotChimio(euler, seq(0, 60, by = 1))

#' La seconde, avec un pas de 0.5 donne un résultat qui n'est pas correct, et dévie également ensuite :
plotChimio(euler, seq(0, 100, by = 0.5))

#' Les deux suivantes fournissent un tracé plus cohérent :

plotChimio(euler, seq(0, 130, by = 0.1))
plotChimio(euler, seq(0, 130, by = 0.01))

#' Un pas de 0.001 est trop petit pour être exécuté dans un temps raisonnable donc nous nous contentons de ceux plus grands.
#' 
#' 
#' On reproduit maintenant cela avec la fonction ``ode`` :
outR2 = ode(Nini, times = seq(0, 130, by = 0.5), gompertz_chimiotherapie, pars)
plot(outR2)

outR3 = ode(Nini, times = seq(0, 130, by = 0.1), gompertz_chimiotherapie, pars)
plot(outR3)

#' On constate que les résultats sont similaires.
#' 
#' 

#' **Question 8**
#' 
#' On code une fonction pour appliquer la méthode de Runge-Kutta avec $q = 4$ :
RK4 = function(yini, times, func, pars) {
    n = length(times)
    y1 = c(yini[1])
    y2 = c(yini[2])
    a2 = 1/2; a3 = 1/2; a4 = 1
    b1 = 1/6; b2 = 2/6; b3 = 3/6; b4 = 4/4
    for (i in 2:n) {
        y = c(y1[i-1], y2[i-1])
        h   = times[i] - times[i-1]
        k1  = unlist( func(times[i-1], y, pars) )
        k2  = unlist( func(times[i-1] + a2*h , y + a2*k1*h, pars) )
        k3  = unlist( func(times[i-1] + a3*h , y + a3*k2*h, pars) )
        k4  = unlist( func(times[i-1] + a4*h , y + a4*k3*h, pars) )
        phi = b1*k1 + b2*k2 + b3*k3 + b4*k4
        y1  = c(y1, y1[i-1] + h * phi[1])
        y2  = c(y2, y2[i-1] + h * phi[2])
    }
    return( rbind(y1,y2) )
}

#' et l'on teste notre fonction pour les pas 0.5, 0.1 et 0.01.
#' 
#' Le premier n'est pas probant, et dévie rapidement vers $-\infty$ pour les cellules cancéreuses sur un intervalle plus grand :
plotChimio(RK4, seq(0, 80, by = 0.5))

#' Les deux suivants fournissent une solution assez bonne :
plotChimio(RK4, seq(0, 150, by = 0.1))
plotChimio(RK4, seq(0, 150, by = 0.01))
#' On est proche ici des résultats précédents. Mais pour un même pas on constate que le résultat avec ``RK4`` semble meilleur que celui avec ``euler``.
#' 
#' 

#' **Question 9**
#' 
#' On constate que c'est au moment des injections que se produisent les changements les plus rapides.
#' On peut donc définir un maillage temporel plus serré sur ces moments ci.
#' 
#' 
#' Le pas de base adopté est 0.5, qui ne fonctionne pas normalement, et l'on divise ce pas par 5 autour des moments d'injection.
pas = 0.5
t1  = seq(0, 29, by = pas)
t2  = seq(29 + pas, 32, by = pas/5)
t3  = seq(32 + pas, 54, by = pas)
t4  = seq(54 + pas, 57, by = pas/5)
t5  = seq(57 + pas, 79, by = pas)
t6  = seq(79 + pas, 82, by = pas/5)
t7  = seq(82 + pas, 120, by = pas)
times = c(t1,t2,t3,t4,t5,t6,t7)

plotChimio(RK4, times)

#' On constate donc que l'on trouve un résultat proche de celui que nous avions pour un pas de 0.1, alors que le nombre de points a été ici nettement réduit, et le maillage temporel n'est affiné qu'en quelques endroits.