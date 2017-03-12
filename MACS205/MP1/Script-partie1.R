#' ---
#' title: "Mini-projet - Partie 1"
#' output: html_document
#' author: Régis Gourdel
#' ---
#'
#' *Initialisation de l'environnement :*
rm(list = ls())
load('absc-Val.Rdata')
source('fonctions-miniprojet.R')
source('mesfonctions.R')
a = 2^(-16); b = 1.5

#' **1.1**
#' 
#' On réalise une fonction pour tracer les polynômes d'interpolation quatre par quatre, les lignes de plus haut degré étant les plus foncées :
neval = 1000
z = seq(a,b,length.out = neval)
y = evalBoiteNoire(z)
cols = c("red","orange","green","blue")

intpol = function(deg) {
    plot(z, y, type = 'l', col = 'black')
    for (i in 1:length(deg)) {
        x = seq(a, b, length.out = deg[i])
        f = interpolDividif(x, evalBoiteNoire(x), z)
        lines(z, f, col = cols[i], lty = 2)
    }
}

#' Au bas degré, de 2 à 5, on constate que l'interpolation suit assez mal la valeur réelle.
#' Il n'y a donc pas encore de phénomène de Runge mais l'approximation semble assez grossière.
intpol(c(2,3,4,5))

#' À partir du degré 7, en bleu ci-dessous, on commence à voir apparaître un phénomène de Runge.
#' Le polynôme s'écarte fortement de la courbe sur la fin de l'intervalle.
#' Mais au degré 6, en rouge, l'approximation demeure grossière.
intpol(c(6,7,8,9))

#' Pour des degré plus élevés, jusqu'au 30 ci-dessous, on constate que l'approximation est très fine au centre de l'intervalle mais le phénomène de Runge est très important sur les bords.
intpol(c(10,15,20,30))

#' **1.2**
#' 
#' On crée une fonction qui réalise l'interpolation avec les points de Tchebychev :
intpoltche = function(deg) {
    plot(z, y, type = 'l', col = 'black')
    for (i in 1:length(deg)) {
        m = deg[i]
        rtche = function(x) { cos((2*x - 1)*pi/2/m) }
        x = rtche(seq(1:m)) * (b - a)/2 + (a + b)/2
        f = interpolDividif(x, evalBoiteNoire(x), z)
        lines(z, f, col = cols[i], lty = 2)
    }
}

#' Pour des degrés inférieurs à 18 l'approximation s'améliore mais n'est pas encore très fine :
intpoltche(c(4,8,12,18))

#' On arrive alors sur une plage de degrés, de 20 à 42, que l'on pourra qualifier de raisonnables, et qui fournissent une bonne approximation de la fonction :
intpoltche(c(20,26,34,42))

#' Puis à partir du degré 43 on observe un phénomène de Runge dont l'effet s'accentue pour les degrés encore supérieurs :
intpoltche(c(43,47))

#' **1.3.a**
#' 
#' On emploie la fonction piecewiseInterpolation réalisée dans le TD2 :

M = 10
for (n in c(1,2,3,30)) {
    piecewiseInterpol(n, M, a, b, neval/10, nodes = "equi", evalBoiteNoire, TRUE)
}

#' On constate qu'il est plus intéressant de travailler avec des petits degrés, comme $n = 4$.
#' 
#' 
#' **1.3.b**
#' 
#' On code les données et la fonction utiles pour les deux questions qui suivent :
myfun = function(x) { 1 + (sin(3.4 * x) + sin(6 * x)) * (cos(6 * x) + 1) }
n = 3; neval = 1000
M = seq(2,102,5)

myst = function(FUN, M) {
  maxErr = rep(0,length(M))
  for (i in 1:length(M)) {
    temp = piecewiseInterpol(n, M[i], a, b, ceiling(neval/M[i]), "equi", FUN, FALSE)
    yy = temp[1,]; xx = temp[2,]
    estimErr = max( abs(yy - FUN(xx)) )
    maxErr[i] = estimErr
  }
  plot(log(M), log(maxErr), type = 'p')
  return( log(maxErr) )
}

#' On applique notre fonction avec \texttt{myfun} :
logErr = myst(myfun, M)

#' et on obtient la pente sur le première partie linéaire avec :
lm( logErr[1:12] ~ log(M[1:12]) )
#' La pente est donc ici de $-3{,}762$.
#'
#' **1.3.c**
#' 
logErr = myst(evalBoiteNoire, M)
lm( logErr[2:6] ~ log(M[2:6]) )
##' On retrouve une pente proche de $-3{,}745$.