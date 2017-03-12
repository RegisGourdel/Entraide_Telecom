#' ---
#' title: "Mini-projet - Partie 2"
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

#' **2.1.a**
#' 
M = seq(2,50,1)
simp = sapply(M, function(m) { simpsonInt(evalBoiteNoire, a, b, m) } )
plot(M, simp, type = 'l')

M = seq(10,50,1)
simp = sapply(M, function(m) { simpsonInt(evalBoiteNoire, a, b, m) } )
plot(M, simp, type = 'l')

M = seq(40,80,1)
simp = sapply(M, function(m) { simpsonInt(evalBoiteNoire, a, b, m) } )
plot(M, simp, type = 'l')

#' **2.1.b**
#'
#' On a donc, à $10^{-3}$ près, $I \simeq 4{,}936$.
#' 
#' 
#' **2.2.a**
#' 
#' D'après la proposition 3.5.1, l'ordre de grandeur théorique de l'erreur en fonction de $h$ est $h^{N + 1}$ où $N$ est l'ordre de la méthode.
#' Or ici la méthode de Cavalieri-Simpson est de rang 2, donc d'ordre $N = 3$.
#' 
#' On en déduit l'erreur $EI_{3,M}^{(a,b)} \propto h^4$.
#' 
#' **2.2.b**
#' 
#' On a, d'après la proposition 3.5.1 et le résultat précédent :
#' $$\log \left\vert I - \hat{I}_M^{simp} \right\vert = \log \left\vert EI_{3,M}^{(a,b)} \right\vert = \log \left\vert C^{te}_1 \cdot h^4 \right\vert = C^{te}_2 + 4 \cdot \log(h) = C^{te}_2 + 4 \cdot \log \left( \left( \frac{b - a}{M} \right)^4 \right) = C^{te}_3 - 4 \cdot \log(M)$$
#' d'où un graphe d'erreur log-log avec l'allure suivante :
M = seq(1,80)
plot(log(M), 6 - 4*log(M), type = 'l', main = "log erreur en fonction de log(M)", ylab = "log(Errr)")

#' **2.2.c**
M = seq(7,40)
vectErrorSimpson = rep(0, length(M))
for (i in 1:length(M)) {
  vectErrorSimpson[i] = abs( evalErrSimpson(evalBoiteNoire, a, b, M[i])[1] )
}
logM = log(2*M)
logErr = log(vectErrorSimpson)
plot(logM, logErr, type = 'p', main = 'log erreur a posteriori en fonction de log(M)')

#' On effectue une régression linéaire avec la fonction suivante :
lm( logErr[1:10] ~ logM[1:10] )

#' On trouve donc une pente proche de 4, ce qui est donc cohérent avec le résultat prévu.
#' 
#' Le phénomène observé pour $M$ grand peut s'expliquer par le fait que...
#' TODO
#' 
#' **2.3.a**
#' 
M = 10 ; erreur = 1 ; Mvecteur = 0 ; Qvecteur = 0 ; Evecteur = 0
tolerance = 10^(-5)
while (2*erreur > tolerance) {
  temp = evalErrSimpson(evalBoiteNoire, a, b, M)
  erreur = abs(temp[1])
  Evecteur = c(Evecteur, erreur)
  Qvecteur = c(Qvecteur, temp[2])
  Mvecteur = c(Mvecteur, M)
  M = 2*M
}

Evecteur
Qvecteur
Mvecteur

#' On déduit donc de ces résultats que $I = I_{auto} = 4{,}55676$ à $10^{-5}$ près.
#' Le nombre d'intervalles correspondannt est $M_{auto} = 20$.

Iauto = Qvecteur[3]
Mauto = Mvecteur[3]
#' 
#' **4.a**
#' 
Mvecteur = seq(1,40)
Qvecteur = rep(0,length(Mvecteur))
for (i in Mvecteur) {
  Qvecteur = simpsonInt(evalBoiteNoire, a, b, i)
}
plot(log(Mvecteur), log(abs(Qvecteur - Iauto)))