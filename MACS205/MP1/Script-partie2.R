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
plot(log(M), 6 - 4*log(M), type = 'l')