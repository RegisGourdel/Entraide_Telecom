#' ---
#' title: "Mini-projet"
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