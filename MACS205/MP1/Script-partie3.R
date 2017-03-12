#' ---
#' title: "Mini-projet - Partie 3"
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

#' **3.1**

t = 2^(-3)
delta = 1/2
n = 10
estMyst0 = richardson (evalBoiteNoire,n,t,delta)
estMyst0

#' La meilleur estimation de `BoiteNoire` en 0 est donc de $3{,}141684$.
#'
#' **3.2.a**
#' 
#' On utilise la fonction `romberg`, légèrement modifiée par rapport à celle écrite dans le TP4 pour qu'elle renvoie aussi la liste des valeurs "naïves"  de l'intégrale (i.e. sans méthode de Romberg) :
M = 2
n = 15
temp = romberg (evalBoiteNoire, n, a, b, M)
estRomberg = temp[1,]; Inaif = temp[2,]
estRomberg

#' **3.2.b**
#' 
#' On trace les deux erreurs logarithmiques :
Iauto = 4.55676
plot(0:n, log(abs(estRomberg - Iauto)), type = 'l', col = 'red', ylab = 'erreur logarithmique')
lines(0:n, log(abs(Inaif - Iauto)), type = 'l', col = 'blue')

#' On constate que les deux erreurs sont assez similaires et même identiques pour les abscisses 10 à 14, ce qui forme un palier.
#' Ceci peut s'expliquer par le fait que l'on arrive à des valeurs de $M$ telles que chaque point disponible ou presque est utilisé.
#' Le gain dû à la méthode de Romberg n'est alors plus significatif.