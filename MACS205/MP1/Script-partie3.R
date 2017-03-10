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

#' **3.2.a**
#' 
#' On utilise la fonction romberg écrite dans le TP4 :
M = 2
n = 15
estRomberg = romberg (evalBoiteNoire, n, a, b, M)
estRomberg

#' **3.2.b**