#' ---
#' title: "MACS205 : TD1"
#' output: html_document
#' author: Régis Gourdel
#' ---
#'
#'
#' Pour exécuter ce script, ouvrez le avec RStudio et pressez Ctr+Shift+Entrée.
#' Pour l'exécuter et le convertir (compiler)
#' en notebook html, pressez Ctrl+shift+K (il faut avoir installé le package `rmarkdown`).
#' Pour obtenir un pdf plutôt qu'un html (déconseillé), remplacez
#' `output: html_document` par `output: pdf_document` dans l'en-tête ci-dessus.
#'
#' **Exercice 2**
#'
#'On a $\left\{ \begin{array}{ll} I_n = \frac{1}{n} - 5 I_{n - 1}\\ I_0 = \log \left( \frac{6}{5} \right) \end{array}\right.$.
#'En particulier il vient $\lim_n I_n = 0$.
#'
suiteIn = function(n, init){
  output = rep(0,n + 1)
  output[1] = init
  for (k in 2:(n + 1)){
    output[k] = 1/(k - 1) - 5*output[k - 1]
  }
  return(output)
}

u0 = log(6/5);
v0 = u0 + 1e-3;
n = -50;
U = suiteIn(n,u0)
V = suiteIn(n,v0)
plot(1:(n + 1), U, col = 'blue', type = 'l')
lines(1:(n + 1), V, col = 'red')

#'
#'
#' **Exercice 4**
#' 
#' On a le développement limité $e(x) = 1 + x + \frac{x^2}{2} + o(x^2)$, d'où $f(x) = 1 + \frac{x}{2} + o(x)$.
#' 
f = function(x) (exp(x) - 1)/x

myDiff1 = function(x)
{
  if (x == 0) return(1)
  else return(f(x))
}

myDiff2 = function(x)
{
  y = exp(x)
  if (y == 1) return(1)
  else return( (exp(x) - 1)/log(exp(x)) )
}

dev.set(2)
vectPuiss = seq(-10,-18, by = -0.1)
X = 10^vectPuiss
U = sapply(X, myDiff1)
V = sapply(X, myDiff2)
N = length(vectPuiss)
plot(vectPuiss, U, type = 'l')
lines(vectPuiss, V, col = 'red')
abline(v = log(.Machine$double.eps)/log(10), col = 'blue', lwd = 2)

dev.set(3)
plot(vectPuiss, log(log(exp(X))), type = 'l', col = 'red')
lines(vectPuiss, log(X))