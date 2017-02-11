#' ---
#' title: "MACS205 : TD2"
#' output: html_document
#' author: Régis Gourdel
#' ---
#' 
#' **Exercice 1**
#' 
#' *1/*
#' On montre par récurrence avec le théorème de Rolle que pour tout $k \in \mathbf{N}_p$, $g^{(k)}$ s'annule $p + 1 - k$ fois.
#' En effet, entre chaque paire de zéros de $g^{(k)}$ on trouve un zéro de $g^{(k + 1)}$ d'après le théorème de Rolle.
#' 
#' Le résultat pour l'ordre $p$ est alors celui recherché.
#' 
#' 
#' **Exercice 3**
#' 
#' *1/*
#' Le résultat est immédiat en considérant la décomposition dans le base de Newton puisque $\omega_k$ est le seul polynôme dans la décomposition de degré au moins $k$.
#' 
#' *2/*
#' Un calcul immédiat montre $\forall 1 \leq i \leq k - 1, \tilde{p}_k(x_i) = f(x_i) = p_k(x_i)$.
#' On observe aussi que ce résultat est encore vrai pour $i = 0$ et $i = k$.
#' Donc $\tilde{p}_k$ et $p_k$ coïncident en $k + 1$ points.
#' Puisqu'ils sont de degré $k$ on en déduit $\tilde{p}_k = p_k$.
#' 
#' *3/*
#' Le résultat voulu s'obtient en prenant les coefficients dominant dans le résultat précédent.
#' 
#' *4/*
dividif = function(x,y){
  n = length(x) - 1
  d = y
  for (j in 2:(n + 1)){
    d[j : (n + 1)]
  }
}
#' 
#' 
#' 