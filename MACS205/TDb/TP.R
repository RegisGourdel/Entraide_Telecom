MCU = function(n, phi) {
  return(sum( lapply(runif(n), phi) )/n)
}

MCN = function(n, phi, m, var) {
  return(sum( lapply(rnorm(n,m,var), phi) )/n)
}

VC = function(n, phi, f, m, var) {
  X = rnorm(n,m,var)
  phiVec = lapply(X,phi)
  aux = function(x) { x - mean(x) }
  Z = lapply(fVec, aux)
  Psi = aux(phiVec)
}
