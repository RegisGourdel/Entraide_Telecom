# Fonctions pour la première question

dividif = function(x,y){
  ## Newton's Divided differences
  ## @param x: a vector containing the interpolation nodes 
  ## @param y: a vector of same size as x:
  ##           values of the interpolated function at the nodes
  ## @return : a vector of same size as x:
  ##          the divided differences
  ##          \eqn{f_[x_0, ... x_k]} of order 'length(x) -1'.
  n = length(x) - 1 ## degree of the Lagrange polynomial
  d  = y
  if (n > 0) {
    for (j in 2:(n+1)) {
      d[j : (n + 1)] = (d[j : (n + 1)] - d[(j - 1) : n]) / (x[j  : (n + 1)] - x[1 : (n - j + 2)])
    }
  }
  return(d)
}

hornerNewton = function(a,x,z){
  ## Horner's method: Evaluates  a polynom P at points z, given
  ## nodes x and the coefficients a of P in Newton's basis
  ## @param a : vector: the  coefficients of the polynomial in Newton's basis
  ## @param x : the interpolation nodes. 
  ## @param z : vector of points where the polynom needs to be evaluated. 
  ## @return  : a vector of same size as z: the value of the polynomial at points z.
  n = length(x)
  f = a[n] * rep(1,length(z))
  if (n >= 2){
    for(i in 1:n-1){
      f = f * (z - x[n - i]) + a[n - i]
    }
  }
  return(f)
}

interpolDividif = function(x,y,z){
  ## Efficient Lagrange interpolation using Horner's method with  
  ## Newton basis for evaluation
  ## @param x : vector containing the interpolation nodes 
  ## @param y : vector of same size as x: values of the interpolated
  ##            function at the nodes
  ## @param z : vector of points where the  interpolating polynomial
  ##            needs to be evaluated. 
  ## @return  : vector of same size as z: the value of the
  ##            interpolating polynomial at points z.
  d = dividif(x,y)
  return(hornerNewton(d,x,z))
}