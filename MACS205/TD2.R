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
#' *2/*
#' Pour tout $x \in \{ x_0,\ldots, x_n \}$ on a $f(x) - L_n f(x) = 0$ donc la formule (1) est vérifiée en prenant le $\xi_x$ donné par la question précédente.
#' 
#' *3/*
#' On a
#' $$L_{n + 1}f - L_n f = \sum_{i = 0}^{n + 1} f(x_i) l_i - \sum_{i = 0}^{n} f(x_i) l_i = f(x_{i + 1}) l_{i + 1} = c \omega_{n + 1}\ .$$
#' avec $c = \frac{f(x_{i + 1})}{ \prod_{i = 0}^{n} x_{n + 1} - x_i }$.
#' 
#' *4/*
#' D'après la question 2, la fonction $g$ s'annule en $n + 1$ points.
#' En lui appliquant le résultat du lemme 1 il vient donc $\exists \xi \in [ \min(x_0,x), \max(x_0,x)], g^{(n + 1)}(\xi) = 0$.
#' En remarquant que $g = f - L_n f - c \omega_{n + 1}$ on obtient donc le résultat voulu.
#' 
#' *5/*
#' On sait que $L_n f$ est un polynôme de degré $n$, d'où $(L_n f)^{(n + 1)}(\xi) = 0$.
#' De plus $\omega_{n + 1}$ est de degré $n + 1$, donc $\omega_{n + 1}^{(n + 1)}$ est constant égal à $(n + 1)!$.
#' Il vient donc $c = \frac{f^{(n + 1)}(\xi)}{(n + 1)!}$.
#' 
#' *6/*
#' On a
#' $$E_n f = L_{n + 1} f - L_n f $$
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
  ##  Newton's Divided differences
  ## @param x: a vector containing the interpolation nodes 
  ## @param y: a vector of same size as x:
  ##           values of the interpolated function at the nodes
  ## @return : a vector of same size as x:
  ##          the divided differences
  ##          \eqn{f_[x_0, ... x_k]} of order 'length(x) -1'. 
  
  n = length(x) - 1 ## degree of the Lagrange polynomial
  d  = y
  for (j in 2:(n+1)) {
    d[j : (n + 1)] = (d[j : (n + 1)] - d[(j - 1) : n]) / (x[j  : (n + 1)] - x[1 : (n - j + 2)])
  }
  return(d)
}
#' 
#' 
#' **Exercice 4**
#' 
#' *1/*
eval = function(a,x){
  n = length(a) - 1
  if (n == 0) return(a[1])
  else{
    a[n] = a[n] + x*a[n + 1]
    return(eval(a[1:n], x))
  }
}

#' *2/*
hornerNewton = function(a,x,z){
  ## Horner's method: Evaluates  a polynom P at points z, given
  ## nodes x and the coefficients a of P in Newton's basis
  ##
  ## @param a : vector: the  coefficients of the polynomial in
  ##           Newton's basis
  ## @param x : the interpolation nodes. 
  ## @param z : vector of points where the polynom needs to be
  ##            evaluated. 
  ## @return  : a vector of same size as z: the value of the
  ##            polynomial at points z.
  ## 
  n = length(x);
  f  = a[n] * (z-x[n-1]) + a[n-1]
  for(i in 2:n-1){
    f = f * (z - x[n - i]) + a[n - i]
  }
  return(f)
}

#' **Exercice 5**
interpolDividif=function(x,y,z){
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

myfun = function(x){ cos(5*x/(2 + x)) / (1 + (1 + x)^2) }
n = 20
x = seq(-1,1,length.out = n)
z = seq(-1.4,1.4,length.out = 10*n)
val = interpolDividif(x,myfun(x),z)

plot(z, myfun(z), type = 'l', col="red")
lines(z, val, type = 'l', col="blue", lty = 2)
abline(v = -1)
abline(v = 1)

#' **Exercice 6**
funRunge = function(x){ 1/(1 + (5*x)^2) }

#' *1/*
n = 8
m = 100
x = seq(-1,1,length.out = n)
x
z = seq(-1,1,length.out = m)
val = interpolDividif(x, funRunge(x), z)

Y = funRunge(z)
plot(z, Y, type = 'l', col="red")
lines(z, val, type = 'l', col="blue", lty = 2)
lines(x, funRunge(x), type = 'p', col="blue")

#' *2/*
xtche = sapply( pi/n*seq(0, (n-1)) + pi/(2*n), cos )

val = interpolDividif(xtche, funRunge(xtche), z)
plot(z, Y, type = 'l', col="red")
lines(z, val, type = 'l', col="blue", lty = 2)
lines(xtche, funRunge(xtche), type = 'p', col="blue")

#' **Exercice 7**
#' 
#' *1/*
interpolLagrange = function(n, a, b, neval, nodes = 'equi', FUN, Plot){
    ## Generic Lagrange interpolation, with equidistant or Chebyshev nodes. 
    ## @param n : the degree of the interpolating polynomial on each
    ## subinterval
    ## @param a : left end-point of the interval
    ## @param b : right end-point of the interval
    ## @param neval :number of evaluation points (a regular grid will be
    ## used on [a,b]
    ## @param nodes :string, either "equi" (default) for equidistant
    ## Lagrange interpolation (on each subinterval) or "cheby" for
    ## using Chebyshev nodes.
    ## @param FUN: the function to be interpolated 
    ## @param Plot : logical. Setting 'Plot' to TRUE produces a plot
    ## showing the graph of
    ## the true functions and its interpolation.  
    ## @return : vector of size 'neval': the values of the Lagrange
    ## polynomial on an equi-distant grid.
    
    if (nodes == "equi"){
        x = seq(a,b,length.out = n)
    }
    else if (nodes == "cheby"){
        x = ( sapply( pi/n*seq(0,(n-1)) + pi/(2*n), cos ) + (a + b)/2 ) * (b - a)
    }
    else{stop("the nodes must be either 'equi' or 'cheby'") }
    
    z = seq(a,b,length.out = neval)
    f = interpolDividif(x, FUN(x), z)
    
    if( Plot ){
        if (nodes == "equi"){ methodName = " equidistant "}
        else { methodName = " Chebyshev "}
        
        y = sapply(z,FUN)
        plot(z, y, type="l", ylim=range(c(y,f)) )
        title(main = paste("Lagrange interpolation with ",
                           toString(n+1), methodName,
                           " nodes", sep=""))
        lines(z,f, col = 'blue') 
        legend('topright',
               legend=c('function','interpolation'),
               col = c('black','red'),
               lwd=1)
    }
    return(f)
}

#' *2/*
#' Test :
interpolLagrange(6,-1,1,100,'cheby',myfun,TRUE)
interpolLagrange(8,-1,1,100,'equi',funRunge,TRUE)
