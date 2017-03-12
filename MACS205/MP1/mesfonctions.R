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

piecewiseInterpol = function(n, nInt, a, b, neval, nodes = "equi", FUN, Plot){
    ## @param n : the degree of the interpolating polynomial on each
    ## subinterval
    ## @param nInt :  the number of sub-intervals
    ## @param a, b : endpoints of the interval
    ## @param neval : the number of points on the interpolating grid (on
    ## each subinterval)
    ## @param nodes : string, either "equi" (default) for equidistant
    ## Lagrange interpolation (on each subinterval) or "cheby" for
    ## chebyshev nodes.
    ## @param FUN the function to be interpolated
    ## @param Plot : logical. Should the result be plotted ?
    ## @return : a matrix with 2 rows and neval * nInt -neval + 1:
    ## values of the interpolated funtion on a regular grid (first row)
    ## and the corresponding abscissas (second row).
    
    intEndPoints = seq(a,b,length.out = nInt+1)
    f = c()
    z = c()
    if (nodes == "equi"){
        x = seq(-1,1,length.out = (n + 1))
    }
    else if (nodes == "cheby"){
        rtche = function(t) { cos((2*t - 1)*pi/2/(n + 1)) }
        x = rtche( seq(1:(n + 1)) )
    }
    for (m in 1:nInt){
        A = intEndPoints[m]; B = intEndPoints[m+1]
        xm = x * (B - A)/2 + (A + B)/2
        zm = seq(A, B, length.out = neval)
        fm = interpolDividif(xm, FUN(xm), zm)
        if( m >= 2){
            ## remove first element of zm, fm to avoid
            ## duplicate values of the  interpolating vector
            zm = zm[2:neval]
            fm = fm[2:neval]
        }
        z = c(z,zm)
        f = c(f,fm)
    }
    if (Plot == 1){
        if (nodes == "equi") {methodName = " equidistant "}
        else  {methodName = " Chebyshev "}
        plot(z, sapply(z,FUN),type="l", lwd = 2)
        title(main = paste("Piecewise  Lagrange  interpolation with ",
                           toString(n+1), methodName, " nodes  on ",
                           toString(nInt), " Intervals", sep=""))
        lines(z,f, col='red', lwd=2, lty = 2)
        legend('topright', legend = c('function','interpolation'),
               lwd=c(1,2), col=c('black','red'))
    }
    return(rbind(f,z))
}

trapezeInt = function(FUN,a,b,M){
    ##' TRAPEZOIDAL INTEGRATION RULE (COMPOSITE)
    ##' @param FUN : the function to be integrated
    ##' @param a, b : interval end points 
    ##' @param M : number of intervals (each of size (b-a)/M)
    ##' @return: the value of the composite trapezoidal quadrature. 
    x = seq(a,b, length.out= M+1)
    y = sapply(x, FUN)
    h = (b-a)/M
    q = h * (y[1]/2 + sum(y[2:M]) + y[M + 1]/2)
    return(q)
}

refineTrapeze = function(FUN,a,b,M,q){
    ##' refinement of the subdivision step: incremental method
    ##' @param FUN : the function to be integrated
    ##' @param a, b : interval end points 
    ##' @param M : initial number of intervals (each of size (b-a)/M)
    ##'  having been used to compute q
    ##' @param  q : the value of the trapezoidal  quadrature method
    ##'  of stepsize (b-a)/M
    ##' @return : the value of the quadrature for a stepsize h' = h/2
    h = (b-a)/M
    x = seq(a,b, length.out = M+1)[1:M] + h/2
    ##  x : a vector of size M :
    ##     the additional abscissas where 'fun' must be evaluated.
    y = sapply(x, FUN)
    Q = q/2 + sum(y)*(h/2)
    return(Q)
}

simpsonInt = function(FUN,a,b,M){
    ##' Simpson integration via trapeze rule
    ##' uses the fact that 
    ##' simpson(h) = 4/3(trapeze(h/2) - 1/4 trapeze(h))
    h = (b-a)/M;
    qtrapeze = trapezeInt(FUN,a,b,M)
    qrefined = refineTrapeze(FUN,a,b,M,qtrapeze)
    q = 4/3*qrefined - 1/4*qtrapeze
    return(q)
}

evalErrSimpson = function(FUN,a,b,M){
  ## Computes an approximation E of the error 
  ## for the composite Simpson rule of step h=(b-a)/(2M). 
  ## This requires computing I_M and I_{2M}. 
  ## The value  q = I_{2M} is also returned. 
  qth = trapezeInt(FUN,a,b,M)   ## M +1 evaluations
  qth2 = refineTrapeze ( FUN,a,b,M,qth )  ## M evaluations
  qth4 = refineTrapeze ( FUN,a,b,2*M,qth2 )   ## 2M evaluations
  simps_h  = 4/3*(qth2 - 1/4* qth)
  simps_h2 = 4/3*(qth4 - 1/4* qth2)
  q = simps_h2
  E = (simps_h - simps_h2)/15
  return(c(E,q))
}

richardson = function(FUN,n,t,delta){
    ## Calcule le tableau des differences  divisees en 0 du 
    ## polynome d'interpolation en t,delta t, ... delta^n t
    ## renvoie un vecteur de taille n+1:
    ## le vecteur des A_{k,k}, k= 0 .. n 
    ## (pas la matrice).   
    ## La meilleure approximation est le dernier element A[n+1].
    ##
    lx = log(t)  +  log(delta) *(0:n)
    x = exp(lx) 
    A = sapply(x,FUN) 
    for( j in 2:(n+1)){
        A[j : (n+1)] = (A[j : (n+1)] - delta^j * A[(j-1) : n]) / (1 - delta^j)
    }
    return(A)
}

romberg = function(FUN,n,a,b,M){
    ## methode de Romberg avec n etapes
    ## appliquee sur la fonction FUN sur l'intervalle (a,b), avec un
    ## pas initial h = (b-a)/M
    A = rep(0, n+1)
    A[1] = trapezeInt(FUN, a, b, M);
    Mc = M
    ## initialisation des différences divisees
    for( i in 2:(n+1) ){
        A[i] = refineTrapeze(FUN, a, b, Mc, q = A[i-1])
        Mc = 2*Mc 
    }
    B = A
    delta = 1/4;
    for (j in 2:(n+1)){
        A[j : (n+1)] = (A[j : (n+1)] - delta^j * A[(j-1) : n]) / (1 - delta^j)
    }
    return(rbind(A, B))
}