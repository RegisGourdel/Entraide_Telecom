#' ---
#' title: "MACS205 : TD3"
#' output: html_document
#' author: Régis Gourdel
#' ---
#' 
#' **Exercice 1**
#' 
#' *1/*
#' On a :
#' $$\hat{I}_{a,b,M} = h \left[ \frac{f(a) + f(b)}{2} + \sum_{m = 1}^{M - 1} f(a + mh) \right]\ .$$
#' Donc, on a en particulier $w_0 = w_M = \frac{1}{2}$ et $\forall m \in [\![ 1, M - 1 ]\!], w_m = 1$.
#'
#' *2/*
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

#' Test :
trapezeInt(cos,0,pi/2,10)
trapezeInt(cos,0,pi/2,100)

#' *3/*
#' Montrons que la méthode est d'ordre 1. Par linéarité de la formule il suffit de le montrer pour $1$ et $X$ sur tout intervalle.
#' Soient $a$ et $b$ réels, $a < b$.
#' On a
#' \begin{align*}
#' & \hat{I}_{a,b,1}(1) = h = b - a = \int_a^b 1 \\
#' & \hat{I}_{a,b,1}(X) = h \cdot \frac{a + b}{2} = \frac{(b - a)(b + a)}{2} = \frac{b^2}{2} - \frac{a^2}{2} = \int_a^b X
#' \end{align*}
#' 
#' Mais on a $\hat{I}_{0,1,1}(X^2) = 1 \cdot \frac{0 + 1}{2} = \frac{1}{2} \neq \frac{1}{3} = \int_a^b X^2$, ce qui prouve que la méthode est d'ordre $1$.
#' 

#' **Exercice 2**
#' 
refineTrapeze=function(FUN,a,b,M,q){
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

#' Test
p4 = function(x){x^4}
M = 5
myfun = p4
Qh = trapezeInt(myfun, 0, 1, M)
refineQh = refineTrapeze(myfun,0,1,M, Qh)
Qh2 = trapezeInt(myfun, 0, 1, 2*M)
err = Qh2 - refineQh
err


#' **Exercice 3**
#' 
#' *1/*
#' On a, par changement de variable affine $x = \frac{a + b}{2} + \frac{b - a}{2}t$ :
#' $$\lambda_i = \frac{1}{b - a} \int_a^b l_i(x)\ \mathrm{d}x = \frac{1}{2} \int_{-1}^1 l_i \left( \frac{a + b}{2} + \frac{b - a}{2}t \right) \mathrm{d}t$$
#' Notons $(t_i)_{0 \leq i \leq n}$ les $n + 1$ nœuds équidistants sur $[-1,1]$, de sorte que $x_i = \frac{a + b}{2} + \frac{b - a}{2}t$ et $x_i - x_j = \frac{b - a}{2}(t_i - t_j)$.
#' Alors on a
#' $$l_i \left( \frac{a + b}{2} + \frac{b - a}{2}t \right) = \prod_{j \neq i} \frac{b - a}{2} \cdot \frac{t - t_j}{x_i - x_j} = \prod_{j \neq i} \frac{t - t_j}{t_i - t_j}$$
#' d'où $\lambda_i = \frac{1}{2} \int_{-1}^1 \tilde{l}_i(u)\ \mathrm{d}u$.
#' 
#' *2/*
#' Dans le cas $n = 3$ on a $t_0 =-1, t_1 = 0, t_2 = 1$.
#' Alors $\lambda_0 = \frac{1}{2} \int_{-1}^1 \frac{x(x - 1)}{2} = \frac{1}{4} \left[ \frac{x^3}{3} - \frac{x^2}{2} \right]_{-1}^1 = \frac{1}{6}$.
#' Par symétrie il vient $\lambda_2 = \frac{1}{6}$ et on a
#' $$\sum \lambda_i = \frac{1}{2} \int_{-1}^1 \sum_i \tilde{l}_i(u)\ \mathrm{d}u = \frac{1}{2} \int_{-1}^1 L_3 \mathbf{1} = 1$$
#' donc $\lambda_1 = \frac{4}{6}$.
#' 
#' 
#' **Exercice 4**
#' 
#' *1/*
#' On pose $\forall i \in [\![ 0,M ]\!], x_i = a + ih$.
#' On a :
#' \begin{align*}
#' \hat{I}_{S,h} & = \sum_{i = 0}^{M - 1} \left[ \frac{1}{6} f(x_i) + \frac{2}{3} f \left( \frac{x_i + x_{i + 1}}{2} \right) + \frac{1}{6} f(x_{i + 1}) \right] \cdot h \\
#'               & = h \cdot \left( \frac{f(x_0) + f(x_M)}{6}  + \frac{1}{3} \sum_{i = 1}^{M - 1} f(x_i) + \frac{2}{3} \sum_{i = 0}^{M - 1} f \left( \frac{x_i + x_{i + 1}}{2} \right) \right) \\
#'               & = h \left( \frac{1}{6} f(a) + \frac{2}{3} f(a + h/2) + \sum_{i = 1}^{M - 1} \left[ \frac{1}{3} f(x_i) + \frac{2}{3} f(x_i + h/2) \right] + \frac{1}{6} f(b) \right)
#' \end{align*}
#' 
#' Or on a également :
#' $$\hat{I}_{T,h} = h \left( \frac{1}{2} f(a) + \sum_{i = 1}^{M - 1} f(x_i) + \frac{1}{2} f(b) \right)$$
#' et
#' $$\hat{I}_{T,h/2} = h \left( \frac{1}{4} f(a) + \frac{1}{2} f(a + h/2) + \sum_{i = 1}^{M - 1} \left[ \frac{1}{2} f(x_i) + \frac{1}{2} f(x_i + h/2) \right] + \frac{1}{4} f(b) \right)$$
#' d'où
#' $$\hat{I}_{S,h} = \frac{4}{3} \hat{I}_{T,h/2} - \frac{1}{3} \hat{I}_{T,h}\ .$$
#' 
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

#' *2/*
#' Test
d = 3
a = 1; b=5;

MyPolynome = function(x){x^d}
trueInt = (b^(d+1) - a^(d+1))/(d+1);

MM = 5:20
Resultats = rep(0,length(MM));
for (i in  1:length(MM)){
    Resultats[i]  = simpsonInt(MyPolynome,a,b,MM[i]);
    ## pi/2, 2*pi+pi/2, MM(i));
}

plot(MM,trueInt-Resultats);
title("error as a function of M")

if(max(abs(trueInt-Resultats)) > trueInt * .Machine$double.eps*10 ){
    logerr = log(abs(trueInt-Resultats));
    plot(log(MM),logerr);
    title("log error as  a function of log(M)")
    neval = length(logerr)
    pente = (logerr[neval] - logerr[1])/(log(MM[neval]) - log(MM[1]))
    pente
}

#' Plaçons nous sur un intervalle $[a,b]$ de taille $h = \frac{1}{M}$.
#' Notons $\epsilon$ l'epsilon machine.
#' On cherche
#' $$\min \left\{ d \in \mathbf{N} \mid \left\vert \frac{I(X^d) - \hat{I}(X^d)}{I(X^d)} \right\vert \geq 10 \cdot \epsilon \right\}$$
#' 
#' On emploie une méthode de Newton-Cotes de rang 2 donc d'ordre $N = 3$.
#' Donc l'erreur d'intégration est nulle pour tout $d \in [\![ 0,3 ]\!]$.
#' Donc la condition a une chance d'être vérifié à partir de $d = 4$ où l'erreur ainsi que l'erreur relative seront forcément non nulle sur certains intervalles.
#' 
#' 

#' **Exercice 5**
#' 
#' *1/*
#' Avec l'hypothèse $f^{(4)}(\xi_M) \simeq f^{(4)}(\xi) = C^{te}$ on a :
#' $$E_{2M}/E_M \simeq \frac{K f^{(4)}(\xi) (h/2)^4}{K f^{(4)}(\xi) h^4} = \frac{1}{16}\ .$$
#' 
#' *2/*
#' On a
#' $$I_{2M} - I_M = E_{2M} - E_M \simeq \frac{1}{16} E_M - E_M = -\frac{15}{16} E_M = - 15 \cdot E_{2M}$$
#' d'où
#' $$E_{2M} \simeq \frac{I_M - I_{2M}}{15}\ .$$

evalErrSimpson = function(FUN,a,b,M){
  ## Computes an approximation E of the error 
  ## for the composite Simpson rule of step h=(b-a)/(2M). 
  ##This requires computing I_M and I_{2M}. 
  ##The value  q = I_{2M} is also returned. 
  qth = trapezeInt(FUN,a,b,M)   ## M +1 evaluations
  qth2 = refineTrapeze ( FUN,a,b,M,qth )  ## M evaluations
  qth4 = refineTrapeze ( FUN,a,b,2*M,qth2 )   ## 2M evaluations
  simps_h  = 4/3*(qth2 - 1/4* qth)
  simps_h2 = 4/3*(qth4 - 1/4* qth2)
  q = simps_h2
  E = (simps_h - simps_h2)/15
  return(c(E,q))
}

#' *3/*
#'Test de l'évaluation de l'erreur :
d = 13.4 ; M=9
a = 0.5; b=3;
pol = function(x){x^d}
trueInt = (b^(d+1) - a^(d+1))/(d+1)
vv= evalErrSimpson(pol,a,b, M)
q = vv[2] ; E = vv[1]
estErr = E
trueErr = q - trueInt
relativeErrorOnError = (trueErr - estErr)/trueErr
relativeErrorOnError


#' **Exercice 6**
#'
#' *1/*
#' 
#' On a $p_{m,k}(x) = L_{k - m}^{x_=mM,\ldots,x_k} A(x)$.
#' Or on remarque par des calculs immédiat que $\frac{(x - x_m) p_{m + 1,k}(x) - (x - x_k) p_{m,k - 1}(x)}{x_k - x_m}$ coïncide avec $p_{m,k}(x)$ sur tous les $x_m,\ldots,x_k$, soit en $k + 1$ points disctincts.
#' Puisqu'ils sont de degré $k$ on en déduit l'égalité voulue.
#' 
#' *2/*
#' D'après le résultat précédent il vient :
#' \begin{align*}
#' q_{m,k} & = \frac{-x_m q_{m + 1,k} + x_k q_{m,k - 1}}{x_k - x_m} \\
#'         & = \frac{q_{m + 1,k} - \frac{x_k}{x_m} q_{m,k - 1}}{1 - \frac{x_k}{x_m}} \\
#'         & = \frac{q_{m + 1,k} - \delta^{k - m} q_{m,k - 1}}{1 - \delta^{k - m}}
#' \end{align*}
#' 
#' *3/*
#' En posant $i = k$ et $j = i - m$ on déduit directement du résultat précédent :
#' $$A_{i,j} = q_{i - j,i} = \frac{q_{i - j + 1,i} - \delta^{j} q_{i - j,i - 1}}{1 - \delta^j} = \frac{A_{i,j - 1} - \delta^{j} A_{i - 1,j - 1}}{1 - \delta^j}\ .$$
#' 
#' 
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

#' Test :
myfun = function(x){sin(x)+ (cos(x)-1)*sin(2*x)}
n =10 ; t =pi/4 ;  delta = 1/4
A = richardson(myfun,n,t,delta)
lx = log(t) +  log(delta)*(0:n)
x = exp(lx);
y = sapply(x,myfun);

plot(0:n, y,col='blue', type = "l")
lines( 0:n,A,col='red');
legend('topright', legend=c('erreur naive', 'erreur Richardson'),
       col=c('blue', 'red'), lwd=2)

#'
dev.new()
LerrRich = log(abs(A - myfun(0))) ;
LerrNaive = log(abs(y - myfun(0)));
plot(-lx, LerrNaive,col='blue',type='l')
lines(-lx, LerrRich,col='red')
grid()
legend('topright', legend=c('log-erreur naive','log-erreur Richardson'),
       col=c('blue', 'red'), lwd=2)


#' **Exercice 7**

romberg = function(FUN,n,a,b,M){## methode de Romberg avec n etapes
  ## appliquee sur la fonction FUN sur l'intervalle (a,b), avec un
  ## pas initial h = (b-a)/M
  h= (b-a)/M 
  A = rep(0, n+1)
  A[1] = trapezeInt(FUN,a,b,M);
  Mc = M
  ## initialisation des differences divisees
  for( i in 2:(n+1)){
    A[i] = refineTrapeze( FUN,a,b, Mc, q= A[i-1])
    Mc = 2*Mc 
  }
  delta = 1/4;
  for (j in 2:(n+1)){
    ## A[j : (n+1) ] = ## completer sur le modele de richardson
  }
  return(A)
}


#' test
# d=6
# myfun = function(x){x^d}
# n =10; M = 5; a= 0 ; b = 1.5
# A = romberg(myfun, n, a, b , M)
# B = rep(0,n+1)
# B[1] = trapezeInt(myfun,a,b,M) 
# Mc = M
# for( i in  2:(n+1)){
#     B[i] = refineTrapeze( myfun,a,b, Mc, B[i-1])
#     Mc = 2* Mc 
# }
# 
# plot(0:n, B,col='blue', type='l')
# lines( 0:n, A,col = 'red');
# 
# LerrRich = log(abs(A - (b^(d+1) - a^(d+1))/(d+1))) ;  
# LerrNaive = log(abs(B - (b^(d+1) - a^(d+1))/(d+1)));
# plot((0:n), LerrNaive/log(4),col='blue',type='l')
# lines( (0:n), LerrRich/log(4), col='red')
# grid()
# legend('topright', legend=c('log-erreur naive', 'log-erreur Romberg'),
#        col=c('blue', 'red'), lwd=2)