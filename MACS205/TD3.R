#' ---
#' title: "MACS205 : TD3"
#' output: html_document
#' author: Régis Gourdel
#' ---
#' 
#' **Exercice 1**
#' 
#' *1/*

trapezeInt =function(FUN,a,b,M){
    ##' TRAPEZOIDAL INTEGRATION RULE (COMPOSITE)
    ##' @param FUN : the function to be integrated
    ##' @param a, b : interval end points 
    ##' @param M : number of intervals (each of size (b-a)/M)
    ##' @return: the value of the composite trapezoidal quadrature. 
    x = seq(a,b, length.out= M+1)
    y = sapply(x, FUN)
    h = (b-a)/M 
    ## q = ## Complete the code here. Use vector operations, not a loop !
    return(q)
}

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
    ## x =  ## complete here : 
        ##  x : a vector of size M :
        ##     the additional abscissas where 'fun' must be evaluated.
        y = sapply(x, FUN)
    ## Q = ##  complete here : a function of y, h and q. 
    return(Q)
}

#' Test
# p4 = function(x){x^4}
# M = 5
# myfun = p4 
# Qh = trapezeInt(myfun, 0, 1, M)
# refineQh = refineTrapeze(myfun,0,1,M, Qh)
# Qh2 = trapezeInt(myfun, 0, 1, 2*M)
# err = Qh2 - refineQh
# err


#' **Exercice 4**
simpsonInt = function(FUN,a,b,M){
    ##' Simpson integration via trapeze rule
    ##' uses the fact that 
    ##' simpson(h) = 4/3(trapeze(h/2) - 1/4 trapeze(h))
    h = (b-a)/M;
    ## qtrapeze = ## Complete  
    ## qrefined = ## Complete 
    ## q =  ## Complete 
    return(q)
}

#' Test
# d = 4
# a = 1; b=5;
# 
# MyPolynome = function(x){x^d}
# trueInt = (b^(d+1) - a^(d+1))/(d+1);
# 
# MM = 5:20
# Resultats = rep(0,length(MM));
# for (i in  1:length(MM)){
#     Resultats[i]  = simpsonInt(MyPolynome,a,b,MM[i]);
#     ##pi/2,2*pi+pi/2,MM(i));
# }
# 
# plot(MM,trueInt-Resultats);
# title("error as a function of M")
# 
# if(max(abs(trueInt-Resultats)) > trueInt * .Machine$double.eps*10 ){
#     logerr = log(abs(trueInt-Resultats));
#     plot(log(MM),logerr);
#     title("log error as  a function of log(M)")
#     neval = length(logerr)
#     pente = (logerr[neval] - logerr[1])/(log(MM[neval]) - log(MM[1]))
#     pente
# }


#' **Exercice 5**
evalErrSimpson = function(FUN,a,b,M){
    ## Computes an approximation E of the error 
    ## for the composite Simpson rule of step h=(b-a)/(2M). 
    ##This requires computing I_M and I_{2M}. 
    ##The value  q = I_{2M} is also returned. 
    qth = trapezeInt(FUN,a,b,M)   ## M +1 evaluations
    qth2 = refineTrapeze ( FUN,a,b,M,qth )  ## M evaluations
    qth4 = refineTrapeze ( FUN,a,b,2*M,qth2 )   ## 2M evaluations
    simps_h =   4/3*(qth2 - 1/4* qth ) 
    ## simps_h2 =  ## Completer 
        q = simps_h2  
    ## E = ## Completer 
        return(c(E,q))
}

#' Test
# d = 13.4 ; M=9
# a = 0.5; b=3;
# pol = function(x){x^d}
# trueInt = (b^(d+1) - a^(d+1))/(d+1)
# vv= evalErrSimpson(pol,a,b, M)
# q = vv[2] ; E = vv[1]
# estErr = E
# trueErr = q - trueInt 
# relativeErrorOnError = (trueErr - estErr)/trueErr
# relativeErrorOnError

#' **Exercice 6**
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
        ## A[j : (n+1) ] =  ## Completer le code 
    }
    return(A)
}

#' test
# myfun = function(x){sin(x)+ (cos(x)-1)*sin(2*x)}
# n =10 ; t =pi/4 ;  delta = 1/4
# A = richardson(myfun,n,t,delta)
# lx = log(t) +  log(delta)*(0:n)
# x = exp(lx);
# y = sapply(x,myfun);
# 
# plot(0:n, y,col='blue', type = "l")
# lines( 0:n,A,col='red');
# legend('topright', legend=c('erreur naive', 'erreur Richardson'),
#        col=c('blue', 'red'), lwd=2)
# 
# dev.new()
# LerrRich = log(abs(A - myfun(0))) ;  
# LerrNaive = log(abs(y - myfun(0)));
# plot(-lx, LerrNaive,col='blue',type='l')
# lines(-lx, LerrRich,col='red')
# grid()
# legend('topright', legend=c('log-erreur naive','log-erreur Richardson'),
#        col=c('blue', 'red'), lwd=2)


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