# Fonctions prédéfinies

evalID = function (x){
  #### usage: id = evalX (x)
  #### returns the index ID such that r = ID * 2^(-BETA) is the 
  #### best approximation of x in the set { k * 2^(-BETA), k=1 ... 1.5*2^BETA}
     if (any(x < min(ABSCISSES) - 2^(-BETA - 1)  )|| any(x > max(ABSCISSES)  + 2^(-BETA - 1)) ){
         warning(' x should be within  the range of ABSCISSES ');
     }
     id = round(2^BETA * x)
     id = sapply(id, function(x){max(x,1)})
     ## max(id,1)
     id = sapply(id, function(x){min(x, length(ABSCISSES))})
     ##min(id, length(ABSCISSES))
     
     return(id)
}


evalX = function(x){
  #### usage: absc = evalX (x)
  #### returns the closest number to x which writes r = id * 2^beta , 1 <= id <= 1.5*2^BETA
  #### where id is an integer. 
  id = evalID(x) ; 
  absc  = id * 2^(-BETA) ;
return(absc)
}


evalBoiteNoire = function(x){
  ## USAGE : f = evalBoiteNoire (x)
  ## x:   a vector [1, n] : the abscissas where one wants to evaluate
  ##    the function BoiteNoire
  ## VALUE : 
  ## f:   a vector of same size than x : the values of BoiteNoire at 
  ##      at the closest points in the 'ABSCISSES' vector. 
  id = evalID(x)
  f = VALEURS[id]
    return(f)
}