library(fda)
library(coda)
library(truncdist)
library(mvtnorm)
library(dplyr)
library(dfosr)

#MGP generator, based on Durante,Daniele.2017."A Note On The Multiplicative Gamma Process". 
MGP = function(nu,a1 = 2, a2 = 3, k, P){
  deltas = c(rgamma(1,a1,1),rgamma(k-1, a2, 1))
  taus = exp(cumsum(log(deltas))) 
  Phi = matrix(rgamma(k*P,nu/2,nu/2),nrow = P, ncol = k,byrow = TRUE) #p by k
  Lambda = t(apply(Phi,1, function(phi){return(rmvnorm(1,rep(0,k), diag(1/(phi*taus))))}))
  return(list(Lambda,Phi,taus))

  }
get_Bmat = function(Y){
  T = ncol(Y); n = nrow(Y)
  # Knots
  num_int_knots = max(20,min(ceiling(T/4), 150))
  quantile_locations = seq(0,1, length=(num_int_knots+2))[-c(1,(num_int_knots+2))]
  knots = quantile(seq(0,1,length.out = T), quantile_locations)
  
  # Basis functions: B-splines
  b = create.bspline.basis(rangeval = c(0,1),breaks = knots,norder = 4)

  #Time Grid
  Lt_f <- lapply(seq(n), function(i) which(!is.na(Y[i,]))/T) 
  Bmat = list()
  
  #Full basis matrix:
  B = eval.basis(seq(T)/T, b) #T by p
  #A list of  partial basis matrix:
  Bmat = lapply(seq(n),function(j){eval.basis(Lt_f[[j]],b) })
  
  return(list(B = B, Bmat = Bmat))
}
MakeLong <- function(wide){
  D <- nrow(wide)
  
  long <- as.data.frame(wide) %>% 
    mutate(.index = 1:D) %>% 
    gather(.id, .observed, -.index) %>%
    mutate(.id = paste0("I", .id)) %>% 
    as_tibble()
  long
}
plot_F = function(wide,alpha = .2,col = FALSE){
  if(col){t(wide) %>% MakeLong() %>%
      ggplot(aes(x = .index, y = .observed, group = .id,color = .id)) +
      geom_line(alpha = alpha)}
  else{t(wide) %>% MakeLong() %>%
    ggplot(aes(x = .index, y = .observed, group = .id)) +
    geom_line(alpha = alpha)}
}
