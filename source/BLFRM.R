source(here::here("source", "BLFRM_Utilities.R"))
BLFRM = function(X,Y,nu = 5,k = 5,a_sigma = 0.5, b_sigma = 0.5,a_psi = 0.5, b_psi = 0.5,a1 = 2, a2 = 3,
                 S = 10000){
# inputs:
#   X--data.frame with external covariates
#   Y--functional data.frame 
#   nu,a_1,a_2--parameters of the multiplicative gamma process
#   a_sigma,b_sigma--parameters for sigma term prior
#   a_psi,b_psi--parameters for psi term prior
#   k--how many FPCs/functional prototypes to estimate

# returns a list with elements
#   post_Sigma: posterior for Sigma
#   post_Delta: posterior for Delta
#   post_Beta: posterior for Beta
#   post_Eta: posterior for Eta
#   post_Theta: posterior for Theta
#   post_Yhat: posterior for Y_hat
#   B: the full basis matrix

n = nrow(Y);T = ncol(Y);r = ncol(X)

#Non-na entries in Y
Ts = sapply(seq(nrow(Y)), function(i) sum(!is.na(Y[i,])))  

# Basis matrix:
Bm = get_Bmat(Y)$B
Bmat =get_Bmat(Y)$Bmat
p = ncol(Bm)

#Initialize Lambda, Phi, Tau with MGP
MGP_init = MGP(nu, a1 = a1, a2 = a2, k = k ,P =p) 
Lambda = MGP_init[[1]] #k by p
Phi = MGP_init[[2]];#k by p
Tau = MGP_init[[3]] #k by 1

post_Sigma = array(dim = c(p,p,S))
post_Delta = array(dim = c(n,k,S))
post_Beta = array(dim = c(k,r,S))
post_Eta = array(dim = c(n,k,S))
post_Theta = array(dim = c(p,n,S))
post_Yhat = array(dim = c(n,T,S))

#Initialize other parameters
psi2 = c(1/rgamma(1,a_psi, b_psi),rep(0,S-1)) #1 by S
post_Sigma[,,1] = diag(1/(rgamma(p,a_sigma,b_sigma))); #p by p
Zeta =  rmvnorm(n,sigma = post_Sigma[,,1]) #n by p
post_Delta[,,1] = rmvnorm(n, mean = rep(0,k),sigma = diag(k)) #n by k
post_Beta[,,1] =  rmvnorm(k,sigma = diag(rgamma(r,1/2,1/2)))#k by r 
post_Eta[,,1] = t(post_Beta[,,1]%*%t(X)) + post_Delta[,,1] #n by k
post_Theta[,,1] = sapply(seq(n), function(i) Lambda%*% post_Eta[i,,1] + Zeta[i,]) #p by n
post_Yhat[,,1] =  t(sapply(seq(n),function(i)rmvnorm(1,Bm%*%post_Theta[,i,1],psi2[1]*diag(T)))) #n by T

for(i in seq(2,S)){
  for(j in seq(p)){
    #1): Update Lambda
    #2): Update Phi
    D_inv = diag(Phi[j,]*Tau)
    sigma_j = post_Sigma[j,j,i - 1]
    Eta_T_Eta = crossprod(post_Eta[,,i-1])
    Q_lambda = D_inv + sigma_j^(-1)*Eta_T_Eta
    l_lambda =sigma_j^(-1) * t(post_Eta[,,i-1])%*%post_Theta[j,,i - 1]
    Lambda[j,] <- rmvnorm(1, crossprod(solve(Q_lambda), l_lambda), solve(Q_lambda));
    Phi[j,] = sapply(1:k,function(h) rgamma(1,(nu+1)/2, nu/2 + (Tau*Lambda[j,h]^2/2)))
  }
  #3): Update Tau
  Deltas = array(dim = k)
  tau_h_1 = 1
  for(h in 1:k){
    tau = c(rep(1,h - 1),rep(tau_h_1/Tau[h],k - h + 1))*Tau
    beta1 = 1 + 1/2*sum(sapply(h:k,function(l)tau[l]*(Phi[,l]%*%Lambda[,l]^2)))
    Deltas[h] = rgamma(1, a2 + (p/2)*(k - h +1),beta1)
    Tau[h] = tau_h_1*Deltas[h]
    tau_h_1 = Tau[h]
  }

  #2): Update sigma2
  Sigma = array(dim = p)
  for(j in 1:p){
    b_sig1 = b_sigma+sum(sapply(1:n,function(l)crossprod(post_Theta[j,l,i - 1]-Lambda[j,]%*%post_Eta[l,,i - 1])/2))
    a_sig1 = n/2 + a_sigma
    Sigma[j] = rgamma(1, a_sig1,b_sig1)
  }
  post_Sigma[,,i] =  diag(1/Sigma)
  
  #3): Update psi
  y = as.vector(Y)%>% na.exclude()
  y_hat = sapply(seq(n),function(j) Bmat[[j]] %*% post_Theta[,j,i - 1]) %>% unlist()
  psi2[i]  = 1/rgamma(1,a_psi + n*T/2,crossprod(y - y_hat)/2 + b_psi)
  #4)Update beta and w
  for(j in 1:k){
   # w_l = sapply(1:r,function(l) rgamma(1,1/2*(1+post_Beta[j,l,i - 1]^2)))
    w_l = rgamma(r,rep(1,r),1/2*(1+post_Beta[j,,i - 1]^2))
    Q_beta = crossprod(X) + solve(diag(1/w_l))
    l_beta = t(X)%*%post_Eta[,j,i - 1]
    post_Beta[j,,i] = rmvnorm(1,solve(Q_beta)%*%l_beta,solve(Q_beta))
  }
  
  #5) Update eta
  for(j in 1:n){
    #recurring term
    M = t(Lambda)%*%t(Bmat[[j]])%*%solve(psi2[i]*diag(Ts[j]) + Bmat[[j]]%*%post_Sigma[,,i]%*%t(Bmat[[j]]))
    A = M %*%Bmat[[j]]%*%Lambda + diag(k)
    B = post_Beta[,,i]%*%X[j,] + M %*% na.exclude(Y[j,])
    post_Eta[j,,i] = rmvnorm(1,solve(A)%*%B, solve(A))
  }
  
  #Update theta
  for(j in 1:n){
    Q_theta = 1/psi2[i]*crossprod(Bmat[[j]]) + solve(post_Sigma[,,i])
    l_theta = 1/psi2[i]*t(Bmat[[j]])%*% na.exclude(Y[j,])+  solve(post_Sigma[,,i])%*%Lambda%*%post_Eta[j,,i]
    post_Theta[,j,i] = rmvnorm(1,solve(Q_theta)%*%l_theta,solve(Q_theta))
  }
  
  #Update y_hat
  post_Yhat[,,i] =  t(sapply(seq(n),function(j)rmvnorm(1,Bm%*%post_Theta[,j,i],psi2[i]*diag(T)))) 
}
return(list(post_Sigma = post_Sigma,
            post_Delta = post_Delta,
            post_Beta = post_Beta,
            post_Eta = post_Eta,
            post_Theta = post_Theta,
            post_Yhat = post_Yhat,
            B = Bm))
}
