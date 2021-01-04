##################################################################
############   Simulate Data for Heterogenous case    ############

#### B-spline
P = 4
day = 30
t = seq(0,1,length.out=day)
basis = create.bspline.basis(c(0,1),nbasis = P)
s = getbasismatrix(t,basis)
#### orthogonalization
b = qr.Q(qr(s))

######   mean functions   ######

mu = list()
mu[[1]] = exp(t)*cos(t)/2
mu[[2]] = cos(5*pi/2*t)
cluster.num = length(mu)

lambda_true <- list()
for(i in 1:cluster.num){
  lambda_true[[i]] = solve(t(b)%*%b,t(b)%*%mu[[i]])
}

#### 
Theta_true = list()
Theta_true[[1]] = matrix(c(0,1,0,0,0,0,1,0),,2)
Theta_true[[2]] = matrix(c(0,0,0,1,1,0,0,0),,2)
Theta_true[[1]] = qr.orth(Theta_true[[1]])
Theta_true[[2]] = qr.orth(Theta_true[[2]])


FPC = list()
FPC[[1]] = b %*% Theta_true[[1]]
FPC[[2]] = b %*% Theta_true[[2]]

####
sigma_gamma_true = matrix(c(4,2,1,0.5),2,2)
sigma_true = 0.4

phi_true = 1


###########      location     ###########

station = read.csv('../../../data/city_location_China.csv',sep=',',fileEncoding  = 'GBK')  
station = station[ 107<station$lon & station$lon<125 & 28<station$lat & station$lat<43, ]
location = station[,3:4]

curve.num = dim(location)[1]


##########    Simulate markov random field
source('../homo/SimulateMRF.R')

###############   Data Generation


data_list = list()

Nsimu = 100
for( ite in 1:Nsimu){
  cluster_true = mrf[,ite]
  dd = matrix(NA,curve.num,day)
  
  #######  random effect  ######
  
  Gamma = list()
  Gamma[[1]] = list()
  Gamma[[1]][[1]] =  sigma_gamma_true[1,1] * Exp.cov(location[cluster_true==1,],theta=phi_true)
  Gamma[[1]][[2]] =  sigma_gamma_true[1,2] * Exp.cov(location[cluster_true==1,],theta=phi_true)
  Gamma[[2]] = list()
  Gamma[[2]][[1]] =  sigma_gamma_true[2,1] * Exp.cov(location[cluster_true==2,],theta=phi_true)
  Gamma[[2]][[2]] =  sigma_gamma_true[2,2] * Exp.cov(location[cluster_true==2,],theta=phi_true)
  
  gamma =list()
  gamma[[1]] = matrix(NA,dim(Gamma[[1]][[1]])[1],2)
  gamma[[2]] = matrix(NA,dim(Gamma[[2]][[1]])[1],2)
  
  set.seed(312)
  gamma[[1]][,1] = mvrnorm(n=1,rep(0,dim(Gamma[[1]][[1]])[1]),Gamma[[1]][[1]])
  set.seed(123)
  gamma[[1]][,2] = mvrnorm(n=1,rep(0,dim(Gamma[[1]][[1]])[1]),Gamma[[1]][[2]])
  set.seed(231)
  gamma[[2]][,1] = mvrnorm(n=1,rep(0,dim(Gamma[[2]][[1]])[1]),Gamma[[2]][[1]])
  set.seed(321)
  gamma[[2]][,2] = mvrnorm(n=1,rep(0,dim(Gamma[[2]][[1]])[1]),Gamma[[2]][[2]])
  
  location$gamma1 = rep(NA, curve.num)
  location$gamma1[cluster_true==1] = gamma[[1]][,1]
  location$gamma1[cluster_true==2] = gamma[[2]][,1]
  
  location$gamma2 = rep(NA, curve.num)
  location$gamma2[cluster_true==1] = gamma[[1]][,2]
  location$gamma2[cluster_true==2] = gamma[[2]][,2]
  
  for( i in 1:curve.num){
    set.seed(132)
    epsilon = rnorm(n=day,0,sigma_true)
    dd[i,] = mu[[cluster_true[i]]] + 
      as.vector(b %*% Theta_true[[cluster_true[i]]] %*% t(location[i,c('gamma1','gamma2')])) +
      epsilon
  }  
  data_list[[ite]] = dd
}

## save.image("hetero_data.Rdata") 
## It would be better to save these simulated data first



