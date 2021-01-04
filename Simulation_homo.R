##########################################################################
##  Spatial-Functional Mixture Models: Simulation of Homoscedastic Case ##
##  Outputs: Tables 1 and 2 of the manuscript.                          ##
##########################################################################
rm(list=ls())


##################################################################
##  Load associated required packages
library(fda)
library(funcy)
library(fields)
library(mclust)
library(MASS)
library(zoo)


##################################################################
##  Source files for the algorithm functions 
source('../../functions.R')
source('homo_iid_multi.R')
source('homo_iid_mrf.R')
source('homo_spatial_multi.R')
source('homo_spatial_mrf.R')



##################################################################
##  Simulate Homoscedastic Data             
##################################################################

##  location    
station = read.csv('../../../data/city_location_China.csv',sep=',',fileEncoding  = 'GBK')  
station = station[ 107<station$lon & station$lon<125 & 28<station$lat & station$lat<43, ]
location = station[,3:4]

curve.num = dim(station)[1]

##  B-spline
day = 30
t = seq(0,1,length.out=day)
basis = create.bspline.basis(c(0,1),nbasis=4)
s = getbasismatrix(t,basis)
b = qr.Q(qr(s)) ## orthogonalization

##  mean functions   
mu = list()
mu[[1]] = exp(t)*cos(t)/2
mu[[2]] = cos(5*pi/2*t)
cluster.num = length(mu)

lambda_true <- list()
for(i in 1:cluster.num){
  lambda_true[[i]] = solve(t(b)%*%b,t(b)%*%mu[[i]])
}

##  principal components  
Theta_ori= matrix(c(1,0,0,0,0,1,0,0),,2)
FPC1=sin(2*pi*t)       
Theta_ori[,1] = solve(t(b)%*%b,t(b)%*%FPC1) 

Theta_true = abs(qr.orth(Theta_ori))
FPC = b %*% Theta_true

##  parameters
sigma_gamma_true = c(7,2)
sigma_true = 0.4
phi_true = 1

##  Simulate markov random field
source('SimulateMRF.R')

knn = 5
dist=as.matrix(dist(location))
dist.ord=apply(dist,2,order)
dist.order=t(dist.ord)
nb=dist.order[,2:(knn+1)]

y.neigh = list()
for(i in 1:curve.num){
  y.neigh[[i]] = nb[i,]
}

##################################################################
##  simulation   
##################################################################

Nsimu = 100

## simulation result records 
randIndex = matrix(NA,Nsimu,7)
colnames(randIndex)=c('kmeans','Jiang','James','FMM','FMM-MRF','SFMM','SFMM-MRF')
RMSE = matrix(NA,Nsimu,7)
colnames(randIndex)=c('kmeans','Jiang','James','FMM','FMM-MRF','SFMM','SFMM-MRF')
phi_list = matrix(NA,Nsimu,2)
theta_list = matrix(NA,Nsimu,2)
sigma_gamma_list = matrix(NA,Nsimu,2)
data=list()
cluster=list()

P = 6
aa <- Sys.time()
for(ite in 1:Nsimu){
  print(paste("Simulation No.",ite,sep=''))
#  Sys.sleep(30)
  
  cluster_true = mrf[,ite]
  cluster[[ite]] = cluster_true
  
  Gamma = list()
  Gamma[[1]] =  sigma_gamma_true[1] * Exp.cov(location,theta=phi_true)
  Gamma[[2]] =  sigma_gamma_true[2] * Exp.cov(location,theta=phi_true)
  
  dd = matrix(NA,curve.num,day)

  set.seed(321)  
  gamma = cbind(mvrnorm(n=1,rep(0,dim(Gamma[[1]])[1]),Gamma[[1]]),mvrnorm(n=1,rep(0,dim(Gamma[[2]])[1]),Gamma[[2]]))
  
  for( i in 1:curve.num){
    ## noise 
    set.seed(231)  
    epsilon = rnorm(n=day,0,sigma_true)
    
    dd[i,] = mu[[cluster_true[i]]] + as.vector(b%*%Theta_true%*%gamma[i,]) + epsilon
  }
  
  # ########### plot ##########
  # 
  # plot(dd[1,],type='l',col=cluster_true[1],ylim=c(-7,7))
  # for(i in 2:curve.num){
  #   lines(dd[i,],col=cluster_true[i])
  # }
  
  data[[ite]]=dd
  
  
  ####################    k-means    ##############
  fit1 = kmeans(dd,cluster.num,iter.max=100)
  randIndex[ite,1] = round(adjustedRandIndex(fit1$cluster,cluster_true),4)
  
  muhat = list()
  rms = 0
  for(k in 1:cluster.num){
    muhat[[k]]=fit1$centers[k,]
    tmp = rep(NA,cluster.num)
    for(j in 1:cluster.num){
      tmp[j] = rmse(muhat[[k]],mu[[j]])
    }
    rms=rms+min(tmp)
  }
  RMSE[ite,1] = round(rms,4)
  
  ####################    Serban (Jiang)    #################
  
  fit2 = funcit(data=t(dd),method='fscm',knn=5,k=cluster.num,
                location=location,clusters=cluster_true)
  randIndex[ite,2] = fit2@randIndex
  
  muhat = list()
  rms = 0
  for(k in 1:cluster.num){
    muhat[[k]]=fit2@models$fscm@centers[,k]
    tmp = rep(NA,cluster.num)
    for(j in 1:cluster.num){
      tmp[j] = rmse(muhat[[k]],mu[[j]])
    }
    rms=rms+min(tmp)
  }
  RMSE[ite,2] = round(rms,4)
  
  
  ######################  James  ####################
  
  fit3=funcit(data=t(dd),method='fitfclust',clusters=cluster_true,k=cluster.num)
  
  randIndex[ite,3] = round(adjustedRandIndex(fit3@votedCluster,cluster_true),4)
  
  muhat = list()
  rms = 0
  for(k in 1:cluster.num){
    muhat[[k]]=fit3@models$fitfclust@centers[,k]
    tmp = rep(NA,2)
    for(j in 1:cluster.num){
      tmp[j] = rmse(muhat[[k]],mu[[j]])
    }
    rms=rms+min(tmp)
  }
  RMSE[ite,3] = round(rms,4)
  

  
  ######################  fcm (FMM)  ####################
  
  fit4 = homo_iid_multi(data=t(dd),cluster.num=cluster.num,P)
  
  randIndex[ite,4] = round(adjustedRandIndex(fit4$cluster,cluster_true),4)
  
  muhat2 = list()
  rms = 0
  for(k in 1:cluster.num){
    muhat2[[k]] = apply(dd[fit4$cluster==k,],2,mean)
    tmp = rep(NA,cluster.num)
    for(j in 1:cluster.num){
      tmp[j] = rmse(muhat2[[k]],mu[[j]])
    }
    rms=rms+min(tmp)
  }
  RMSE2 = round(rms,4)
  
  RMSE[ite,4] = RMSE2

  ######################  mfcm (FMM-MRF) ####################
  
  fit5 = homo_iid_mrf(data=t(dd),cluster.num=cluster.num,neigh=y.neigh,P)
  
  randIndex[ite,5] = round(adjustedRandIndex(fit5$cluster,cluster_true),4)
  
  muhat2 = list()
  rms = 0
  for(k in 1:cluster.num){
    muhat2[[k]] = apply(dd[fit5$cluster==k,],2,mean)
    tmp = rep(NA,cluster.num)
    for(j in 1:cluster.num){
      tmp[j] = rmse(muhat2[[k]],mu[[j]])
    }
    rms=rms+min(tmp)
  }
  RMSE2 = round(rms,4)
  
  RMSE[ite,5] = RMSE2
  theta_list[ite,1] = fit5$theta
  
  ######################  fcm2  ####################
  
  fit6 = homo_iid_multi2(data=t(dd),cluster.num=cluster.num,fit4=fit4,P)
  
  # randIndex[ite,6] = round(adjustedRandIndex(fit6$cluster,cluster_true),4)
  # 
  # muhat2 = list()
  # rms = 0
  # for(k in 1:cluster.num){
  #   muhat2[[k]] = apply(dd[fit6$cluster==k,],2,mean)
  #   tmp = rep(NA,cluster.num)
  #   for(j in 1:cluster.num){
  #     tmp[j] = rmse(muhat2[[k]],mu[[j]])
  #   }
  #   rms=rms+min(tmp)
  # }
  # RMSE2 = round(rms,4)
  # 
  # RMSE[ite,6] = RMSE2
  # sigma_list[ite,6] = fit6$sigma    
  
  ######################  fscm (SFMM)  ####################
  
  fit7 = homo_spatial_multi2(data=t(dd),cluster.num=cluster.num,location=location,fit4=fit4,fit6=fit6,P)
  
  randIndex[ite,6] = round(adjustedRandIndex(fit7$cluster,cluster_true),4)
  
  muhat2 = list()
  rms = 0
  for(k in 1:cluster.num){
    muhat2[[k]] = apply(dd[fit7$cluster==k,],2,mean)
    tmp = rep(NA,cluster.num)
    for(j in 1:cluster.num){
      tmp[j] = rmse(muhat2[[k]],mu[[j]])
    }
    rms=rms+min(tmp)
  }
  RMSE2 = round(rms,4)
  
  RMSE[ite,6] = RMSE2
  phi_list[ite,1]=fit7$phi

  #######################################################
  ###############     mfscm (SFMM-MRF)    ###############
  #######################################################
  
  fit8 = homo_spatial_mrf2(data=t(dd),cluster.num=cluster.num,location=location,neighbor=y.neigh,fit6=fit6,P)
  
  randIndex[ite,7] = round(adjustedRandIndex(fit8$cluster,cluster_true),4)
  
  muhat2 = list()
  rms = 0
  for(k in 1:cluster.num){
    muhat2[[k]] = apply(dd[fit8$cluster==k,],2,mean)
    tmp = rep(NA,cluster.num)
    for(j in 1:cluster.num){
      tmp[j] = rmse(muhat2[[k]],mu[[j]])
    }
    rms=rms+min(tmp)
  }
  RMSE2 = round(rms,4)
  
  RMSE[ite,7] = RMSE2
  phi_list[ite,2]=fit8$phi
  theta_list[ite,2]=fit8$theta
  sigma_gamma_list[ite,1]=fit8$sigma_gamma[1]
  sigma_gamma_list[ite,2]=fit8$sigma_gamma[2]
  
  
  print(Sys.time()-aa)
  Sys.sleep(30)
  aa=Sys.time()
}



##################################################################
## Reproduce Table 1 and Table 2 in the manuscript
##################################################################

#-----------------------------------------------------------------
# Table 1: 
# Means and standard deviations of adjusted Rand index and RMSE 
#-----------------------------------------------------------------

round(apply(randIndex,2,mean),3)   ##  RandIndex
round(apply(randIndex,2,sd),3)   

round(apply(RMSE,2,mean),3)  ##  RMSE
round(apply(RMSE,2,sd),3)

#-----------------------------------------------------------------
# Table 2: 
# Means and standard deviations of parameters using SFMM-MRF
#-----------------------------------------------------------------

para_list = cbind(phi_list[,2],theta_list[,2],sigma_gamma_list)
colnames(para_list) = c('phi','theta','sigma_gamma_1','sigma_gamma_2')

round( apply(para_list,2,mean), 3)
round( apply(para_list,2,sd), 3)










