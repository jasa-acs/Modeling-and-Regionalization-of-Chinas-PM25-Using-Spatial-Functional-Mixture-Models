##########################################################################
##  Spatial-Functional Mixture Models: Simulation of Homoscedastic Case ##
##  Outputs: Table 3, Table 4, and Figure 4  of the manuscript.         ##
##########################################################################
rm(list=ls())

##################################################################
##  Load associated required packages
library(fda)
library(fields)
library(MASS)
library(mclust)
library(zoo)
library(doParallel)
library(foreach)

##################################################################
##  Source files for the algorithm functions 
source('../../functions.R')

## log likelihood for hfscm (HSFMM) 
loglike = function(y.data,O1, gamma_list,Z_list,lambda,Theta,phi,sigma_gamma,sigma,station,prob) {
  
  curve.num = length(y.data)
  day = length(y.data[[1]])
  cluster.num = dim(sigma_gamma)[1]
  Q = dim(sigma_gamma)[2]
  MCC = dim(Z_list)[2]
  
  zero = function(x){ze = rep(0,curve.num);ze[x]=1;return(ze)}
  
  s1_list = rep(NA,MCC)   
  s2_list = rep(NA,MCC)
  s3_list = rep(NA,MCC)
  for(t in 1:MCC){	
    z = Z_list[,t]
    gamma = gamma_list[,t]
    
    s1 = 0; s3 = 0
    
    #################    s1    ###############
    for(i in 1:curve.num){
      for(k in 1:cluster.num){
        s1 = s1 - 0.5 * z_ik(z)[k,i] * (day*log(sigma) + t( y.data[[i]] - b %*% (lambda[[k]] + Theta[[k]] %*% gamma[((i-1)*Q+1):(i*Q)]) ) %*% (y.data[[i]] - b %*% (lambda[[k]] + Theta[[k]] %*% gamma[((i-1)*Q+1):(i*Q)])) / (sigma) )
      }
    }
    
    #################    s2    ###############	
    
    mem = which(z==1)
    for( i in 2:cluster.num){ mem = c(mem,which(z==i))}
    
    o2 = matrix(unlist(lapply(as.list(mem),'zero')),curve.num,)
    O2 = as.matrix(do.call('bdiag',rep(list(o2),Q)))
    O = O1 %*% O2
    
    Gamma1 = as.matrix(do.call('bdiag',mapply(FUN='*',as.list(sigma_gamma),lapply(as.list(rep(1:cluster.num,2)),function(x){return(P.func(station[which(z==x),],phi = phi))}),SIMPLIFY=F)))
    
    Gamma = O %*% Gamma1 %*% t(O)          
    
    s2 =  -0.5*(log(det(Gamma))+t(gamma)%*%solve(Gamma)%*%gamma)  
    
    #################    s3    ###############	
    for(i in 1: curve.num){
      for(k in 1:cluster.num){ s3 = s3 + z_ik(z)[k,i] * log(prob[k]) }
    }    
    
    s1_list[t] = s1
    s2_list[t] = s2
    s3_list[t] = s3
  }
  s11 = mean(s1_list)
  s22 = mean(s2_list)
  s33 = mean(s3_list)
  print(c(s1,s2,s3))
  
  return(s11+s22+s33)
}

## log likelihood for hmfscm (HSFMM-MRF) 
loglike_mrf = function(y.data,O1, gamma_list,Z_list,lambda,Theta,phi,sigma_gamma,sigma,station,theta,y.neigh) {
  
  curve.num = length(y.data)
  day = length(y.data[[1]])
  cluster.num = dim(sigma_gamma)[1]
  Q = dim(sigma_gamma)[2]
  MCC = dim(Z_list)[2]
  
  
  zero = function(x){ze = rep(0,curve.num);ze[x]=1;return(ze)}
  
  s1_list = rep(NA,MCC)   
  s2_list = rep(NA,MCC)
  s3_list = rep(NA,MCC)
  for(t in 1:MCC){	
    z = Z_list[,t]
    gamma = gamma_list[,t]
    
    s1 = 0; s3 = 0
    
    #################    s1    ###############
    for(i in 1:curve.num){
      for(k in 1:cluster.num){
        s1 = s1 - 0.5 * z_ik(z)[k,i] * (day*log(sigma) + t( y.data[[i]] - b %*% (lambda[[k]] + Theta[[k]] %*% gamma[((i-1)*Q+1):(i*Q)]) ) %*% (y.data[[i]] - b %*% (lambda[[k]] + Theta[[k]] %*% gamma[((i-1)*Q+1):(i*Q)])) / (sigma) )
      }
    }
    
    #################    s2    ###############	
    
    mem = which(z==1)
    for( i in 2:cluster.num){ mem = c(mem,which(z==i))}
    
    o2 = matrix(unlist(lapply(as.list(mem),'zero')),curve.num,)
    O2 = as.matrix(do.call('bdiag',rep(list(o2),Q)))
    O = O1 %*% O2
    
    Gamma1 = as.matrix(do.call('bdiag',mapply(FUN='*',as.list(sigma_gamma),lapply(as.list(rep(1:cluster.num,2)),function(x){return(P.func(station[which(z==x),],phi = phi))}),SIMPLIFY=F)))
    
    Gamma = O %*% Gamma1 %*% t(O)          
    
    s2 =  -0.5*(log(det(Gamma))+t(gamma)%*%solve(Gamma)%*%gamma)  
    
    #################    s3    ###############	
    for(i in 1: curve.num){
      for(k in 1:cluster.num){ s3 = s3 + z_ik(z)[k,i] * ( theta * sum(z_ik(z)[k,y.neigh[[i]]]) - log( sum(apply(as.matrix(z_ik(z)[,y.neigh[[i]]]),1,function(t){exp(theta*sum(t))})) )) }
    }    
    
    s1_list[t] = s1
    s2_list[t] = s2
    s3_list[t] = s3
  }
  s11 = mean(s1_list)
  s22 = mean(s2_list)
  s33 = mean(s3_list)
  print(c(s1,s2,s3))
  
  return(s11+s22+s33)
}

##################################################################
##  Simulate heteroscedastic data
##################################################################

source('SimulateHeteroData.R')
## It would be better to use: load ("hetero_data.Rdata") 
## and then perform the following simulations separately 
## after save.image("hetero_data.Rdata") in the "Simulate_HeteroData.R"


#################################################################
#################   fit hmfscm (HSFMM-MRF) models   #############
#################################################################

Nsimu = 100

#detectCores()
cl = makeCluster(50) 
registerDoParallel(cl)
options(warn=-1)

results_HSFMM_MRF = foreach(ite=1:Nsimu,.packages=c('foreach','Matrix','MASS','parallel','fields','fda','funcy','mclust','zoo'))%dopar%
{
  
  cluster_true = mrf[,ite]
  
  dd = data_list[[ite]]
  
  fit1 = kmeans(dd,cluster.num,iter.max=100)  ####  kmeans
  
  #############     Hmfscm (HSFMM-MRF)    ###############

  tryCatch({ 
    x = t(dd)
    curve.num = dim(x)[2]
    
    zero = function(x){ze = rep(0,curve.num);ze[x]=1;return(ze)}
    
    y.data = as.list(as.data.frame(na.approx(x)))
    Y = as.matrix(unlist(y.data))
    day = dim(x)[1]
    
    t = seq(0,1,length.out=day)
    basis = create.bspline.basis(c(0,1),nbasis=4)
    s = getbasismatrix(t,basis)
    
    P = 4
    Q = 2
    O1 = NULL
    for(j in 1:curve.num){
      ttmp = list()
      for(i in 1:Q){
        a = rep(0,curve.num*Q)
        a[curve.num*(i-1)+j] = 1
        ttmp[[i]] = a
      }
      tmp0 = as.matrix(ttmp[[1]])
      for(i in 2:Q){ tmp0 = cbind(tmp0,as.matrix(ttmp[[i]])) }
      tmp0=t(tmp0)
      O1 = rbind(O1,tmp0)
    }
    #### orthogonalization
    b = qr.Q(qr(s))
    B = b
    for(i in 2:curve.num){ B = as.matrix(bdiag(B,b))}
    
    ######## neigh ########
    
    knn=5
    dist=as.matrix(dist(location))
    dist.ord=apply(dist,2,order)
    dist.order=t(dist.ord)
    nb=dist.order[,2:(knn+1)]
    
    y.neigh = list()
    for(i in 1:curve.num){
      y.neigh[[i]] = nb[i,]
    }
    
    ###########   initial values  ###########

    theta = 1      
    
    con.prob <- z_ik(fit1$cluster)
    
    max.con.prob <- apply(con.prob,2,which.max)
    
    mar.prob <- matrix(NA,cluster.num,curve.num)
    for(i in 1:curve.num){
      for(k in 1:cluster.num){
        tmp <- apply(as.matrix(con.prob[,y.neigh[[i]]]),1,function(t){exp(theta*sum(t))})
        mar.prob[k,i]<- tmp[k]/sum(tmp)
      }
    }    
    ##### random Para : lambda, sigma, theta
    
    lambda <- list()
    for(i in 1:cluster.num){
      lambda[[i]] = solve(t(b)%*%b,t(b)%*%fit1$centers[i,])
    }
    
    sigma <- 1
    
    ############
    Theta = NULL
    for(k in 1:cluster.num){
      Theta[[k]] = rbind(diag(Q),matrix(0,P-Q,Q))
    }         
    
    
    naive <- c(2,1.5,1,0.5,0.3,0.2,0.1,0.05,0.02,0.01)
    sigma_gamma = matrix(rep(naive[1:Q],cluster.num),cluster.num,Q,byrow=T)
    phi = 1
    
    
    #####################
    a <- Sys.time()
    LOOP = 35
    loglike_list = rep(NA,LOOP)
    conv = rep(0,8)
    for(loop in 1:LOOP){
      print(paste('iteration',loop))
      ########## E-step ##########

      #### Monto-Carlo EM
      MC = 100
      MCC = 75
      
      ######## gamma
      gamma_list = matrix(NA,curve.num*Q,MC)
      Z_list = matrix(NA,curve.num,MC)
      
      ######## initial z and its sons
      z = max.con.prob
      
      e.Alpha = do.call('rbind',lapply(as.list(z),function(x){unlist(lambda[[x]])}))
      
      e.S = as.matrix(do.call('bdiag',lapply(as.list(z),function(x){b%*%Theta[[x]]})))
      
      e.mem = which(z==1)
      for( i in 2:cluster.num){ e.mem = c(e.mem,which(z==i))}
      
      e.o2 = matrix(unlist(lapply(as.list(e.mem),'zero')),curve.num,)
      e.O2 = as.matrix(do.call('bdiag',rep(list(e.o2),Q)))
      e.O = O1 %*% e.O2
      
      e.Gamma1 = as.matrix(do.call('bdiag',mapply(FUN='*',as.list(sigma_gamma),lapply(as.list(rep(1:cluster.num,2)),function(x){return(P.func(location[which(z==x),],phi = phi))}),SIMPLIFY=F)))
      
      e.Gamma = e.O %*% e.Gamma1 %*% t(e.O)
      
      t2 = Sys.time()
      print( paste('Monta-Carlo EM') )     
      for( t in 1:MC){
        
        ############ gamma
        
        #### var(Y)
        v_Y = e.S %*%e.Gamma%*%t(e.S) + diag(sigma,curve.num*day)  ### computation
        #solve(v_Y)
        ### Woodbury inverse
        inv_v_Y = diag(1/sigma,curve.num * day)-1/(sigma^2) * e.S %*% solve((solve(e.Gamma)+t(e.S) %*% e.S/sigma))%*%t(e.S)
        
        #### expectation
        E.gamma = e.Gamma %*% t(e.S) %*% inv_v_Y %*% (Y - B %*% e.Alpha)
        #### variance
        V.gamma = e.Gamma - e.Gamma %*% t(e.S) %*% inv_v_Y %*%  e.S %*% e.Gamma
        
        #### Gibbs sampling
        gamma = mvrnorm(1,E.gamma,V.gamma)
        gamma_list[,t] = gamma
        
        ############## z
        
        for(i in 1:curve.num){
          gamma.i = gamma[((i-1)*Q+1):(i*Q)]
          tmp = rep(0,cluster.num)
          p = rep(0,cluster.num)
          tmp1 =  mapply(function(x,y){-0.5*t(y.data[[i]]-b%*%x-b%*%y%*%gamma.i)%*%
              (y.data[[i]]-b%*%x-b%*%y%*%gamma.i)/sigma},lambda,Theta)
          tmp1 = unlist(tmp1)
          if((max(tmp1)-min(tmp1))>500){tmp[which.max(tmp1)]=1}else{      tmp = exp(tmp1)*mar.prob[,i]}
          for(k in 1:cluster.num){  p[k] <- tmp[k]/sum(tmp) }  
          z[i] = sample(1:cluster.num,1, prob = p)  
        }
        Z_list[,t] = z
        
        e.Alpha = do.call('rbind',lapply(as.list(z),function(x){unlist(lambda[[x]])}))
        
        e.S = as.matrix(do.call('bdiag',lapply(as.list(z),function(x){b%*%Theta[[x]]})))
        
        e.mem = which(z==1)
        for( i in 2:cluster.num){ e.mem = c(e.mem,which(z==i))}
        
        e.o2 = matrix(unlist(lapply(as.list(e.mem),'zero')),curve.num,)
        e.O2 = as.matrix(do.call('bdiag',rep(list(e.o2),Q)))
        e.O = O1 %*% e.O2
        
        e.Gamma1 = as.matrix(do.call('bdiag',mapply(FUN='*',as.list(sigma_gamma),lapply(as.list(rep(1:cluster.num,2)),function(x){return(P.func(location[which(z==x),],phi = phi))}),SIMPLIFY=F)))
        
        e.Gamma = e.O %*% e.Gamma1 %*% t(e.O)          
        
        #t2=Sys.time()
        #print(Sys.time()-t2)  
      }
      print(Sys.time()-t2)                  
      
#      MCC = 30
      
      index = apply(Z_list,2,function(x){length(levels(factor(x)))})==cluster.num
      Z_list = Z_list[,index]
      Z_list = Z_list[,(dim(Z_list)[2]-(MCC-1)):dim(Z_list)[2]]
      gamma_list = gamma_list[,index]
      gamma_list = gamma_list[,(dim(gamma_list)[2]-(MCC-1)):dim(gamma_list)[2]]
      
      
      ##### con.prob : z_ik | Y_i : maybe marginal from Monto Carlo EM
      
      
      for(k in 1:cluster.num){
        for(i in 1:curve.num){
          con.prob[k,i] = sum(Z_list[i,]==k)/MCC
        }
      }	
      print(con.prob)  
      
      #############  M-step ############

      max1 = apply(con.prob,2,which.max)
      print(paste('z','diffrence',norm(max.con.prob-max1)))
      max.con.prob <- max1
      
      ######## theta
      theta_func <- function(x){
        tmp <- matrix(NA, cluster.num, curve.num)
        for(i in 1:curve.num){
          for( k in 1:cluster.num){
            tmp1 <- apply(as.matrix(con.prob[,y.neigh[[i]]]),1,function(t){exp(x*sum(t))*sum(t)})
            tmp2 <- apply(as.matrix(con.prob[,y.neigh[[i]]]),1,function(t){exp(x*sum(t))})
            tmp[k,i] <- con.prob[k,i]*(sum(con.prob[k,y.neigh[[i]]])-sum(tmp1)/sum(tmp2))
          }
        }
        return(sum(tmp))
      }
      theta1 <-  tryCatch({ uniroot(theta_func,c(0,5))$root },error=function(e){cat(conditionMessage(e));
        theta=theta})
      
      a4 = abs(theta1- theta)/(theta+0.001)
      theta=theta1
      
      mar.prob <- matrix(NA,cluster.num,curve.num)
      for(i in 1:curve.num){
        for(k in 1:cluster.num){
          tmp <- apply(as.matrix(con.prob[,y.neigh[[i]]]),1,function(t){exp(theta*sum(t))})
          mar.prob[k,i]<- tmp[k]/sum(tmp)
        }
      }
      
      
      ########   lambda_k
      
      for(k in 1:cluster.num){
        tmp1 = sum(apply(Z_list,2,function(x){return(sum(z_ik(x)[k,]))}))
        
        tmp = matrix(0,P,1)
        for(t in 1:(MCC)){
          for(i in 1:curve.num){
            tmp = tmp + z_ik(Z_list[,t])[k,i]*t(b)%*%(y.data[[i]]-b%*%Theta[[k]]%*%gamma_list[((i-1)*Q+1):(i*Q),t])
          }
        }
        tmp = solve(tmp1*t(b)%*%b)%*%tmp
        a1 = sqrt(sum((tmp-lambda[[k]])^2)/P)/sqrt(sum(lambda[[k]]^2)/P+0.001)
        # a1 = norm(tmp-lambda[[k]])
        print(paste('lambda',k,'difference',a1))
        lambda[[k]] <- tmp
      }
      
      ########   sigma
      
      tmp1 = rep(0,(MCC))
      tmp2 = rep(0,(MCC))
      
      for(t in 1:(MCC)){
        for(i in 1:curve.num){
          for(k in 1:cluster.num){
            tmp1[t] = tmp1[t] + z_ik(Z_list[,t])[k,i]*t(y.data[[i]]-b%*%(lambda[[k]])-b%*%Theta[[k]]%*%gamma_list[((i-1)*Q+1):(i*Q),t]) %*% (y.data[[i]]-b%*%(lambda[[k]])-b%*%Theta[[k]]%*%gamma_list[((i-1)*Q+1):(i*Q),t])
            tmp2[t] = tmp2[t] + z_ik(Z_list[,t])[k,i]*day
          }
        }
      }
      sigma1 <- mean(tmp1)/mean(tmp2)
      a2 = abs(sigma-sigma1)/(abs(sigma)+0.001)
      print(paste('sigma','difference',a2))
      sigma = sigma1
      
      #########  Theta_k
      for(k in 1:cluster.num){
        for(q in 1:Q){
          tmp1 = 0
          tmp = rep(0,day)
          for(t in 1:(MCC)){
            for(i in 1:curve.num){
              tmp1 = tmp1 + z_ik(Z_list[,t])[k,i]*(gamma_list[((i-1)*Q+q),t])^2
              tmp  = tmp   + z_ik(Z_list[,t])[k,i]*(gamma_list[((i-1)*Q+q),t])* (y.data[[i]]-b%*%(lambda[[k]])-apply((b%*%Theta[[k]][,-q]) %*% (diag(gamma_list[((i-1)*Q+1):(i*Q),t][-q],nrow=length(gamma_list[((i-1)*Q+1):(i*Q),t][-q]))),1,sum) )
            }
          }
          tmp = solve(tmp1*t(b)%*%b) %*% t(b) %*% tmp
          print(paste('Theta',k,q,'difference',norm(Theta[[k]][,q]-tmp)))
          Theta[[k]][,q] = tmp
        }
        Theta[[k]] = qr.orth(Theta[[k]])
        
      }        
      
      ######## sigma_gamma
      
      for(k in 1:cluster.num){
        for(q in 1:Q){
          tmp = rep(0,(MCC))
          for(t in 1:(MCC)){
            #print(t)
            z = Z_list[,t]
            Nz = rep(NA,cluster.num)
            for(i in 1:cluster.num){ Nz[i] = sum(z==i)}
            mem = which(z==1)
            for( i in 2:cluster.num){ mem = c(mem,which(z==i))}
            o2 = matrix(unlist(lapply(as.list(mem),'zero')),curve.num,)
            O2 = as.matrix(do.call('bdiag',rep(list(o2),Q)))
            O = O1 %*% O2
            gamma1 = ((solve(O)%*%gamma_list[,t])[(curve.num*(q-1)+1):(curve.num*q),1])[(if(k>1) sum(Nz[1:(k-1)])+1 else 1):sum(Nz[1:k])] 
            tmp[t] = t(as.matrix(gamma1)) %*% solve(P.func(phi,location[which(z==k),])) %*% as.matrix(gamma1)/Nz[k]			
          }
          print(paste('sigma_gamma',k,q,'difference',abs(mean(tmp)-sigma_gamma[k,q])))
          sigma_gamma[k,q] = mean(tmp)
        }
      }
      
      
      ######## phi        
      d.Obj <- function(phi){
        latent = rbind(Z_list,gamma_list)
        
        obj = function(x){
          z = x[1:curve.num] ; gamma = x[(curve.num+1):length(x)]
          Nz = rep(NA,cluster.num)
          for(i in 1:cluster.num){ Nz[i] = sum(z==i)}
          mem = which(z==1)
          for( i in 2:cluster.num){ mem = c(mem,which(z==i))}
          o2 = matrix(unlist(lapply(as.list(mem),'zero')),curve.num,)
          O2 = as.matrix(do.call('bdiag',rep(list(o2),Q)))
          O = O1 %*% O2
          
          tmp1 = 0
          tmp2 = 0
          for(k in 1:cluster.num){
            tmp1 = tmp1 + tr(solve(P.func(phi,location[which(z==k),]))%*%d.P.func(phi,location[which(z==k),]))
            for(q in 1:Q){
              gamma1 = ((solve(O)%*%gamma)[(curve.num*(q-1)+1):(curve.num*q),1])[(if(k>1) sum(Nz[1:(k-1)])+1 else 1):sum(Nz[1:k])] 
              
              tmp2 = tmp2 + tr(solve(P.func(phi,location[which(z==k),]))%*% d.P.func(phi,location[which(z==k),]) %*% solve(P.func(phi,location[which(z==k),])) %*% as.matrix(gamma1) %*% t(as.matrix(gamma1)))/sigma_gamma[k,q]
            }
          }
          return(Q * tmp1 - tmp2)            
        }
        return(mean(apply(latent,2,obj)))
      }
      
      
      print('optimization for phi')
      t3 = Sys.time()  
      phi1 <- tryCatch({  uniroot(d.Obj,c(0.01,5))$root},error=function(e){cat(conditionMessage(e));print('sb');return(phi)})
      print(Sys.time()-t3)
      a3 = abs(phi-phi1)/(abs(phi)+0.001)
      print(paste('phi','difference',a3))
      phi = phi1      
      
      ######### log likelihood
      loglike_list[loop] = loglike_mrf(y.data,O1, gamma_list,Z_list,lambda,Theta,phi,sigma_gamma,sigma,location,theta,y.neigh)    
      
      
      #     if(norm(max.con.prob-max1)==0 & a1<0.005 &a2<0.005 & a3<0.005 &a4<0.005){conv = conv+1}
      if(a1<0.004){conv[1]=conv[1]+1}else{conv[1]=0}
      if(a2<0.001){conv[2]=conv[2]+1}else{conv[2]=0}
      if(a3<0.001){conv[3]=conv[3]+1}else{conv[3]=0}
      if(a4<0.001){conv[4]=conv[4]+1}else{conv[4]=0}
      if(a1<0.01){conv[5]=conv[5]+1}else{conv[5]=0}
      if(a2<0.005){conv[6]=conv[6]+1}else{conv[6]=0}
      if(a3<0.005){conv[7]=conv[7]+1}else{conv[7]=0}
      if(a4<0.005){conv[8]=conv[8]+1}else{conv[8]=0}
      
      #     if(max(a1,a2,a3,a4)<0.005) break
      print(Sys.time()-a)
      a <- Sys.time()
      print(paste('sigma=',sigma,'phi=',phi,'loglike=',loglike_list[loop],'conv='))
      print(conv)
      print(sigma_gamma)      
      
      print(round(adjustedRandIndex(max.con.prob,cluster_true),4))
      
      Sys.sleep(30)
      
      write.csv(randIndex,paste(ite,'-',loop,'.csv',sep=''))
    }
    
    randIndex = round(adjustedRandIndex(max.con.prob,cluster_true),4)
    cluster=max.con.prob
    
    muhat = list()
    rms = 0
    for(k in 1:cluster.num){
      muhat[[k]]=as.vector(b%*%lambda[[k]])
      tmp = rep(NA,cluster.num)
      for(j in 1:cluster.num){
        tmp[j] = rmse(muhat[[k]],mu[[j]])
      }
      rms=rms+min(tmp)
    }
    RMSE1 = round(rms,4)
    
    muhat2 = list()
    rms = 0
    for(k in 1:cluster.num){
      muhat2[[k]] = apply(dd[max.con.prob==k,],2,mean)
      tmp = rep(NA,cluster.num)
      for(j in 1:cluster.num){
        tmp[j] = rmse(muhat2[[k]],mu[[j]])
      }
      rms=rms+min(tmp)
    }
    RMSE2 = round(rms,4)
    
    
    hmfscm = list(conv=conv,cluster=cluster,phi=phi,sigma=sigma,sigma_gamma=sigma_gamma,
                  con.prob=con.prob,gamma_list=gamma_list,Z_list=Z_list,Theta=Theta,lambda=lambda,
                  theta=theta,y.neigh=y.neigh,loglike_list=loglike_list)
    
    return(list(randIndex=randIndex,RMSE=RMSE1,RMSE2=RMSE2,
                hmfscm = hmfscm,cluster_true=cluster_true))
    
    
  },error=function(e){cat(conditionMessage(e));return(NA)})
  
}  

stopCluster(cl)
#save.image("results\\Hmfscm_results1.Rdata")


#################################################################
#################     fit hfscm (HSFMM) models     ##############
#################################################################

Nsimu = 100

#detectCores()
cl = makeCluster(3)
registerDoParallel(cl)
options(warn=-1)

results_HSFMM =foreach(ite=index,.packages=c('foreach','Matrix','MASS','parallel','fields','fda','funcy','mclust','zoo'))%dopar%
{
  
  cluster_true = mrf[,ite]
  
  dd = data_list[[ite]]
  
  fit1 = kmeans(dd,cluster.num,iter.max=100)  ####  kmeans
  
  ####################     Hfscm     ####################

  tryCatch({ 
    x = t(dd)
    curve.num = dim(x)[2]
    
    zero = function(x){ze = rep(0,curve.num);ze[x]=1;return(ze)}
    
    y.data = as.list(as.data.frame(na.approx(x)))
    Y = as.matrix(unlist(y.data))
    day = dim(x)[1]
    
    t = seq(0,1,length.out=day)
    basis = create.bspline.basis(c(0,1),nbasis=4)
    s = getbasismatrix(t,basis)
    
    P = 4
    Q = 2
    O1 = NULL
    for(j in 1:curve.num){
      ttmp = list()
      for(i in 1:Q){
        a = rep(0,curve.num*Q)
        a[curve.num*(i-1)+j] = 1
        ttmp[[i]] = a
      }
      tmp0 = as.matrix(ttmp[[1]])
      for(i in 2:Q){ tmp0 = cbind(tmp0,as.matrix(ttmp[[i]])) }
      tmp0=t(tmp0)
      O1 = rbind(O1,tmp0)
    }
    #### orthogonalization
    b = qr.Q(qr(s))
    B = b
    for(i in 2:curve.num){ B = as.matrix(bdiag(B,b))}
    
    ###########   initial values  ###########

    prob = rep(1/cluster.num,cluster.num)
    
    con.prob <- z_ik(fit1$cluster)
    
    max.con.prob <- apply(con.prob,2,which.max)
    
    ##### random Para : lambda, sigma, theta
    
    lambda <- list()
    for(i in 1:cluster.num){
      lambda[[i]] = solve(t(b)%*%b,t(b)%*%fit1$centers[i,])
    }
    
    sigma <- 1
    
    ############
    Theta = NULL
    for(k in 1:cluster.num){
      Theta[[k]] = rbind(diag(Q),matrix(0,P-Q,Q))
    }         
    
    ############  spatial Para ##############
    
    naive <- c(2,1.5,1,0.5,0.3,0.2,0.1,0.05,0.02,0.01)
    sigma_gamma = matrix(rep(naive[1:Q],cluster.num),cluster.num,Q,byrow=T)
    phi = 1
    
    
    #####################
    a <- Sys.time()
    LOOP = 60
    loglike_list = rep(NA,LOOP)
    conv = rep(0,8)
    for(loop in 1:LOOP){
      print(paste('iteration',loop))
      ########## E-step ##########

      
      #### Monto-Carlo EM
      MC = 150
      MCC = 100
      
      ######## gamma
      gamma_list = matrix(NA,curve.num*Q,MC)
      Z_list = matrix(NA,curve.num,MC)
      
      ######## initial z and its sons
      z = max.con.prob
      
      e.Alpha = do.call('rbind',lapply(as.list(z),function(x){unlist(lambda[[x]])}))
      
      e.S = as.matrix(do.call('bdiag',lapply(as.list(z),function(x){b%*%Theta[[x]]})))
      
      e.mem = which(z==1)
      for( i in 2:cluster.num){ e.mem = c(e.mem,which(z==i))}
      
      e.o2 = matrix(unlist(lapply(as.list(e.mem),'zero')),curve.num,)
      e.O2 = as.matrix(do.call('bdiag',rep(list(e.o2),Q)))
      e.O = O1 %*% e.O2
      
      e.Gamma1 = as.matrix(do.call('bdiag',mapply(FUN='*',as.list(sigma_gamma),lapply(as.list(rep(1:cluster.num,2)),function(x){return(P.func(location[which(z==x),],phi = phi))}),SIMPLIFY=F)))
      
      e.Gamma = e.O %*% e.Gamma1 %*% t(e.O)
      
      t2 = Sys.time()
      print( paste('Monta-Carlo EM') )     
      for( t in 1:MC){
        
        ############ gamma
        
        #### var(Y)
        v_Y = e.S %*%e.Gamma%*%t(e.S) + diag(sigma,curve.num*day)  ### computation
        #solve(v_Y)
        ### Woodbury inverse
        inv_v_Y = diag(1/sigma,curve.num * day)-1/(sigma^2) * e.S %*% solve((solve(e.Gamma)+t(e.S) %*% e.S/sigma))%*%t(e.S)
        
        #### expectation
        E.gamma = e.Gamma %*% t(e.S) %*% inv_v_Y %*% (Y - B %*% e.Alpha)
        #### variance
        V.gamma = e.Gamma - e.Gamma %*% t(e.S) %*% inv_v_Y %*%  e.S %*% e.Gamma
        
        #### Gibbs sampling
        gamma = mvrnorm(1,E.gamma,V.gamma)
        gamma_list[,t] = gamma
        
        ############## z
        
        for(i in 1:curve.num){
          gamma.i = gamma[((i-1)*Q+1):(i*Q)]
          tmp = rep(0,cluster.num)
          p = rep(0,cluster.num)
          tmp1 =  mapply(function(x,y){-0.5*t(y.data[[i]]-b%*%x-b%*%y%*%gamma.i)%*%
              (y.data[[i]]-b%*%x-b%*%y%*%gamma.i)/sigma},lambda,Theta)
          tmp1 = unlist(tmp1)
          if((max(tmp1)-min(tmp1))>500){tmp[which.max(tmp1)]=1}else{      tmp = exp(tmp1)*prob}
          for(k in 1:cluster.num){  p[k] <- tmp[k]/sum(tmp) }  
          z[i] = sample(1:cluster.num,1, prob = p)  
        }
        Z_list[,t] = z
        
        e.Alpha = do.call('rbind',lapply(as.list(z),function(x){unlist(lambda[[x]])}))
        
        e.S = as.matrix(do.call('bdiag',lapply(as.list(z),function(x){b%*%Theta[[x]]})))
        
        e.mem = which(z==1)
        for( i in 2:cluster.num){ e.mem = c(e.mem,which(z==i))}
        
        e.o2 = matrix(unlist(lapply(as.list(e.mem),'zero')),curve.num,)
        e.O2 = as.matrix(do.call('bdiag',rep(list(e.o2),Q)))
        e.O = O1 %*% e.O2
        
        e.Gamma1 = as.matrix(do.call('bdiag',mapply(FUN='*',as.list(sigma_gamma),lapply(as.list(rep(1:cluster.num,2)),function(x){return(P.func(location[which(z==x),],phi = phi))}),SIMPLIFY=F)))
        
        e.Gamma = e.O %*% e.Gamma1 %*% t(e.O)          
        
        #t2=Sys.time()
        #print(Sys.time()-t2)  
      }
      print(Sys.time()-t2)                  
      
      index = apply(Z_list,2,function(x){length(levels(factor(x)))})==cluster.num
      Z_list = Z_list[,index]
      Z_list = Z_list[,(dim(Z_list)[2]-(MCC-1)):dim(Z_list)[2]]
      gamma_list = gamma_list[,index]
      gamma_list = gamma_list[,(dim(gamma_list)[2]-(MCC-1)):dim(gamma_list)[2]]
      
      
      ##### con.prob : z_ik | Y_i : maybe marginal from Monto Carlo EM
      
      
      for(k in 1:cluster.num){
        for(i in 1:curve.num){
          con.prob[k,i] = sum(Z_list[i,]==k)/MCC
        }
      }	
      print(con.prob)  
      
      #############  M-step ############

      ########    prob
      
      prob1 = apply(con.prob,1,mean)   
      a4 = abs(prob[1]- prob1[1])/(prob[1]+0.001)
      print(paste('pi','difference',a4))
      prob = prob1
      
      ########   lambda_k
      
      for(k in 1:cluster.num){
        tmp1 = sum(apply(Z_list,2,function(x){return(sum(z_ik(x)[k,]))}))
        
        tmp = matrix(0,P,1)
        for(t in 1:(MCC)){
          for(i in 1:curve.num){
            tmp = tmp + z_ik(Z_list[,t])[k,i]*t(b)%*%(y.data[[i]]-b%*%Theta[[k]]%*%gamma_list[((i-1)*Q+1):(i*Q),t])
          }
        }
        tmp = solve(tmp1*t(b)%*%b)%*%tmp
        a1 = sqrt(sum((tmp-lambda[[k]])^2)/P)/sqrt(sum(lambda[[k]]^2)/P+0.001)
        # a1 = norm(tmp-lambda[[k]])
        print(paste('lambda',k,'difference',a1))
        lambda[[k]] <- tmp
      }
      
      ########   sigma
      
      tmp1 = rep(0,(MCC))
      tmp2 = rep(0,(MCC))
      
      for(t in 1:(MCC)){
        for(i in 1:curve.num){
          for(k in 1:cluster.num){
            tmp1[t] = tmp1[t] + z_ik(Z_list[,t])[k,i]*t(y.data[[i]]-b%*%(lambda[[k]])-b%*%Theta[[k]]%*%gamma_list[((i-1)*Q+1):(i*Q),t]) %*% (y.data[[i]]-b%*%(lambda[[k]])-b%*%Theta[[k]]%*%gamma_list[((i-1)*Q+1):(i*Q),t])
            tmp2[t] = tmp2[t] + z_ik(Z_list[,t])[k,i]*day
          }
        }
      }
      sigma1 <- mean(tmp1)/mean(tmp2)
      a2 = abs(sigma-sigma1)/(abs(sigma)+0.001)
      print(paste('sigma','difference',a2))
      sigma = sigma1
      
      #########  Theta_k
      for(k in 1:cluster.num){
        for(q in 1:Q){
          tmp1 = 0
          tmp = rep(0,day)
          for(t in 1:(MCC)){
            for(i in 1:curve.num){
              tmp1 = tmp1 + z_ik(Z_list[,t])[k,i]*(gamma_list[((i-1)*Q+q),t])^2
              tmp  = tmp   + z_ik(Z_list[,t])[k,i]*(gamma_list[((i-1)*Q+q),t])* (y.data[[i]]-b%*%(lambda[[k]])-apply((b%*%Theta[[k]][,-q]) %*% (diag(gamma_list[((i-1)*Q+1):(i*Q),t][-q],nrow=length(gamma_list[((i-1)*Q+1):(i*Q),t][-q]))),1,sum) )
            }
          }
          tmp = solve(tmp1*t(b)%*%b) %*% t(b) %*% tmp
          print(paste('Theta',k,q,'difference',norm(Theta[[k]][,q]-tmp)))
          Theta[[k]][,q] = tmp
        }
        Theta[[k]] = qr.orth(Theta[[k]])
        
      }        
      
      ######## sigma_gamma
      
      for(k in 1:cluster.num){
        for(q in 1:Q){
          tmp = rep(0,(MCC))
          for(t in 1:(MCC)){
            #print(t)
            z = Z_list[,t]
            Nz = rep(NA,cluster.num)
            for(i in 1:cluster.num){ Nz[i] = sum(z==i)}
            mem = which(z==1)
            for( i in 2:cluster.num){ mem = c(mem,which(z==i))}
            o2 = matrix(unlist(lapply(as.list(mem),'zero')),curve.num,)
            O2 = as.matrix(do.call('bdiag',rep(list(o2),Q)))
            O = O1 %*% O2
            gamma1 = ((solve(O)%*%gamma_list[,t])[(curve.num*(q-1)+1):(curve.num*q),1])[(if(k>1) sum(Nz[1:(k-1)])+1 else 1):sum(Nz[1:k])] 
            tmp[t] = t(as.matrix(gamma1)) %*% solve(P.func(phi,location[which(z==k),])) %*% as.matrix(gamma1)/Nz[k]			
          }
          print(paste('sigma_gamma',k,q,'difference',abs(mean(tmp)-sigma_gamma[k,q])))
          sigma_gamma[k,q] = mean(tmp)
        }
      }
      
      ######## phi        
      
      d.Obj <- function(phi){
        latent = rbind(Z_list,gamma_list)
        
        obj = function(x){
          z = x[1:curve.num] ; gamma = x[(curve.num+1):length(x)]
          Nz = rep(NA,cluster.num)
          for(i in 1:cluster.num){ Nz[i] = sum(z==i)}
          mem = which(z==1)
          for( i in 2:cluster.num){ mem = c(mem,which(z==i))}
          o2 = matrix(unlist(lapply(as.list(mem),'zero')),curve.num,)
          O2 = as.matrix(do.call('bdiag',rep(list(o2),Q)))
          O = O1 %*% O2
          
          tmp1 = 0
          tmp2 = 0
          for(k in 1:cluster.num){
            tmp1 = tmp1 + tr(solve(P.func(phi,location[which(z==k),]))%*%d.P.func(phi,location[which(z==k),]))
            for(q in 1:Q){
              gamma1 = ((solve(O)%*%gamma)[(curve.num*(q-1)+1):(curve.num*q),1])[(if(k>1) sum(Nz[1:(k-1)])+1 else 1):sum(Nz[1:k])] 
              
              tmp2 = tmp2 + tr(solve(P.func(phi,location[which(z==k),]))%*% d.P.func(phi,location[which(z==k),]) %*% solve(P.func(phi,location[which(z==k),])) %*% as.matrix(gamma1) %*% t(as.matrix(gamma1)))/sigma_gamma[k,q]
            }
          }
          return(Q * tmp1 - tmp2)            
        }
        return(mean(apply(latent,2,obj)))
      }
      
      
      print('optimization for phi')
      t3 = Sys.time()  
      phi1 <- tryCatch({  uniroot(d.Obj,c(0.01,5))$root},error=function(e){cat(conditionMessage(e));print('sb');return(phi)})
      print(Sys.time()-t3)
      a3 = abs(phi-phi1)/(abs(phi)+0.001)
      print(paste('phi','difference',a3))
      phi = phi1      
      
      ######### log likelihood
      loglike_list[loop] = loglike(y.data,O1, gamma_list,Z_list,lambda,Theta,phi,sigma_gamma,sigma,location,prob)    
      
      
      max1 = apply(con.prob,2,which.max)
      print(paste('z','diffrence',norm(max.con.prob-max1)))
      
      print(Sys.time()-a)
      a <- Sys.time()
      print(paste('sigma=',sigma,'phi=',phi,'loglike=',loglike_list[loop],'conv='))
      print(conv)
      print(sigma_gamma)
      
      #     if(norm(max.con.prob-max1)==0 & a1<0.005 &a2<0.005 & a3<0.005 &a4<0.005){conv = conv+1}
      if(a1<0.004){conv[1]=conv[1]+1}else{conv[1]=0}
      if(a2<0.001){conv[2]=conv[2]+1}else{conv[2]=0}
      if(a3<0.001){conv[3]=conv[3]+1}else{conv[3]=0}
      if(a4<0.001){conv[4]=conv[4]+1}else{conv[4]=0}
      if(a1<0.01){conv[5]=conv[5]+1}else{conv[5]=0}
      if(a2<0.005){conv[6]=conv[6]+1}else{conv[6]=0}
      if(a3<0.005){conv[7]=conv[7]+1}else{conv[7]=0}
      if(a4<0.005){conv[8]=conv[8]+1}else{conv[8]=0}
      
      #     if(max(a1,a2,a3,a4)<0.005) break
      
      max.con.prob <- max1
      print(round(adjustedRandIndex(max.con.prob,cluster_true),4))
      
      Sys.sleep(30)
      
      #write.csv(randIndex,paste(ite,'-',loop,'.csv',sep=''))
      
    }
    
    randIndex = round(adjustedRandIndex(max.con.prob,cluster_true),4)
    cluster=max.con.prob
    
    muhat = list()
    rms = 0
    for(k in 1:cluster.num){
      muhat[[k]]=as.vector(b%*%lambda[[k]])
      tmp = rep(NA,cluster.num)
      for(j in 1:cluster.num){
        tmp[j] = rmse(muhat[[k]],mu[[j]])
      }
      rms=rms+min(tmp)
    }
    RMSE1 = round(rms,4)
    
    muhat2 = list()
    rms = 0
    for(k in 1:cluster.num){
      muhat2[[k]] = apply(dd[max.con.prob==k,],2,mean)
      tmp = rep(NA,cluster.num)
      for(j in 1:cluster.num){
        tmp[j] = rmse(muhat2[[k]],mu[[j]])
      }
      rms=rms+min(tmp)
    }
    RMSE2 = round(rms,4)
    
  
    hfscm = list(conv=conv,cluster=cluster,phi=phi,sigma=sigma,sigma_gamma=sigma_gamma,
                 Theta=Theta,lambda=lambda,prob=prob,loglike_list=loglike_list)
    
    return(list(randIndex=randIndex,RMSE1=RMSE1,RMSE2=RMSE2,
                hfscm = hfscm,cluster_true=cluster_true,location=location))
    
    
  },error=function(e){cat(conditionMessage(e));return(NA)})
  
}   

stopCluster(cl)
#save.image("results\\Hfscm_results1.Rdata")


#################################################################
################      fit homoscedastic models      #############
#################################################################

##################################################################
##  Source files for homoscedastic case
library(funcy)

source('../../functions.R')
source('../homo/homo_iid_multi.R')
source('../homo/homo_iid_mrf.R')
source('../homo/homo_spatial_multi.R')
source('../homo/homo_spatial_mrf.R')


knn=5

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

##### simulation result records #####
randIndex = matrix(NA,Nsimu,7)
colnames(randIndex)=c('kmeans','Jiang','James','FMM','FMM-MRF','SFMM','SFMM-MRF')
RMSE = matrix(NA,Nsimu,7)
colnames(randIndex)=c('kmeans','Jiang','James','FMM','FMM-MRF','SFMM','SFMM-MRF')
phi_list = matrix(NA,Nsimu,2)
theta_list = matrix(NA,Nsimu,2)
sigma_gamma_list = matrix(NA,Nsimu,2)

P=6
aa <- Sys.time()
for(ite in 1:Nsimu){
  print(paste("Simulation No.",ite,sep=''))
  
  dd = data_list[[ite]]
  cluster_true = mrf[,ite]
  
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
  
  fit5 = homo_iid_mrf(data=t(dd),cluster.num=cluster.num,neighbor=y.neigh,P)
  
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
  
  ###############     mfscm (SFMM-MRF)    ###############
  
  
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
## Reproduce Figure 3 in the manuscript
##################################################################

Nsimu = 100

library(ggplot2)
library(fda)

reverse <- function(x){
  n=length(x)
  y=rep(NA,n)
  for(i in 1:n) y[i]=x[n+1-i]
  return(y)
}

index = rep(NA,Nsimu) ## to identify two clusters

day = 30
t = seq(0,1,length.out=day)
basis = create.bspline.basis(c(0,1),nbasis = P)
s = getbasismatrix(t,basis)
b = qr.Q(qr(s))

#-----------------------------------------------------------------
# Figure 3: 
# Top: Estimated Two Mean functions 
# Middle: Estimated Functional Principal Components in cluster 1
# Bottom: Estimated Functional Principal Components in cluster 2
#-----------------------------------------------------------------

#### mean 

cluster.num = 2
P = 4

lambda_list = NULL
lambda_list[[1]] = matrix(NA,Nsimu,P)
lambda_list[[2]] = matrix(NA,Nsimu,P)

for(i in 1:Nsimu){
  for(k in 1:cluster.num){
    tmp0 = rep(NA,cluster.num)
    tmp = results_HSFMM_MRF[[i]]$hmfscm$lambda[[k]]
    for(j in 1:cluster.num){
      tmp0[j]=norm(lambda_true[[j]]-tmp)
    }
    if(k==1) index[i]=which.min(tmp0)
    lambda_list[[which.min(tmp0)]][i,]=tmp
  }	
}

mu1_true = b%*%lambda_true[[1]]
mu1_hat = apply(b%*%t(lambda_list[[1]]),1,mean)
mu1_sd = apply(b%*%t(lambda_list[[1]]),1,sd)
mu1_upper = mu1_hat + 1.96*mu1_sd
mu1_lower = mu1_hat - 1.96*mu1_sd

mu1_list = data.frame(t = 1:day/day,mu1_true,mu1_hat,mu1_lower,mu1_upper)

plot_mu1 <-
  ggplot()+
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=14)
  ) +
  geom_line(data = mu1_list,
            aes(x=t,y=mu1_true),
            size=1) +
  geom_line(data = mu1_list,
            aes(x=t,y=mu1_hat),
            size=1,linetype="dashed") +
  geom_ribbon(data=mu1_list,aes(x=t,ymin=mu1_lower,ymax=mu1_upper),alpha=0.3) +
  ylab("Estimation of Mean Function 1") +
  xlab("Time") 

plot_mu1
ggsave(filename = 'Simulation_EstimatedMean1.pdf', 
       height = 3,
       width = 8)   

####################

mu2_true = b%*%lambda_true[[2]]
mu2_hat = apply(b%*%t(lambda_list[[2]]),1,mean)
mu2_sd = apply(b%*%t(lambda_list[[2]]),1,sd)
mu2_upper = mu2_hat + 1.96*mu2_sd
mu2_lower = mu2_hat - 1.96*mu2_sd

mu2_list = data.frame(t = 1:day/day,mu2_true,mu2_hat,mu2_lower,mu2_upper)

plot_mu2 <-
  ggplot()+
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=14)
  ) +
  geom_line(data = mu2_list,
            aes(x=t,y=mu2_true),
            size=1) +
  geom_line(data = mu2_list,
            aes(x=t,y=mu2_hat),
            size=1,linetype="dashed") +
  geom_ribbon(data=mu2_list,aes(x=t,ymin=mu2_lower,ymax=mu2_upper),alpha=0.3) +
  ylab("Estimation of Mean Function 2") +
  xlab("Time") 

plot_mu2
ggsave(filename = 'Simulation_EstimatedMean2.pdf', 
       height = 3,
       width = 8)   


##########  sigma_gamma  &  FPC  ###########

Theta_list = list()
Theta_list[[1]] = list()
Theta_list[[2]] = list()

Theta_list[[1]][[1]] = matrix(NA,Nsimu,4)
Theta_list[[1]][[2]] = matrix(NA,Nsimu,4)
Theta_list[[2]][[1]] = matrix(NA,Nsimu,4)
Theta_list[[2]][[2]] = matrix(NA,Nsimu,4)

Q = 2
sigma_gamma_1 = matrix(NA,Nsimu,Q)
sigma_gamma_2 = matrix(NA,Nsimu,Q)


for( i in 1:Nsimu){
  if(index[i]==1){
    tmp = results_HSFMM_MRF[[i]]$hmfscm$sigma_gamma[1,]
    if(tmp[1]>tmp[2]){
      sigma_gamma_1[i,] =  tmp
      Theta_list[[1]][[1]][i,] = results_HSFMM_MRF[[i]]$hmfscm$Theta[[1]][,1]
      Theta_list[[1]][[2]][i,] = results_HSFMM_MRF[[i]]$hmfscm$Theta[[1]][,2]
    }else{
      sigma_gamma_1[i,] =  reverse(tmp)
      Theta_list[[1]][[1]][i,] = results_HSFMM_MRF[[i]]$hmfscm$Theta[[1]][,2]
      Theta_list[[1]][[2]][i,] = results_HSFMM_MRF[[i]]$hmfscm$Theta[[1]][,1]
    }
    tmp = results_HSFMM_MRF[[i]]$hmfscm$sigma_gamma[2,]
    if(tmp[1]>tmp[2]){
      sigma_gamma_2[i,] =  tmp
      Theta_list[[2]][[1]][i,] = results_HSFMM_MRF[[i]]$hmfscm$Theta[[2]][,1]
      Theta_list[[2]][[2]][i,] = results_HSFMM_MRF[[i]]$hmfscm$Theta[[2]][,2]
    }else{
      sigma_gamma_2[i,] =  reverse(tmp)
      Theta_list[[2]][[1]][i,] = results_HSFMM_MRF[[i]]$hmfscm$Theta[[2]][,2]
      Theta_list[[2]][[2]][i,] = results_HSFMM_MRF[[i]]$hmfscm$Theta[[2]][,1]
    }
  }else{
    tmp = results_HSFMM_MRF[[i]]$hmfscm$sigma_gamma[1,]
    if(tmp[1]>tmp[2]){
      sigma_gamma_2[i,] =  tmp
      Theta_list[[2]][[1]][i,] = results_HSFMM_MRF[[i]]$hmfscm$Theta[[1]][,1]
      Theta_list[[2]][[2]][i,] = results_HSFMM_MRF[[i]]$hmfscm$Theta[[1]][,2]
    }else{
      sigma_gamma_2[i,] =  reverse(tmp)
      Theta_list[[2]][[1]][i,] = results_HSFMM_MRF[[i]]$hmfscm$Theta[[1]][,2]
      Theta_list[[2]][[2]][i,] = results_HSFMM_MRF[[i]]$hmfscm$Theta[[1]][,1]
    }
    tmp = results_HSFMM_MRF[[i]]$hmfscm$sigma_gamma[2,]
    if(tmp[1]>tmp[2]){
      sigma_gamma_1[i,] =  tmp
      Theta_list[[1]][[1]][i,] = results_HSFMM_MRF[[i]]$hmfscm$Theta[[2]][,1]
      Theta_list[[1]][[2]][i,] = results_HSFMM_MRF[[i]]$hmfscm$Theta[[2]][,2]
    }else{
      sigma_gamma_1[i,] =  reverse(tmp)
      Theta_list[[1]][[1]][i,] = results_HSFMM_MRF[[i]]$hmfscm$Theta[[2]][,2]
      Theta_list[[1]][[2]][i,] = results_HSFMM_MRF[[i]]$hmfscm$Theta[[2]][,1]
    }
  }
}  ## for identifiability


####  FPC1,1

FPC11_true = FPC[[1]][,1]
FPC11_hat = apply(b%*%t(Theta_list[[1]][[1]]),1,mean)
FPC11_sd = apply(b%*%t(Theta_list[[1]][[1]]),1,sd)
FPC11_upper = FPC11_hat + 1.96*FPC11_sd
FPC11_lower = FPC11_hat - 1.96*FPC11_sd

FPC11_list = data.frame(t = 1:day/day,FPC11_true,FPC11_hat,FPC11_lower,FPC11_upper)

plot_FPC11 <-
  ggplot()+
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=14)
  ) +
  geom_line(data = FPC11_list,
            aes(x=t,y=FPC11_true),
            size=1) +
  geom_line(data = FPC11_list,
            aes(x=t,y=FPC11_hat),
            size=1,linetype="dashed") +
  geom_ribbon(data = FPC11_list,
              aes(x=t,ymin=FPC11_lower,ymax=FPC11_upper),alpha=0.3) +
  ylab("Estimation of FPC1,1") +
  xlab("Time") 

plot_FPC11
ggsave(filename = 'Simulation_EstimatedFPC11.pdf', 
       height = 3,
       width = 8)   


####  FPC1,2

FPC12_true = FPC[[1]][,2]
FPC12_hat = apply(b%*%t(Theta_list[[1]][[2]]),1,mean)
FPC12_sd = apply(b%*%t(Theta_list[[1]][[2]]),1,sd)
FPC12_upper = FPC12_hat + 1.96*FPC12_sd
FPC12_lower = FPC12_hat - 1.96*FPC12_sd

FPC12_list = data.frame(t = 1:day/day,FPC12_true,FPC12_hat,FPC12_lower,FPC12_upper)

plot_FPC12 <-
  ggplot()+
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=14)
  ) +
  geom_line(data = FPC12_list,
            aes(x=t,y=FPC12_true),
            size=1) +
  geom_line(data = FPC12_list,
            aes(x=t,y=FPC12_hat),
            size=1,linetype="dashed") +
  geom_ribbon(data = FPC12_list,
              aes(x=t,ymin=FPC12_lower,ymax=FPC12_upper),alpha=0.3) +
  ylab("Estimation of FPC1,2") +
  xlab("Time") 

plot_FPC12
ggsave(filename = 'Simulation_EstimatedFPC12.pdf', 
       height = 3,
       width = 8) 


####  FPC2,1

FPC21_true = FPC[[2]][,1]
FPC21_hat = apply(b%*%t(Theta_list[[2]][[1]]),1,mean)
FPC21_sd = apply(b%*%t(Theta_list[[2]][[1]]),1,sd)
FPC21_upper = FPC21_hat + 1.96*FPC21_sd
FPC21_lower = FPC21_hat - 1.96*FPC21_sd

FPC21_list = data.frame(t = 1:day/day,FPC21_true,FPC21_hat,FPC21_lower,FPC21_upper)

plot_FPC21 <-
  ggplot()+
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=18)
  ) +
  geom_line(data = FPC21_list,
            aes(x=t,y=FPC21_true),
            size=1) +
  geom_line(data = FPC21_list,
            aes(x=t,y=FPC21_hat),
            size=1,linetype="dashed") +
  geom_ribbon(data = FPC21_list,
              aes(x=t,ymin=FPC21_lower,ymax=FPC21_upper),alpha=0.3) +
  ylab("Estimation of FPC2,1") +
  xlab("Time") 

plot_FPC21
ggsave(filename = 'Simulation_EstimatedFPC21.pdf', 
       height = 3,
       width = 8)   


####  FPC2,2

FPC22_true = FPC[[2]][,2]
FPC22_hat = apply(b%*%t(Theta_list[[2]][[2]]),1,mean)
FPC22_sd = apply(b%*%t(Theta_list[[2]][[2]]),1,sd)
FPC22_upper = FPC22_hat + 1.96*FPC22_sd
FPC22_lower = FPC22_hat - 1.96*FPC22_sd

FPC22_list = data.frame(t = 1:day/day,FPC22_true,FPC22_hat,FPC22_lower,FPC22_upper)

plot_FPC22 <-
  ggplot()+
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=18)
  ) +
  geom_line(data = FPC22_list,
            aes(x=t,y=FPC22_true),
            size=1) +
  geom_line(data = FPC22_list,
            aes(x=t,y=FPC22_hat),
            size=1,linetype="dashed") +
  geom_ribbon(data = FPC22_list,
              aes(x=t,ymin=FPC22_lower,ymax=FPC22_upper),alpha=0.3) +
  ylab("Estimation of FPC2,2") +
  xlab("Time") 

plot_FPC22
ggsave(filename = 'Simulation_EstimatedFPC22.pdf', 
       height = 3,
       width = 8)  



##################################################################
## Reproduce Table 3, Table 4 and in the manuscript
##################################################################

#-----------------------------------------------------------------
# Table 1: 
# Means and standard deviations of adjusted Rand index and RMSE 
#-----------------------------------------------------------------

RandIndex_list = cbind(randIndex[,1:5],matrix(NA,Nsimu,2))
colnames(RandIndex_list)[6:7] = c('HSFMM','HSFMM-MRF')
RMSE_list = cbind(RMSE[,1:5],matrix(NA,Nsimu,2))
colnames(RMSE_list) = colnames(RandIndex_list)

####  HSFMM
for(i in 1:Nsimu){
  RandIndex_list[i,6] = results_HSFMM[[i]]$randIndex
  RMSE_list[i,6] = results_HSFMM[[i]]$RMSE2
}
####  HSFMM-MRF
for(i in 1:Nsimu){
  RandIndex_list[i,7] = results_HSFMM_MRF[[i]]$randIndex
  RMSE_list[i,7] = results_HSFMM_MRF[[i]]$RMSE2
}

round(apply(RandIndex_list,2,mean),3)   ##  RandIndex
round(apply(RandIndex_list,2,sd),3)   

round(apply(RMSE_list,2,mean),3)  ##  RMSE
round(apply(RMSE_list,2,sd),3)


#-----------------------------------------------------------------
# Table 2: 
# Means and standard deviations of parameters using SFMM-MRF
#-----------------------------------------------------------------

para_list = matrix(NA,Nsimu,2)
for(i in 1:Nsimu){
  para_list[i,1] = results_HSFMM_MRF[[i]]$hmfscm$phi
  para_list[i,2] = results_HSFMM_MRF[[i]]$hmfscm$theta
}

para_list = cbind(para_list,sigma_gamma_1,sigma_gamma_2)
colnames(para_list) = c('phi','theta','sigma_gamma_11',
                        'sigma_gamma_21','sigma_gamma_12','sigma_gamma_22')

round(apply(para_list,2,mean),3) 
round(apply(para_list,2,sd),3)



