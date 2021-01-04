##    fcm (FMM)

homo_iid_multi <- function(data,cluster.num,P){
	x = na.approx(data)
	curve.num = dim(x)[2]
    y.data <- list()
    for (i in 1:curve.num){
      y.data[[i]] <- na.approx(x[,i])
      #  print(i)
    }
    day = dim(x)[1]
	
	t = seq(1,day,length.out=day)
    basis = create.bspline.basis(c(1,day),nbasis=P)
    s = getbasismatrix(t,basis)
    
    ###########   initial values  ###########

    prob = rep(1/cluster.num,cluster.num)
  
    fit1 = kmeans(t(data),cluster.num,iter.max=100)
    con.prob <- z_ik(fit1$cluster)

    lambda <- list()
    for(i in 1:cluster.num){
      lambda[[i]] = solve(t(s)%*%s,t(s)%*%fit1$centers[i,])
    }
      
    sigma <- 1
    Gamma <- diag(1,ncol(s))  
    
    #####################
    a <- Sys.time()
    loglike_list = rep(NA,60)
    for(loop in 1:60){
    
      tmp.mat <- solve(sigma*diag(1,nrow(s))+s%*%Gamma%*%t(s))
    
    for(i in 1:curve.num){
        tmp = lapply(lambda,function(tt){-0.5*t(y.data[[i]]-s%*%tt)%*%
                        tmp.mat%*%
                        (y.data[[i]]-s%*%tt)})  
        tmp = unlist(tmp)
        if((max(tmp)-min(tmp))>500){
        	con.prob[,i]=rep(0,cluster.num)
        	con.prob[which.max(tmp),i]=1
        }else{ 
        tmp = mapply(FUN="*",exp((tmp)-tmp[1]),prob)
        con.prob[,i] <- tmp/sum(tmp)
        }
    }
    prob=apply(con.prob,1,mean)
    
    v.con.g <- solve(solve(Gamma)+t(s)%*%s/sigma)
    tmp.mat <- solve(sigma*solve(Gamma)+t(s)%*%s)%*%t(s)
    for(k in 1:cluster.num){
      tmp <- matrix(0,ncol(s),1)
      for(i in 1:curve.num){
        tmp <- tmp + con.prob[k,i]*t(s)%*%(y.data[[i]]-s%*%tmp.mat%*%(y.data[[i]]-s%*%lambda[[k]]))
      }
      lambda[[k]] <- solve(sum(con.prob[k,])*t(s)%*%s)%*%tmp
    }
    
    tmp <- matrix(0,ncol(s),ncol(s))
    for(k in 1:cluster.num){
      for(i in 1:curve.num){
        e.con.g <- tmp.mat%*%(y.data[[i]]-s%*%lambda[[k]])
        tmp <- tmp + con.prob[k,i]*(v.con.g+e.con.g%*%t(e.con.g))
      }
    }
    Gamma <- tmp/curve.num
    
    tmp <- 0
    v.con.g <- solve(solve(Gamma)+t(s)%*%s/sigma)
    tmp.mat <- solve(sigma*solve(Gamma)+t(s)%*%s)%*%t(s)
    for(k in 1:cluster.num){
      for(i in 1:curve.num){
        e.con.g <- tmp.mat%*%(y.data[[i]]-s%*%lambda[[k]])
        tmp.vec <- y.data[[i]]-s%*%lambda[[k]]-s%*%e.con.g
        tmp <- tmp + con.prob[k,i]*(t(tmp.vec)%*%tmp.vec+sum(diag(s%*%v.con.g%*%t(s))))
      }
    }
    sigma <- as.numeric(tmp/(curve.num*day))
    
    loglike_list[loop] = loglike_fcm(y.data,con.prob,lambda,s,Gamma,sigma,prob) 
    
#    print(Sys.time()-a)
#    a <- Sys.time()
#    print(c(loop,round(adjustedRandIndex(apply(con.prob,2,which.max),cluster_true),3)))
    # print(round(con.prob,2))
  }

  cluster0 = apply(con.prob,2,which.max)
  
  muhat = list()
  for(k in 1:cluster.num){
    muhat[[k]]=as.vector(s%*%lambda[[k]])
  }

  BIC = -2 * loglike_list[loop] + log(curve.num*day) * cluster.num * P
  
  fcm = list(data=data,cluster.num=cluster.num,cluster=cluster0,muhat=muhat,BIC=BIC,
            sigma=sigma,Gamma=Gamma,lambda=lambda,prob=prob,loglike_list=loglike_list)

  return(fcm)

}




homo_iid_multi2 <- function(data,cluster.num,fit4,P){
  x = na.approx(data)
  curve.num = dim(x)[2]
  y.data <- list()
  for (i in 1:curve.num){
    y.data[[i]] <- na.approx(x[,i])
    #  print(i)
  }
  day = dim(x)[1]
  
  Q = 2
  
  t = seq(0,1,length.out=day)
  basis = create.bspline.basis(c(0,1),nbasis=P)
  s = getbasismatrix(t,basis)
  
  
  
  
  b = qr.Q(qr(s))
  
  #########################################	
  ###########   initial values  ###########
  #########################################
  
  prob = rep(1/cluster.num,cluster.num)
  
  
  lambda <- list()
  for(i in 1:cluster.num){
    lambda[[i]] = solve(t(b)%*%b,t(b)%*%fit4$muhat[[i]])
  }
  
  
  con.prob <- z_ik(fit4$cluster)
  
  #      max.con.prob <- apply(con.prob,2,which.max)
  
  ##### random Para : lambda, sigma, theta
  
  # lambda <- list()
  # for(i in 1:cluster.num){
  # lambda[[i]] = solve(t(b)%*%b,t(b)%*%y.data[[which(max.con.prob==i)[1]]])
  # }
  # lambda <- list()
  # for(i in 1:cluster.num){
  #   lambda[[i]] = solve(t(b)%*%b,t(b)%*%fit4$muhat[[i]])
  # }
  
  sigma <- 1
  
  ############
  Theta = rbind(diag(Q),matrix(0,P-Q,Q))
  
  naive <- c(2,1.5,1,0.5,0.3,0.2,0.1,0.05,0.02,0.01)
  sigma_gamma = naive[1:Q]
  #      phi = 0.5
  #####################
  a <- Sys.time()
  loglike_list = rep(NA,60)
  for(loop in 1:60){
    #    print(paste('iteration',loop))
    ############################
    ########## E-step ##########
    ############################
    
    s = b %*% Theta
    
    ############ z_ik
    
    tmp.mat <- solve(sigma*diag(1,nrow(s))+s%*%diag(sigma_gamma)%*%t(s))    
    for(i in 1:curve.num){
      tmp = lapply(lambda,function(tt){-0.5*t(y.data[[i]]-b%*%tt)%*%
          tmp.mat%*%
          (y.data[[i]]-b%*%tt)})  
      tmp = unlist(tmp)
      if((max(tmp)-min(tmp))>500){
        con.prob[,i]=rep(0,cluster.num)
        con.prob[which.max(tmp),i]=1
      }else{ 
        tmp = mapply(FUN="*",exp((tmp)-tmp[1]),prob)
        con.prob[,i] <- tmp/sum(tmp)
      }
    }
    max1 = apply(con.prob,2,which.max)
    #    print(paste('z','diffrence',norm(max.con.prob-max1)))
    max.con.prob <- max1
    
    ######### gamma | z_ik
    
    s = b %*% Theta
    tmp.mat <- solve(sigma*solve(diag(sigma_gamma))+t(s)%*%s)%*%t(s)   
    
    e.con.g = NULL
    for(i in 1:curve.num){
      e.con.g[[i]] = matrix(NA,cluster.num,Q)
      for(k in 1:cluster.num){
        e.con.g[[i]][k,] <- tmp.mat%*%(y.data[[i]]-b%*%lambda[[k]])
      }
    }
    v.con.g <- solve(solve(diag(sigma_gamma))+t(s)%*%s/sigma)
    
    ##################################
    #############  M-step ############
    ##################################
    
    ########    prob
    
    prob1 = apply(con.prob,1,mean)   
    a4 = abs(prob[1]- prob1[1])/(prob[1]+0.001)
    #    print(paste('pi','difference',a4))
    prob = prob1
    
    ########   lambda_k
    
    for(k in 1:cluster.num){
      tmp <- matrix(0,ncol(b),1)
      for(i in 1:curve.num){
        tmp <- tmp + con.prob[k,i]*t(b)%*%(y.data[[i]]-s%*%e.con.g[[i]][k,])
      }
      tmp <- solve(sum(con.prob[k,])*t(b)%*%b)%*%tmp
      a1 = sqrt(sum((tmp-lambda[[k]])^2)/P)/sqrt(sum(lambda[[k]]^2)/P+0.001)
      #      print(paste('lambda',k,'difference',a1))
      lambda[[k]] <- tmp
    }
    
    ########   sigma
    
    tmp <- 0
    for(k in 1:cluster.num){
      for(i in 1:curve.num){
        tmp.vec <- y.data[[i]]-b%*%lambda[[k]]-s%*%e.con.g[[i]][k,]
        tmp <- tmp + con.prob[k,i]*(t(tmp.vec)%*%tmp.vec+tr(s%*%v.con.g%*%t(s)))
      }
    }
    sigma1 <- as.numeric(tmp/(curve.num*day))
    a2 = abs(sigma-sigma1)/(abs(sigma)+0.001)
    #    print(paste('sigma','difference',a2))
    sigma = sigma1
    
    #########  Theta
    for(q in 1:Q){
      tmp1 = rep(0,cluster.num)
      tmp = rep(0,day)
      for(k in 1:cluster.num){
        for(i in 1:curve.num){
          e2.con.g = v.con.g + e.con.g[[i]][k,] %*% t(e.con.g[[i]][k,])
          tmp1[k] = tmp1[k] + con.prob[k,i]*(e2.con.g[q,q])
          tmp = tmp + con.prob[k,i] * (y.data[[i]] - b %*% lambda[[k]]) * e.con.g[[i]][k,][q] - con.prob[k,i] * b %*% Theta[,-q] %*% e2.con.g[q,-q]
        }
      }
      tmp = solve(sum(tmp1)*t(b)%*%b) %*% t(b) %*% tmp   
      #      print(paste('Theta',q,'difference',norm(Theta[,q]-tmp)))
      Theta[,q] = tmp       
    }
    Theta = qr.orth(Theta)
    
    ######## sigma_gamma
    
    for(q in 1:Q){
      tmp =0
      for(k in 1:cluster.num){
        for(i in 1:curve.num){
          tmp <- tmp + con.prob[k,i]*(v.con.g[q,q]+(e.con.g[[i]][k,][q])^2 )
        }
      }
      tmp <- tmp/curve.num
      #      print(paste('sigma_gamma',q,'difference',abs(tmp-sigma_gamma[q])))
      sigma_gamma[q] = tmp
    }
    
    loglike_list[loop] = loglike_fcm2(b,y.data,e.con.g,v.con.g,con.prob,lambda,Theta,sigma_gamma,sigma,prob)  
    
    #     print(Sys.time()-a)
    #     a <- Sys.time()
    #     print(paste('sigma=',sigma,'phi=',phi,'loglike=',loglike_list[loop]))
    #    print(sigma_gamma)
    
    #     print(round(adjustedRandIndex(max.con.prob,cluster_true),4))
    
  }
  
  cluster0 = apply(con.prob,2,which.max)
  
  muhat = list()
  for(k in 1:cluster.num){
    muhat[[k]]=as.vector(b%*%lambda[[k]])
  }
  
  BIC = -2 * loglike_list[loop] + log(curve.num*day) * cluster.num * P
  
  fcm = list(data=data,cluster.num=cluster.num,cluster=cluster0,muhat=muhat,BIC=BIC,
             sigma=sigma,sigma_gamma=sigma_gamma,Theta=Theta,lambda=lambda,prob=prob,loglike_list=loglike_list)
  
  return(fcm)
  
}


## log likelihood (iid,multinomial)##
loglike_fcm = function(y.data,con.prob,lambda,s,Gamma,sigma,prob) {
  
  curve.num = length(y.data)
  day = length(y.data[[1]])
  cluster.num = length(lambda)
  
  z_ik <- function(z){
    L<-length(z)
    tmp=matrix(0,max(z),L)
    indx <- 1:L
    tmp[matrix(c(z,indx),L)] <-1
    return(tmp)
  }
  
  tr <- function(x){if(is.matrix(x)) sum(diag(x))}	
  
  tmp.mat <- solve(sigma*solve(Gamma)+t(s)%*%s)%*%t(s)
  v.gamma <- solve(solve(Gamma)+t(s)%*%s/sigma)
  
  s1 = 0; s2=0; s3 = 0
  
  #################    s1    ###############
  
  for(i in 1:curve.num){
    for(k in 1:cluster.num){
      e.gamma <- tmp.mat%*%(y.data[[i]]-s%*%lambda[[k]])
      s1 = s1 - 0.5 * con.prob[k,i] * (day*log(sigma) + (t( y.data[[i]] - s %*% (lambda[[k]] + e.gamma) ) %*% ( y.data[[i]] - s %*% (lambda[[k]] + e.gamma) ) + tr(s %*% v.gamma %*% t(s))) / (sigma) )
    }
  }
  
  #################    s2    ###############	
  for(i in 1:curve.num){
    for( k in 1:cluster.num){
      e.gamma <- tmp.mat%*%(y.data[[i]]-s%*%lambda[[k]])
      s2 = s2 - 0.5 * con.prob[k,i] * ( log(det(Gamma))+tr( solve(Gamma) %*% (v.gamma+e.gamma%*%t(e.gamma)) ) )
    }        
  }
  #################    s3    ###############	
  for(i in 1: curve.num){
    for(k in 1:cluster.num){ s3 = s3 + con.prob[k,i] * log(prob[k])}
  }    
  #	print(c(s1,s2,s3))
  return(s1 + s2 + s3)
  
}

## log likelihood (iid,multinomial)##
loglike_fcm2 = function(b,y.data,e.con.g,v.con.g,con.prob,lambda,Theta,sigma_gamma,sigma,prob) {
  
  curve.num = length(y.data)
  day = length(y.data[[1]])
  cluster.num = length(lambda)
  Q = length(sigma_gamma)	
  s = b %*% Theta
  
  zero = function(x){ze = rep(0,curve.num);ze[x]=1;return(ze)}
  
  
  s1 = 0; s2=0; s3 = 0
  
  #################    s1    ###############
  
  for(i in 1:curve.num){
    for(k in 1:cluster.num){
      s1 = s1 - 0.5 * con.prob[k,i] * (day*log(sigma) + (t( y.data[[i]] - b %*% (lambda[[k]] + Theta %*% e.con.g[[i]][k,]) ) %*% ( y.data[[i]] -  b %*% (lambda[[k]] + Theta %*% e.con.g[[i]][k,]) ) + tr(s %*% v.con.g %*% t(s))) / (sigma) )
    }
  }
  
  #################    s2    ###############	
  
  Gamma = diag(sigma_gamma)        
  for(i in 1:curve.num){
    for( k in 1:cluster.num){
      s2 = s2 - 0.5 * con.prob[k,i] * ( log(det(Gamma))+tr( solve(Gamma) %*% (v.con.g+e.con.g[[i]][k,]%*%t(e.con.g[[i]][k,])) ) )
    }        
  }
  
  #################    s3    ###############	
  for(i in 1: curve.num){
    for(k in 1:cluster.num){ s3 = s3 + con.prob[k,i] * log(prob[k]) }
  }    
  
  return(s1 + s2 + s3)
  
}







