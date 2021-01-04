##   Regionalization of China's PM2.5
regionalize_China <- function(data,station,cluster.num,neighbor,P){
	
  x = na.approx(data)
  curve.num = dim(x)[2]
  
  y.data <- list()
  for (i in 1:curve.num){
    y.data[[i]] <- na.approx(x[,i])
  }
    day = dim(x)[1]
  	
  t = seq(0,1,length.out=day)
  basis = create.bspline.basis(c(0,1),nbasis=P)
  s = getbasismatrix(t,basis)
  
  y.neigh = neighbor
  
  ###########   initial values  ###########
  
  theta <- 1 
  
  lambda00 <- list(solve(t(s)%*%s,t(s)%*%y.data[[8]]),  # Baoding
                   solve(t(s)%*%s,t(s)%*%y.data[[37]]), # Shenyang
                   solve(t(s)%*%s,t(s)%*%y.data[[310]]), # Xining
                   solve(t(s)%*%s,t(s)%*%y.data[[176]]), # Wuhan
                   solve(t(s)%*%s,t(s)%*%y.data[[14]]), # Taiyuan
                   solve(t(s)%*%s,t(s)%*%y.data[[195]]), # Guangzhou
                   solve(t(s)%*%s,t(s)%*%y.data[[73]]), # Shanghai
                   solve(t(s)%*%s,t(s)%*%y.data[[233]]), # Chengdu
                   solve(t(s)%*%s,t(s)%*%y.data[[263]]), # Kunming
                   solve(t(s)%*%s,t(s)%*%y.data[[134]]), # Qingdao
                   solve(t(s)%*%s,t(s)%*%y.data[[113]]) # Fuzhou                                    
  )
  lambda <- lambda00[1:cluster.num]
  
  benchcurve.mat00 = rbind(station[8,],station[37,],station[310,],station[176,],
                           station[14,],station[195,],station[73,],station[233,],
                           station[263,],station[134,],station[113,] )
  benchcurve.mat = benchcurve.mat00[1:cluster.num,]
  
  station_list = list()
  for(i in 1:curve.num){
    station_list[[i]] = station[i,]
  }
  cluster0 = unlist( lapply( station_list, function(x){ which.min( rdist(x,benchcurve.mat) ) } ) )
  
  con.prob <- z_ik(cluster0)
  
  mar.prob <- matrix(NA,cluster.num,curve.num)
  for(i in 1:curve.num){
    for(k in 1:cluster.num){
      tmp <- apply(as.matrix(con.prob[,y.neigh[[i]]]),1,function(t){exp(theta*sum(t))})
      mar.prob[k,i]<- tmp[k]/sum(tmp)
    }
  }
  
  sigma <- 1
  Gamma <- diag(1,ncol(s))         
  
  #####################
  a <- Sys.time()
  LOOP = 60
  loglike_list = rep(NA,LOOP)
  for(loop in 1:LOOP){
  #  print(loop)
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
          tmp = mapply(FUN="*",exp((tmp)-tmp[1]),mar.prob[,i])
          con.prob[,i] <- tmp/sum(tmp)
          }
      }
    
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
    theta <-  tryCatch({ uniroot(theta_func,c(0,5))$root },error=function(e){cat(conditionMessage(e));
      theta=theta})
    
    mar.prob <- matrix(NA,cluster.num,curve.num)
    for(i in 1:curve.num){
      for(k in 1:cluster.num){
        tmp <- apply(as.matrix(con.prob[,y.neigh[[i]]]),1,function(t){exp(theta*sum(t))})
        mar.prob[k,i]<- tmp[k]/sum(tmp)
      }
    }
    
    
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
    
  
  
    loglike_list[loop] = loglike(y.data,con.prob,lambda,s,Gamma,sigma,theta,y.neigh)
    
  #  print(Sys.time()-a)
  #  a <- Sys.time()
  #  print(c(loop,round(adjustedRandIndex(apply(con.prob,2,which.max),cluster_true),3)))
    # print(round(con.prob,2))
  }

  cluster0 = apply(con.prob,2,which.max)
  
  muhat = list()
  for(k in 1:cluster.num){
    muhat[[k]]=as.vector(s%*%lambda[[k]])
  }
  
  eigensystem = eigen(Gamma)
  O = eigensystem$vectors
  phi = s %*% t(O)
  eigenvalues = eigensystem$values
  

  BIC = -2 * loglike_list[loop] + log(curve.num*day) * cluster.num * P
  
  mfcm = list(data=data,cluster.num=cluster.num,cluster=cluster0,muhat=muhat,
              y.neigh=y.neigh,sigma=sigma,lambda=lambda,
              Gamma=Gamma,FPC = phi,eigenvalues = eigenvalues,
              theta=theta,loglike_list=loglike_list,BIC=BIC)

  return(mfcm)
}




## log likelihood ##
loglike = function(y.data,con.prob,lambda,s,Gamma,sigma,theta,y.neigh) {
  
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
    for(k in 1:cluster.num){ s3 = s3 + con.prob[k,i] * ( theta * sum(con.prob[k,y.neigh[[i]]]) - log( sum(apply(as.matrix(con.prob[,y.neigh[[i]]]),1,function(t){exp(theta*sum(t))})) )) }
  }     
  #	print(c(s1,s2,s3))
  return(s1 + s2 + s3)
  
}



