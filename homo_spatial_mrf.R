##   mfscm (SFMM-MRF)

homo_spatial_mrf2 <- function(data,cluster.num,location,neighbor,fit6,P){

  	x = na.approx(data)
  	curve.num = dim(x)[2]
      
    y.data = as.list(as.data.frame(na.approx(x)))
    Y = as.matrix(unlist(y.data))
    day = dim(x)[1]
    
    Q = 2

    t = seq(0,1,length.out=day)
    basis = create.bspline.basis(c(0,1),nbasis=P)
    s = getbasismatrix(t,basis)
    

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

    y.neigh = neighbor
    
    ###########   initial values  ###########

    theta = 1      
    con.prob <- z_ik(fit6$cluster)

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
      lambda[[i]] = solve(t(b)%*%b,t(b)%*%fit6$muhat[[i]])
    }
          
    sigma <- 1

    ############
    Theta = fit6$Theta

    ############  spatial Para ##############
    
    naive <- c(2,1.5,1,0.5,0.3,0.2,0.1,0.05,0.02,0.01)
    sigma_gamma = fit6$sigma_gamma
    phi = 0.5
    
    #####################
    a <- Sys.time()
    loglike_list = rep(NA,60)
    for(loop in 1:60){
    # print(paste('iteration',loop))
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
        tmp = mapply(FUN="*",exp((tmp)-tmp[1]),mar.prob[,i])
        con.prob[,i] <- tmp/sum(tmp)
        }
    }

    max1 = apply(con.prob,2,which.max)
#    print(paste('z','diffrence',norm(max.con.prob-max1)))
    max.con.prob <- max1
    
    ######### gamma | z_ik
    
    tmp.mat <- solve(sigma*solve(diag(sigma_gamma))+t(s)%*%s)%*%t(s)   

    e.con.g = NULL
    for(i in 1:curve.num){
      e.con.g[[i]] = matrix(NA,cluster.num,Q)
      for(k in 1:cluster.num){
        e.con.g[[i]][k,] <- tmp.mat%*%(y.data[[i]]-b%*%lambda[[k]])
      }
    }
    v.con.g <- solve(solve(diag(sigma_gamma))+t(s)%*%s/sigma)

    ######### E2(gamma) | z
    
    z = max.con.prob

    e.Alpha = do.call('rbind',lapply(as.list(z),function(x){unlist(lambda[[x]])}))
       
    e.S = as.matrix(do.call('bdiag',rep(list(b %*% Theta),curve.num)))
                
    e.Gamma1 = as.matrix(do.call('bdiag',mapply(FUN='*',as.list(sigma_gamma),rep(list(P.func(location,phi = phi)),Q),SIMPLIFY=F)))
        
    e.Gamma = O1 %*% e.Gamma1 %*% t(O1)
    
    #### var(Y)
    v_Y = e.S %*%e.Gamma%*%t(e.S) + diag(sigma,curve.num*day)  ### computation
    #solve(v_Y)
    ### Woodbury inverse
    inv_v_Y = diag(1/sigma,curve.num * day)-1/(sigma^2) * e.S %*% solve((solve(e.Gamma)+t(e.S) %*% e.S/sigma))%*%t(e.S)
          
    #### expectation
    E.gamma = e.Gamma %*% t(e.S) %*% inv_v_Y %*% (Y - B %*% e.Alpha)
    #### variance
    V.gamma = e.Gamma - e.Gamma %*% t(e.S) %*% inv_v_Y %*%  e.S %*% e.Gamma
    
    E2.gamma = E.gamma %*% t(E.gamma) + V.gamma
    E2.gamma1 = solve(O1) %*% E2.gamma %*% t(solve(O1))

    ##################################
    #############  M-step ############
    ##################################
    
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
 #   print(paste('sigma','difference',a2))
    sigma = sigma1

     #########  Theta
     for(q in 1:Q){
         tmp1 = rep(0,cluster.num)
         tmp = matrix(0,cluster.num,day)
         for(k in 1:cluster.num){
             for(i in 1:curve.num){
                 e2.con.g = v.con.g + e.con.g[[i]][k,] %*% t(e.con.g[[i]][k,])
                 tmp1[k] = tmp1[k] + con.prob[k,i]*((e.con.g[[i]][k,][q])^2+v.con.g[q,q])
                 tmp[k,]  = tmp[k,]   + con.prob[k,i] * ((y.data[[i]]-b%*%lambda[[k]])* e.con.g[[i]][k,][q]
                                                         - b %*% Theta[,-q] %*% e2.con.g[q,-q] )
             }
         }
         tmp = solve(sum(tmp1)*t(b)%*%b) %*% t(b) %*% apply(tmp,2,sum)      
#         print(paste('Theta',q,'difference',norm(Theta[,q]-tmp)))
         Theta[,q] = tmp       
     }
     Theta = qr.orth(Theta)
    
    

    ######## sigma_gamma
   
    for(q in 1:Q){
      tmp = tr( solve(P.func(phi,location)) %*% E2.gamma1[((q-1)*curve.num+1):(q*curve.num),((q-1)*curve.num+1):(q*curve.num)] )/curve.num	
#      print(paste('sigma_gamma',q,'difference',abs(tmp-sigma_gamma[q])))
      sigma_gamma[q] = tmp
    }

    ######## phi        
    d.Obj <- function(phi){          	        
      tmp1 = tr(solve(P.func(phi,location))%*%d.P.func(phi,location))
      tmp2 = 0
      for(q in 1:Q){
        tmp2 = tmp2 + tr(solve(P.func(phi,location))%*% d.P.func(phi,location) %*% solve(P.func(phi,location)) %*%  E2.gamma1[(curve.num*(q-1)+1):(curve.num*q),(curve.num*(q-1)+1):(curve.num*q)])/sigma_gamma[q]
      }
      return(Q * tmp1 - tmp2)            
    }

#    print('optimization for phi')
#    t3 = Sys.time()  
    phi1 <- tryCatch({  uniroot(d.Obj,c(0.01,20))$root},error=function(e){cat(conditionMessage(e));return(phi)})
#    print(Sys.time()-t3)
    a3 = abs(phi-phi1)/(abs(phi)+0.001)
#    print(paste('phi','difference',a3))
    phi = phi1   

    # d.Obj <- function(phi){          	        
      # tmp1 = tr(solve(P.func(phi,location))%*%d.P.func(phi,location))
      # tmp2 = 0
      # for(q in 1:Q){
        # tmp2 = tmp2 + tr(solve(P.func(phi,location))%*% d.P.func(phi,location) %*% solve(P.func(phi,location)) %*%  gamma[,q] %*% t(gamma[,q]))/sigma_gamma[q]
      # }
      # return(Q * tmp1 - tmp2)            
    # }

######### log likelihood
     loglike_list[loop] = loglike_mfscm(b,y.data,e.con.g,v.con.g,con.prob,E2.gamma,lambda,Theta,phi,sigma_gamma,sigma,location,theta,y.neigh)  

 #    print(Sys.time()-a)
     a <- Sys.time()
 #    print(paste('sigma=',sigma,'phi=',phi,'theta=',theta,'loglike=',loglike_list[loop]))
#    print(sigma_gamma)

#     if(max(a1,a2,a3,a4)<0.005) conv=conv+1

#     print(round(adjustedRandIndex(max.con.prob,cluster_true),4))

#     Sys.sleep(5)
}
	

 	cluster0 = apply(con.prob,2,which.max)

  muhat = list()
  for(k in 1:cluster.num){
    muhat[[k]]=as.vector(b%*%lambda[[k]])
  }

  BIC = -2 * loglike_list[loop] + log(curve.num*day) * cluster.num * P
  
    
  mfscm = list(data=data,cluster.num=cluster.num,cluster=cluster0,muhat=muhat,y.neigh=y.neigh,theta=theta,
  sigma=sigma,sigma_gamma=sigma_gamma,phi=phi,lambda=lambda,loglike_list=loglike_list,BIC=BIC)


  return(mfscm)


}

## log likelihood (spatial,mrf)##
loglike_mfscm = function(b,y.data,e.con.g,v.con.g,con.prob,E2.gamma,lambda,Theta,phi,sigma_gamma,sigma,station,theta,y.neigh) {
  
  curve.num = length(y.data)
  day = length(y.data[[1]])
  cluster.num = length(lambda)
  Q = length(sigma_gamma)	
  s = b %*% Theta
  
  
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
  
  s1 = 0; s3 = 0
  
  #################    s1    ###############
  
  for(i in 1:curve.num){
    for(k in 1:cluster.num){
      s1 = s1 - 0.5 * con.prob[k,i] * (day*log(sigma) + (t( y.data[[i]] - b %*% (lambda[[k]] + Theta %*% e.con.g[[i]][k,]) ) %*% ( y.data[[i]] -  b %*% (lambda[[k]] + Theta %*% e.con.g[[i]][k,]) ) + tr(s %*% v.con.g %*% t(s))) / (sigma) )
    }
  }
  
  #################    s2    ###############	
  
  
  Gamma1 = as.matrix(do.call('bdiag',mapply(FUN='*',as.list(sigma_gamma),rep(list(P.func(station,phi = phi)),Q),SIMPLIFY=F)))
  
  Gamma = O1 %*% Gamma1 %*% t(O1)          
  
  s2 =  -0.5*( log(det(Gamma))+tr(solve(Gamma) %*% E2.gamma) )
  
  #################    s3    ###############	
  for(i in 1: curve.num){
    for(k in 1:cluster.num){ s3 = s3 + con.prob[k,i] * ( theta * sum(con.prob[k,y.neigh[[i]]]) - log( sum(apply(as.matrix(con.prob[,y.neigh[[i]]]),1,function(t){exp(theta*sum(t))})) )) }
  }  
  
  return(s1 + s2 + s3)
  
  
}





