tr <- function(x){if(is.matrix(x)) sum(diag(x))}  ## trace

rmse <- function(x,y){sqrt(sum((x-y)^2)/sum(y^2))}  ## RMSE

P.func <- function(phi,location){Exp.cov(location,theta=phi)} ## exponential covariance
        
d.P.func <- function(phi,location){
  exp(log(rdist(location,location))-rdist(location,location)/phi)/(phi^2) 
} ## exponential covariance, derivative 

qr.orth <- function(x){
	y = qr.Q(qr(x))
	for(i in 1:dim(y)[2]){
		if(y[1,i]<0){ y[,i] = -y[,i] }
	}
	return(y)
} ## orthogonalize with identification

fve <- function(x){
  d=length(x)
  FVE = x/sum(x)
  for(i in 2:d){FVE[i]= FVE[i-1] + FVE[i]}
  return(FVE)
}  ## FVE

z_ik <- function(z){
  L<-length(z)
  tmp=matrix(0,max(z),L)
  indx <- 1:L
  tmp[matrix(c(z,indx),L)] <-1
  return(tmp)
}  ## identity 

