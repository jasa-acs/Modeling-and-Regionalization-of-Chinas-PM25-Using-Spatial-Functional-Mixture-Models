############  Simulate markov random fields

library(PottsUtils)

cluster.num = 2
theta_true = 0.5

station = read.csv('../../../data/city_location_China.csv',sep=',',fileEncoding  = 'GBK')  
station = station[ 107<station$lon & station$lon<125 & 28<station$lat & station$lat<43, ]

location = station[,3:4]
curve.num = dim(location)[1]

### markov random fields

knn=5
dist=as.matrix(dist(location))
dist.ord=apply(dist,2,order)
dist.order=t(dist.ord)
nb=dist.order[,2:(knn+1)]

edges = NULL
for(i in 1:curve.num){
  for(j in (i+1):curve.num){
    if(sum(nb[i,]==j)>0){edges = rbind(edges,c(i,j))}
  }
}

Nsimu = 100

set.seed(123)

mrf = SW(n=Nsimu, nvertex=curve.num, ncolor=cluster.num, edges, beta=theta_true)

#write.csv(mrf,'mrf_data.csv',row.names = F)


