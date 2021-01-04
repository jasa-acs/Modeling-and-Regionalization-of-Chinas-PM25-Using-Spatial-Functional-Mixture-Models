######  For national 338neighbor selection

####   elevation data
library(rgdal)

ele0 = readGDAL('../../data/China_elevation_1km.tif')

coor0 = SpatialPoints(ele0)
#ele0@grid@cells.dim = 7397 5669
lon.mat = matrix(coor0@coords[,1],ele0@grid@cells.dim[1],)
lat.mat = matrix(coor0@coords[,2],ele0@grid@cells.dim[1],)

ele.mat = matrix(as.vector(ele0@data[,1]),ele0@grid@cells.dim[1],)

grid.x = dim(lon.mat)[1]
grid.y = dim(lat.mat)[2]

####  moubtain block function

M.height = 1000

lonlat.index <- function(location){
  bb = (lon.mat - as.numeric(location)[1])^2 + (lat.mat - as.numeric(location)[2])^2
  return( which( bb == bb[which.min(bb)], arr.ind = T ) )
}


############################################

station0 = read.csv('../../data/city_location_China.csv',sep=",",fileEncoding  = 'GBK' )
location = station0[,3:4]
curve.num = dim(location)[1]

neigh.num <- 10
colnames(location)=c('lon','lat')
y.location <- location

#distance on earth
earth.dist <- function (location1, location2) 
{
  long1 <- location1[,1]
  lat1 <- location1[,2]
  long2 <- location2[1]
  lat2 <- location2[2]
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.145
  d <- R * c
  return(d)
}

y.dist <- list()
for(i in 1:curve.num){
  y.dist[[i]] <- apply(y.location,1,function(tt){earth.dist(y.location[i,],tt)})
  names(y.dist[[i]]) <- 1:curve.num
}


y.neigh.origin <- list()
for(i in 1:curve.num){
  y.neigh.origin[[i]] <- as.numeric(names(sort(y.dist[[i]],decreasing=F)[2:(1+neigh.num)]))
  tmp <- sapply(y.neigh.origin[[i]],FUN=function(tt){if(y.dist[[i]][tt]>500) return(NA)
    else return(tt)})
  if(sum(is.na(tmp))==neigh.num)  y.neigh.origin[[i]] <-  y.neigh.origin[[i]][1]
  else y.neigh.origin[[i]] <- as.matrix(na.omit(tmp))
}

#save(y.neigh.origin,file='National338_neighbor_origin.Rdata')

###########   mountain block effect

diff.alt <- function(loc1,neigh){
  loc2 <- y.location[neigh,]
  index1 = lonlat.index(loc1)
  index2 = lonlat.index(loc2)
  
  ele.grid = ele.mat[ ( min(index1[1],index2[1]) : max(index1[1],index2[1]) ) , ( min(index1[2],index2[2]) : max(index1[2],index2[2]) ) ]
  
  if( (max(ele.grid, na.rm=T) - min(ele.grid, na.rm=T)) > M.height) return(NA)
  else return(neigh)
}


y.neigh <- list()
for(i in 1:curve.num){
  print(c(i,y.neigh.origin[[i]]))
  tmp <- na.omit(apply(as.matrix(y.neigh.origin[[i]]),1,function(x){diff.alt(y.location[i,],x)}))
  if(length(tmp)==0) y.neigh[[i]] <-  y.neigh.origin[[i]][1]
  else y.neigh[[i]] <- as.matrix(tmp)
  print(c(i,y.neigh[[i]]))
  #  print(station$cityname[i])
  #  print(station$cityname[y.neigh[[i]]])
  Sys.sleep(2)
} 

save(y.neigh,file='../../Regionalization_data/National338_neighbor.Rdata')


