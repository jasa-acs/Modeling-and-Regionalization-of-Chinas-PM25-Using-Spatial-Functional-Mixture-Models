##################################################################
##  Spatial-Functional Mixture Models: Regionalization of China ##
##  Outputs: Figures 1, 4, 5, 6, and 9 of the manuscript.       ##
##################################################################
rm(list = ls())

##################################################################
##  Load associated required packages
library(fda)
library(funcy)
library(fields)
library(mclust)
library(MASS)
library(zoo)
library(ggmap) 
library(maptools)
library(plyr)
library(scales)


##################################################################
##  Source files for the algorithm functions 
source('../functions.R')
source('regionalize_China.R')

  
##################################################################
##  Read original data

station0 = read.csv('../../data/city_location_China.csv',sep=',',fileEncoding  = 'GBK')
d0 = read.csv('../../data/China_PM25_daily.csv',sep=",",fileEncoding  = 'GBK')

##################################################################
## Reproduce Figure 1 in the manuscript
##################################################################
#-----------------------------------------------------------------
# Figure 1: 
# Location maps of all China's cities with times series of 4 cities
#-----------------------------------------------------------------
#### Location maps of all China's cities
#-----------------------------------------------
china_map = readShapePoly('../../data/China_boundary_map/bou2_4p.shp')
x = china_map@data
xs = data.frame(x,id=seq(0:924)-1)
china_map1 = fortify(china_map)
china_map_data0 = join(china_map1,xs,type='full',by='id')

p = ggplot(china_map_data0,aes(x=long,y=lat))
p = p +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none"
  )
p1 = p + geom_polygon(aes(group=group),fill="white",color="grey60")

N0 = dim(station0)[1]
color = rep("grey53",N0)
color[1] = "red1" #### Beijing
color[2] = "green4" #### Tianjin
color[73] = "blue1" #### Shanghai
color[195] = "violet" #### Guangzhou

size = rep(1,N0)
size[c(1,2,73,195)] = 3
alpha = rep(1,N0)
alpha[c(1,2,73,195)] = 1

p2 = p1 + geom_point(data=station0, aes(lon,lat),color = color,size=size,alpha=alpha)

png(filename = 'station_outline.png',width = 700,height = 600)
print(p2)
dev.off()

####  Time series of 4 cities
#-----------------------------------------------
date = seq(as.Date("2015-01-01"), as.Date("2016-12-31"), by="1 day")

city = 1 ## Beijing; 2: Tianjin, 73: Shanghai, 195: Guangzhou
c = data.frame(date=date,y=d0[,city+3])

g1 = ggplot(na.omit(c)) + 
  aes(x=date,y=y) + 
  geom_line(size=2,
            colour=color[city])+
  labs(title = "Daily PM2.5 Concentration in Beijing ",y="",x='') +
  scale_y_continuous(breaks = c(50,200,350),limits=c(0,400)) +
  scale_x_date( date_labels = "%Y/%m", breaks=date_breaks("6 month")) +
  theme(axis.ticks = element_blank(),
        title = element_text(size=40),plot.title=element_text(hjust=0.5),
        axis.text.y=element_text(size=30,hjust=0.5),
        axis.text.x=element_text(size=30,hjust=0.5),
        legend.position = "none"
  )

png(filename = paste('data1516_overview_Beijing.png',sep=""),width = 1200,height = 500)
print(g1)
dev.off()


##################################################################
############            Regionalization         ##################
##################################################################

##################################################################
##  data cleaning

c1 = c(71,279:285,311,323:338) ## remove Xinjiang & Tibet

station = station0[-c1,]
location = station[,c('lon','lat')]

index1 = as.numeric(rownames(station))
data = na.approx(d0[,-(1:3)])
data = data[,-c1]

curve.num = dim(data)[2]

colnames(data) = 1:curve.num
log.data = log(data)

##################################################################
##  neighbor  

load('../../data/National338_neighbor.Rdata')  
## can also use:  source('GetNeighbors.R') 

y.neigh1 = list()       ##   Removes neighbors in Xinjiang & Tibet
for(i in 1:curve.num){
  tmp = as.vector(y.neigh[[index1[i]]])
  y.neigh1[[i]] = apply( matrix(tmp[!(tmp %in% c1)],,1),1,function(xx){which(index1==xx)} )
}

##################################################################
##  regionalization   

cluster.num = 9

fit.China = regionalize_China(log.data, location, cluster.num,neighbor = y.neigh1, P = 20)

cluster = fit.China$cluster


##################################################################
## Reproduce Figure 4 and Figure 5 in the manuscript
##################################################################

#-----------------------------------------------------------------
# Figure 4: 
# Regionalization map of China 
#-----------------------------------------------------------------

china_map_data = china_map_data0[china_map_data0$lat>18.2,]  ##  remove the South China Sea
xxx = c(250,307) ## remove some outliers

## add region label 
region_label = c("North China Plain",
                 "Northeast China Plain", 
                 "Northwest", 
                 "Middle Yangtze River Plain",
                 "Guanzhong Plain",
                 "Pearl River Delta",
                 "Jianghuai Plain",
                 "Sichuan Basin",
                 "Yungui Plateau")
region_func <- function(i){ region_label[as.numeric(i)] }

region_label2 = c("North China Plain",
                  "Northeast China Plain", 
                  "Guanzhong Plain",
                  "Middle Yangtze River Plain",
                  "Jianghuai Plain",
                  "Northwest", 
                  "Sichuan Basin",
                  "Pearl River Delta",
                  "Yungui Plateau")

out_clustermap = station[-xxx,]
out_clustermap$region = factor( region_func(cluster[-xxx]), levels = region_label2 )

plot_clustermap <- 
  ggplot(china_map_data,aes(x=long,y=lat)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.spacing.y = unit(0.01, 'cm'),
        legend.position = "bottom", 
        legend.title = element_blank()
  ) +
  geom_polygon(aes(group=group),fill="white",color="grey60") +
  geom_point(data=out_clustermap,
             aes(lon,lat,colour = region, shape = region),
             size=2 )  +
  scale_colour_manual(
    labels = region_label2,
    values = c("black", "green", "blue", 
               "green","blue", "red", 
               "black","red", "blue")) +
  scale_shape_manual(
    labels = region_label2,
    values = c(16,2,4,
               4,2,16,
               2,4,16))  +
  guides(colour = guide_legend(nrow = 3,ncol=3, byrow = T),
          shape = guide_legend(nrow = 3,ncol=3, byrow = T) ) 


plot_clustermap
ggsave(filename = paste('RealData_ChinaClusterMap',cluster.num,'.pdf',sep=""),
       height = 6,
       width = 8)   


#-----------------------------------------------------------------
# Figure 5: 
# Estimated mean functions of China 
#-----------------------------------------------------------------

rawdata = fit.China$data
day = dim(rawdata)[1]

curve.raw = data.frame(t=rep(1:day,curve.num),
                        date = date,
                        y=as.vector(rawdata),
                        region = rep(cluster,each = day) )

##   estimated cluster mean & points
mu.hat = fit.China$muhat
region = as.factor(cluster)

mean.clustered = matrix(NA,day,cluster.num)
mean.region = matrix(NA,day,cluster.num)
for(i in 1:cluster.num){
  mean.clustered[,i] = mu.hat[[i]] 
  mean.region[,i] = rep(i,day)
}
mean.est = data.frame(t = rep(1:day,cluster.num),
                       date = date,
                       y = as.vector(mean.clustered),
                       region = as.factor(as.vector(mean.region)) )

meanpoint.index = seq(50,650,100)
meanpoint.est = data.frame(date = date[meanpoint.index],
                            y = as.vector(mean.clustered[meanpoint.index,]),
                            region = as.factor(as.vector(mean.region[meanpoint.index,])) )

curve.raw$region = factor(region_func(curve.raw$region), levels = region_label2)
mean.est$region = factor(region_func(mean.est$region), levels = region_label2)
meanpoint.est$region = factor(region_func(meanpoint.est$region), levels =region_label2)


plot_mean <-
  ggplot()+
  theme(legend.position = "none", 
        legend.title = element_blank()
  ) +
  geom_line(aes(x=date,y=y),data=curve.raw,colour="grey60",alpha=0.2,size=1) +
  geom_line(data=mean.est,
            aes(x=date,y=y,colour=region),
            size=1) +
  scale_x_date( date_labels = "%Y/%m", breaks=date_breaks("6 month")) +
  ylab("Estimated Mean Functions") +
  xlab("Time") +
  facet_wrap(~region,
             nrow=4,ncol=3) +
  scale_colour_manual(
    values = c("black", "green", "blue", 
               "green","blue", "red", 
               "black","red", "blue")) +
  geom_point(data=meanpoint.est,
             aes(x=date,y=y,colour=region,shape=region),
             size=2) +
  scale_shape_manual(
    values = c(16,2,4,
               4,2,16,
               2,4,16)) 

plot_mean
ggsave(filename = paste('RealData_ChinaEstimatedMean_',cluster.num,'.png',sep=""), 
       height = 8,
       width = 8)   


#-----------------------------------------------------------------
# Figure 6: 
# Estimated leading FPCs of China
#-----------------------------------------------------------------

FPC = fit.China$FPC

comp.name = c("1st principal component","2nd principal component",
              "3rd principal component","4th principal component",
              "5th principal component","6th principal component")

## check the FVE 
eigenvalues = fit.China$eigenvalues
pn = sum(fve(eigenvalues)<0.8)

comp.est = data.frame(date = date,
                      y = as.vector(FPC[,1:pn]),
                      comp = rep(comp.name[1:pn],each=day) )

ggplot()+
  geom_line(aes(x=date,y=y),data=comp.est)+
  facet_wrap(~comp,nrow=3,ncol=2  )+
  scale_x_date( date_labels = "%Y/%m", breaks=date_breaks("4 month")) +
  ylab("Estimated Principal Components")+
  xlab("Time (Day)")
# which might be slightly different with Figure 6, as the pattern is more important

ggsave(filename = paste('RealData_ChinaEstimatedPC_',cluster.num,'.pdf',sep=""), 
       height = 8,
       width = 8)   

#-----------------------------------------------------------------
# Figure 9: 
# Estimated global and regional mean functions
#-----------------------------------------------------------------

mu0 = apply(mean.clustered,1,mean)
mu0.est = data.frame(   date = date, y = as.vector(mu0) )

plot_globalmean <-
  ggplot()+
  theme(legend.position = "bottom", 
        legend.title = element_blank()
  ) +
  geom_line(data = mu0.est,
            aes(x=date,y=y),
            size = 1) +
  ylab("Glaobal Mean") +
  xlab("") +
  ylim(c(3,4.5)) +
  scale_x_date( date_labels = "%Y/%m", breaks=date_breaks("6 month")) 

plot_globalmean  ## global mean function
ggsave(filename = 'RealData_ChinaDemean_globalmean.png', 
       height = 2,
       width = 6)   


demean.clustered = mean.clustered - mu0
demean.est = data.frame(t = rep(1:day,cluster.num),
                       date = date,
                       y = as.vector(demean.clustered),
                       region = as.factor(as.vector(mean.region)))

demeanpoint.est = data.frame(date = date[meanpoint.index],
                            y = as.vector(demean.clustered[meanpoint.index,]),
                            region = as.factor(as.vector(mean.region[meanpoint.index,])) )

demean.est$region = factor(region_func(demean.est$region), levels = region_label2)
demeanpoint.est$region = factor(region_func(demeanpoint.est$region), levels =region_label2)


plot_demean <-
  ggplot(data = demean.est)+
  theme(legend.position = "none", 
        legend.title = element_blank()
  ) +
  geom_hline(aes(yintercept=0), colour = "brown", linetype="dashed" ) + 
  geom_line(data = demean.est,
            aes(x=date,y=y,colour=region),
            size=1) +
  scale_x_date( date_labels = "%Y/%m", breaks=date_breaks("6 month")) +
  ylab("Estimated Mean Functions") +
  xlab("Time") +
  facet_wrap(~region,
             nrow=4,ncol=3) +
  scale_colour_manual(
    values = c("black", "green", "blue", 
               "green","blue", "red", 
               "black","red", "blue")) +
  geom_point(data=demeanpoint.est,
             aes(x=date,y=y,colour=region,shape=region),
             size=2) +
  scale_shape_manual(
    values = c(16,2,4,
               4,2,16,
               2,4,16)) 

plot_demean  ## regiona mean functions
ggsave(filename = 'RealData_ChinaDemean_regionalmean.png',
       height = 8,
       width = 8)  


##################################################################
##  Additional codes for other regionalization methods

# Kmeans
# fit1 = kmeans(t(na.approx(log.data)),cluster.num,iter.max=100)

# Jiang
#fit2 = funcit(data=na.approx(log.data),method='fscm',knn=5,k=cluster.num,
#              location=station,verbose=T)

# James
# fit3 = funcit(data=na.approx(log.data),method='fitfclust',k=cluster.num)

# FMM
# fit4 = homo_iid_multi(data=log.data,cluster.num)







