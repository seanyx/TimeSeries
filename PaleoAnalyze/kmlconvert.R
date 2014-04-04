## import kml to R
# install.packages('rgdal')
setwd('/Users/yangxiao/Dropbox/Research/PaleoSearch backup/alldata/paleoclimatology dataset')
library(rgdal)
library(stringr)
paleo=readOGR('Paleoclimatology.kml','Paleoclimatology')

# str(paleo)
# Formal class 'SpatialPointsDataFrame' [package "sp"] with 5 slots
  # ..@ data       :'data.frame':	13042 obs. of  2 variables:
  # .. ..$ Name       : Factor w/ 10617 levels "?RSON","\"LD\" LAKE (SEPTILES)",..: 6051 6441 8294 3855 1574 5979 6073 2064 624 6180 ...
  # .. ..$ Description: Factor w/ 12546 levels "<table><tr><td>NETWORK: BOREHOLE</td></tr><tr><td>INVESTIGATOR: ADAMS, J.</td></tr><tr><td>LATITUDE: 31.45</td></tr><tr><td>LON"| __truncated__,..: 2908 6314 12481 5481 2951 1937 2338 9584 10946 2187 ...
  # ..@ coords.nrs : num(0) 
  # ..@ coords     : num [1:13042, 1:3] -125 -50 99.8 -33.5 -76.3 ...
  # .. ..- attr(*, "dimnames")=List of 2
  # .. .. ..$ : NULL
  # .. .. ..$ : chr [1:3] "coords.x1" "coords.x2" "coords.x3"
  # ..@ bbox       : num [1:3, 1:2] -179.8 -90 0 179.5 88.9 ...
  # .. ..- attr(*, "dimnames")=List of 2
  # .. .. ..$ : chr [1:3] "coords.x1" "coords.x2" "coords.x3"
  # .. .. ..$ : chr [1:2] "min" "max"
  # ..@ proj4string:Formal class 'CRS' [package "sp"] with 1 slots
  # .. .. ..@ projargs: chr "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"


# write csv directly is not good. Description has mixed information
#writeOGR(paleo,'paleoNOAA.csv','Paleoclimatology',driver='csv')



t1=paleo@data$Description
t1=as.character(t1)
t2=substring(t1,16,nchar(t1)-18)
t3=str_split(t2,'</td></tr><tr><td>')
t4=unlist(t3)
t5=matrix(t4,ncol=7,byrow=T)
t6=t5[,-6]
t6[,1]=substring(t6[,1],10,nchar(t6[,1]))
t6[,2]=substring(t6[,2],15,nchar(t6[,2]))
t6[,3]=substring(t6[,3],11,nchar(t6[,3]))
t6[,4]=substring(t6[,4],12,nchar(t6[,4]))
t6[,5]=substring(t6[,5],13,nchar(t6[,5]))
t6[,6]=substring(t6[,6],10,nchar(t6[,6])-17)
t7=as.data.frame(t6,stringsAsFactors=F)
t7[,4]=as.numeric(t7[,4])
t7[,3]=as.numeric(t7[,3])
names(t7)=c('NETWORK','INVESTIGATOR','Latitude','Longitude','STUDY NAME','WEBPAGE')
NAME=as.character(paleo@data$Name)
tf=cbind(NAME,t7)
tf[,1]=as.character(tf[,1])

#write.csv(tf,file='paleo.csv',row.names=F)

nw=levels(as.factor(tf$NETWORK))
nwcolor=hcl(h=seq(10,350,length=18),c=90,l=50,alpha=0.7)
nwcol=vector(length=length(tf[,1]))
data=vector('list',length=length(nw))
for (i in 1:length(nw)) {
	nwcol[which(tf$NETWORK==nw[i])]=nwcolor[i]
	data[[i]]=tf[which(tf$NETWORK==nw[i]),]
	#names(data[[i]])=nw[i]
}
names(data)=nw

####plot
library(GEOmap)
data(worldmap)
require('geomapdata')
proj=setPROJ(type=1,0,350)

mdata=tf
mdata[which(mdata$Longitude<0),'Longitude']=mdata[which(mdata$Longitude<0),'Longitude']+360
quartz()
plotGEOmap(worldmap , border=NA , add=FALSE, xaxs='i',lwd=3)
pointsGEOmapXY(mdata$Latitude,mdata$Longitude,proj,col=nwcol,cex=1,pch=20)
xy=locator(1)
legend(xy$x,xy$y,nw,col=nwcolor,pch=20,pt.cex=2,bty='n',horiz=F)

####individual plot
for (i in 1:length(data)) {
	data[[i]][which(data[[i]]$Longitude<0),'Longitude']=data[[i]][which(data[[i]]$Longitude<0),'Longitude']+360
	pdf(file=paste(nw[i],'.pdf',sep=''),width=16,height=8)
	plotGEOmap(worldmap , border=NA , add=FALSE, xaxs='i',lwd=3,main=nw[i])
	pointsGEOmapXY(data[[i]]$Latitude,data[[i]]$Longitude,proj,col=nwcolor[i],cex=1.5,pch=20)
	dev.off()
}
