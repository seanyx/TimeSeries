## This code is used for calculate correlation map for a range of phase shift
## based on filtered anomaly SST data (file: anomfil2to10.RData)
library(RColorBrewer)
source('~/Documents/TimeSeries/createNetCDF/phaseshiftest.R', chdir = TRUE) #
# used to calculate maximum correlation and its phase shift

load('anomfil2to10.RData')
load('SSTXiao.RData')

anomfil=anomfil2to10

ndim=dim(anomfil)
nt=ndim[4]
nlat=ndim[2]
nlon=ndim[1]

inilon=320
inilat=0
ini_ind=c(which(lon==inilon),which(lat==inilat))

nshift=101 ## number of steps in phase shift
shift=seq(0,1,length=nshift) ## shift range 0 to 2*pi
corr=array(NA,dim=c(180,89,nshift)) ## set up the array accepting corr value

print(paste('loop starts at',Sys.time()))
for (p in 1:nshift) {
	for (i in 1:nlon) {
		for (j in 1:nlat) {
			if(is.na(anomfil[i,j,1,1])) next ## NA data
			temp=phshift(anomfil[i,j,1,],shift[p])
			corr[i,j,p]=cor(temp,anomfil[ini_ind[1],
					ini_ind[2],1,],method="pearson")
		}
	}
}
print(paste('loop ends at',Sys.time())) # used ~2min
filename=paste('corr',inilon,'&',inilat,'_anomfil2to10_N',nshift,'.RData',sep='')
save(corr,file=filename)

corrfil=corr

cap=0.5
corrfil[which(corrfil<=cap & corrfil>=-cap,arr.ind=TRUE)]=NA

## building color map using RColorBrewer
nc=11 ## number of colors in the color map
# display.brewer.pal(n=nc,name='RdYlGn') ## preview the color map
# display.brewer.all() ## preview the color map
cols=brewer.pal(n=nc,name='RdBu')
cols=rev(cols)

## building continents outline
outline=corr[,,1]
outline[which(!is.na(corr[,,1]),arr.ind=TRUE)]=1
outline[which(is.na(corr[,,1]),arr.ind=TRUE)]=0 
outline[which(outline[,]==1,arr.ind=TRUE)]=NA

for (i in 1:nshift) {
	png(file=paste('figures/',filename,'_',i,'.png',sep=''),
	    bg='transparent',width=1200,height=800)
	image(lon,lat,corrfil[,,i],col=cols,zlim=c(-1,1),main=paste('corr >',cap))
	image(lon,lat,outline,col='black',add=TRUE) ## add continents for reference
	points(inilon,inilat,pch=22,col='darkgreen',lwd=2,bg='yellow')
	dev.off()
}
