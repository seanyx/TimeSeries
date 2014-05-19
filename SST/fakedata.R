## This code is used for converting time into z (depth) dimension in sst and anom data file
library(RNetCDF) ## used for importing data from .nc file
source('~/Documents/TimeSeries/createNetCDF/CreateNetCDF.R', chdir = TRUE) #
source('~/Documents/TimeSeries/PaleoAnalyze/phaseshift.R')
# used for exporting .nc file

load('SSTXiao.RData')

ndim=dim(sst)
nt=ndim[4]
nlat=ndim[2]
nlon=ndim[1]

## fill the anom matrix with random time series
for (i in 1:nlat) {
	for (j in 1:nlon) {
		if(!is.na(anom[j,i,1,1])) anom[j,i,1,]=runif(nt,-1,1)
	}
}

##
## create the time series
set.seed(1987)
t=seq(1854,2014,length=nt)

y1=sin(2*pi*1/12*t)+0.1*runif(nt,-1,1) ## reference data

y2=y1+0.2*runif(nt,-1,1)  ## high corr (0.98) data at phase 0
y3=y1+0.8*runif(nt,-1,1) ## medium corr (0.83)
y4=y1*0.5+0.6*runif(nt,-1,1) ## low corr (0.69)

## regions for the reference signal
p1=c(176,-46)
p2=c(186,-34)

plant=function(p1,p2,y,lat,lon,anom,ns=0.1) {
	## plant time series y at the region defined by corner coordinates
	## p1 and p2, with noise level added
	indlon1=which(lon==p1[1])
	indlon2=which(lon==p2[1])
	indlat1=which(lat==p1[2])
	indlat2=which(lat==p2[2])
	indlon=c(indlon1,indlon2)
	indlat=c(indlat1,indlat2)
	std=sd(y)
	for (i in min(indlon):max(indlon)) {
		for(j in min(indlat):max(indlat)) {
			if(!is.na(anom[i,j,1,1])) anom[i,j,1,]=y+0.1*std*runif(length(y),-1,1)
		}
	}
	return(anom)
}

## plant the reference data
anom=plant(p1,p2,y1,lat,lon,anom,ns=0.1)

## plant the high correlation region at phase 0
## between North America and South America
p3=c(260,32)
p4=c(330,0)

p5=c(280,24)
p6=c(320,4)

p7=c(296,16)
p8=c(320,10)

anom=plant(p3,p4,y4,lat,lon,anom,ns=0.6)
anom=plant(p5,p6,y3,lat,lon,anom,ns=0.5)
anom=plant(p7,p8,y2,lat,lon,anom,ns=0.2)

## plant the high corr at phase 90
## near Madagascar
y2.90=phshift(y2,0.5)
y3.90=phshift(y3,0.5)
y4.90=phshift(y4,0.5)

p3.90=c(30,-30)
p4.90=c(52,-10)

p5.90=c(36,-20)
p6.90=c(40,-16)

anom=plant(p3.90,p4.90,y4.90,lat,lon,anom,ns=0.5)
anom=plant(p5.90,p6.90,y3.90,lat,lon,anom,ns=0.5)

## plant the high corr at phase 36
## to the southeast of Greenland
y2.36=phshift(y2,0.2)
y3.36=phshift(y3,0.2)
y4.36=phshift(y4,0.2)

p3.36=c(334,60) 
p4.36=c(350,70)

p5.36=c(340,64)
p6.36=c(350,70)

anom=plant(p3.36,p4.36,y4.36,lat,lon,anom,ns=0.5)
anom=plant(p5.36,p6.36,y3.36,lat,lon,anom,ns=0.5)

path='.'
filename='fake_anomaly_corr.nc'
time=seq(1854,2014,length=nt)
height=1
lat=lat
lon=lon
var.names='anom'
variables=list(anom)
null_value=-999
refs='fake ERSST v3b data'

createNetCDF2d(path,filename,time,height,lat,lon,var.names,variables,null_value,refs)
