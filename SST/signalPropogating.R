## This code is used for calculate correlation map for a range of phase shift
library(RNetCDF) ## used for importing data from .nc file
library(RColorBrewer)
source('~/Documents/TimeSeries/createNetCDF/CreateNetCDF.R', chdir = TRUE) #
# used for exporting .nc file
source('~/Documents/TimeSeries/createNetCDF/phaseshiftest.R', chdir = TRUE) #
# used to calculate maximum correlation and its phase shift

## importing sst, anom, and anomfil data set
# data=open.nc("~/Dropbox/netcdf_test/SST_Xiao.nc")
# lat=var.get.nc(data,'lat')
# lon=var.get.nc(data,'lon')
# sst=var.get.nc(data,'sst',collapse=F)
# anom=var.get.nc(data,'anom',collapse=F)
# anomfil=var.get.nc(data,'anomfil',collapse=F)
# close.nc(data)
# 
# save(lat,lon,sst,anom,anomfil,file='SSTXiao.RData')
load('SSTXiao.RData')

ndim=dim(sst)
nt=ndim[4]
nlat=ndim[2]
nlon=ndim[1]

inilon=320
inilat=60
ini_ind=c(which(lon==inilon),which(lat==inilat))

# nshift=10 ## number of steps in phase shift
# shift=seq(0,2,length=nshift) ## shift range 0 to 2*pi
# corr=array(NA,dim=c(180,89,nshift)) ## set up the array accepting corr value
# 
# for (p in 1:nshift) {
#         for (i in 1:nlon) {
#                 for (j in 1:nlat) {
#                         if(is.na(anomfil[i,j,1,1])) next ## NA data
#                         if(any(sst[i,j,1,]==-180)) next ## missing data
#                         temp=phshift(anomfil[i,j,1,],shift[p])
#                         corr[i,j,p]=cor(temp,anomfil[ini_ind[1],
#                                  ini_ind[2],1,],method="pearson")
#                 }
#         }
# }
# save(corr,file='corr_anomfil_N10.RData')
load('corr_anomfil_N10.RData')
corrfil=corr

cap=0.01
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

image(lon,lat,corr[,,1],col=cols,zlim=c(-1,1))
image(lon,lat,outline,col='black',add=TRUE) ## add continents for reference
