## This code is used for converting time into z (depth) dimension in sst and anom data file
library(RNetCDF) ## used for importing data from .nc file
source('~/Documents/TimeSeries/createNetCDF/CreateNetCDF.R', chdir = TRUE) #
# used for exporting .nc file

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

sst3d=array(dim=c(nlon,nlat,nt,1))
temp=anom[,,1,]
sst3d[,,,1]=temp

path='.'
filename='anom3d.nc'
time=1
height=seq(1854,2014,length=nt)
lat=lat
lon=lon
var.names='anom'
variables=list(sst3d)
null_value=-999
refs='ERSST v3b data using depth as time'

createNetCDF2d(path,filename,time,height,lat,lon,var.names,variables,null_value,refs)
