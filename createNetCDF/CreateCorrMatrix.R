library(RNetCDF)
source('~/Dropbox/netcdf_test/CreateNetCDF.R', chdir = TRUE)
source('~/Dropbox/netcdf_test/phaseshiftest.R', chdir = TRUE)


setwd("~/Dropbox/netcdf_test/")
data=open.nc("SST_Xiao.nc")
lat=var.get.nc(data,'lat')
lon=var.get.nc(data,'lon')
sst=var.get.nc(data,'sst',collapse=F)
anom=var.get.nc(data,'anom',collapse=F)
anomfil=var.get.nc(data,'anomfil',collapse=F)
close.nc(data)


path="~/Dropbox/netcdf_test/"
filename="corrMap.nc"
time=1
height=1
lat=lat
lon=lon

var.names=c('Corrsst','Corranom','Corranomfil')
variables=list()
null_value=-999
refs='Extended Reconstructed Sea Surface Temperature (ERSST.v3b) http://www.ncdc.noaa.gov/ersst/'


## caculate the phase map and corr map
inilon=320
inilat=60

inipoint=c(which(lon==inilon),which(lat==inilat))

corrsst=array(NA,dim=c(180,89,1,1))
corranom=corrsst
corranomfil=corrsst

phaseshift=45

for (i in 1:length(lon)) {
	for (j in 1:length(lat)) {
		if(is.na(anom[i,j,1,1])) next
		if(any(sst[i,j,1,]==-180)) next
		tempsst=phshift(sst[i,j,1,],phaseshift/180)
		tempanom=phshift(anom[i,j,1,],phaseshift/180)
		tempanomfil=phshift(anomfil[i,j,1,],phaseshift/180)
		corrsst[i,j,1,1]=cor( sst[inipoint[1],inipoint[2],1,], tempsst , method="pearson")
		corranom[i,j,1,1]=cor( anom[inipoint[1],inipoint[2],1,], tempanom , method="pearson")
		corranomfil[i,j,1,1]=cor( anomfil[inipoint[1],inipoint[2],1,], tempanomfil , method="pearson")
	}
	
}

variables[[1]]=corrsst
variables[[2]]=corranom
variables[[3]]=corranomfil
null_value=c(-999,-999,-999)

filename=paste('lon',inilon,'_lat',inilat,'_phaseshift',phaseshift,'_',filename,sep='')
createNetCDF2d(path,filename,time,height,lat,lon,var.names,variables,null_value=null_value,refs=refs)
