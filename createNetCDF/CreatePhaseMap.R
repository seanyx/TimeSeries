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
filename="PhaseMatrix.nc"
time=1
height=1
lat=lat
lon=lon

var.names=c('phaseanomfil','ampanomfil')
variables=list()
null_value=-999
refs='Extended Reconstructed Sea Surface Temperature (ERSST.v3b) http://www.ncdc.noaa.gov/ersst/'


## caculate the phase map and corr map
inilon=320
inilat=60

inipoint=c(which(lon==inilon),which(lat==inilat))

phaseanomfil=array(NA,dim=c(180,89,1,1))
# phaseanom=phaseanomfil
# phasesst=phaseanomfil
ampanomfil=phaseanomfil

for (i in 1:length(lon)) {
	for (j in 1:length(lat)) {
		if(is.na(anomfil[i,j,1,][1])) {#corrsst[i,j,1,1]=-999; 
			next}
		if(i==inipoint[1] & j==inipoint[2]) {test=list(0,1); next}	
		test=fitPhase(anomfil[inipoint[1],inipoint[2],1,],anomfil[i,j,1,],N=180)
		phaseanomfil[i,j,1,1]=test[[1]]
		ampanomfil[i,j,1,1]=test[[2]]
	}
	
}

variables[[1]]=phaseanomfil
variables[[2]]=ampanomfil

null_value=rep(-999,length(var.names))
filename=paste('lon',inilon,'_lat',inilat,'_',filename,sep='')
createNetCDF2d(path,filename,time,height,lat,lon,var.names,variables,null_value=null_value,refs=refs)
