library(RNetCDF)
source('~/Dropbox/netcdf_test/CreateNetCDF.R', chdir = TRUE)
source('~/Documents/TimeSeries/createNetCDF/phaseshiftest.R', chdir = TRUE)

olddir=getwd()
setwd("~/Dropbox/netcdf_test/")
data=open.nc("SST_Xiao.nc")
lat=var.get.nc(data,'lat')
lon=var.get.nc(data,'lon')
sst=var.get.nc(data,'sst',collapse=F)
anom=var.get.nc(data,'anom',collapse=F)
anomfil=var.get.nc(data,'anomfil',collapse=F)
close.nc(data)
setwd(olddir)

path="."
filename="PhaseMatrix.nc"
time=1
height=1
lat=lat
lon=lon

##### create phase map for the following location
##### only need to change the following 4 lines
da=anomfil
var.names=c('phaseanomfil','ampanomfil')
## caculate the phase map and corr map
inilon=80
inilat=0

variables=list()
null_value=-999
refs='Extended Reconstructed Sea Surface Temperature (ERSST.v3b) http://www.ncdc.noaa.gov/ersst/'



inipoint=c(which(lon==inilon),which(lat==inilat))

phasemap=array(NA,dim=c(180,89,1,1))
# phaseanom=phaseanomfil
# phasesst=phaseanomfil
ampmap=phasemap

for (i in 1:length(lon)) {
	for (j in 1:length(lat)) {
		if(is.na(da[i,j,1,][1])) {#corrsst[i,j,1,1]=-999; 
			next}
		if(i==inipoint[1] & j==inipoint[2]) {test=list(0,1); next}	
		test=fitPhase(da[inipoint[1],inipoint[2],1,],da[i,j,1,],N=18)
		phasemap[i,j,1,1]=test[[1]]
		ampmap[i,j,1,1]=test[[2]]
	}
	
}

variables[[1]]=phasemap
variables[[2]]=ampmap

null_value=rep(-999,length(var.names))
filename=paste('lon',inilon,'_lat',inilat,'_',filename,sep='')
createNetCDF2d(path,filename,time,height,lat,lon,var.names,variables,null_value=null_value,refs=refs)
