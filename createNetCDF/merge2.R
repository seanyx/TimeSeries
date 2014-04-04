library(RNetCDF)
setwd("~/Dropbox/netcdf_test/")

test=create.nc('SST_Xiao1.nc')
#file need to be merged
flist=list.files('~/Desktop/sst_netcdf/netcdf')
Nf=length(flist)

#define dimensions
dim.def.nc(test,"time",Nf)
dim.def.nc(test,"zlev",1)
dim.def.nc(test,"lat",89)
dim.def.nc(test,"lon",180)

#define variables
var.def.nc(test,"time","NC_FLOAT",0)
var.def.nc(test,"zlev","NC_FLOAT",1)
var.def.nc(test,"lat","NC_FLOAT",2)
var.def.nc(test,"lon","NC_FLOAT",3)
var.def.nc(test,"sst","NC_SHORT",c(3,2,1,0))
var.def.nc(test,"anom","NC_SHORT",c(3,2,1,0))
var.def.nc(test,"anomfil","NC_SHORT",c(3,2,1,0))
var.def.nc(test,"magnitude","NC_DOUBLE",c(3,2))
var.def.nc(test,"err","NC_SHORT",c(3,2,1,0))

#add fillValues attribute to variables
att.put.nc(test,'sst','_FillValue',"NC_SHORT",-999)
att.put.nc(test,'anom','_FillValue',"NC_SHORT",-999)
att.put.nc(test,'anomfil','_FillValue',"NC_SHORT",-999)
att.put.nc(test,'err','_FillValue',"NC_SHORT",-999)
att.put.nc(test,'magnitude','_FillValue',"NC_DOUBLE",-999)
att.put.nc(test,'NC_GLOBAL','Data source',"NC_CHAR",'Extended Reconstructed Sea Surface Temperature (ERSST.v3b) http://www.ncdc.noaa.gov/ersst/')



#variable value
time=1:Nf
sst=array(0,dim=c(180,89,1,Nf))
anom=sst
anomfil=sst
err=sst
magnitude=array(0,dim=c(180,89))

for ( i in 1:Nf ) {
	d1=open.nc(paste('~/Desktop/sst_netcdf/netcdf/',flist[i],sep=''))
	sst[,,,i]=var.get.nc(d1,'sst',collapse=F)[,,,1]
	anom[,,,i]=var.get.nc(d1,'anom',collapse=F)[,,,1]
	err[,,,i]=var.get.nc(d1,'err',collapse=F)[,,,1]
	close.nc(d1)
}


#add values to variable
d1=open.nc(paste('~/Desktop/sst_netcdf/netcdf/',flist[1],sep=''))

## calculate the anomfil
fil=1/rep(12,12)

fil.anom=function(anom,fil) {
	## remove linear trend and apply moving average filter using fil
	N=length(anom)
	t=1:N
	anomres=lm(anom~t)$residuals
	anomfil=filter(anomres,fil,circular=T)
	return(anomfil)
}
na.seq=rep(NA,1920)

for (i in 1:length(var.get.nc(d1,'lon'))) {
	for (j in 1:length(var.get.nc(d1,'lat'))) {
		if(is.na(anom[i,j,1,1])) {anomfil[i,j,1,]=na.seq; magnitude[i,j]=NA; next}
		#if((var(anom[i,j,1,])<=0.1)) {anomfil[i,j,1,]=na.seq; next}
		if(any(sst[i,j,1,]==-180)) {
			anomfil[i,j,1,]=na.seq
			magnitude[i,j]=NA
			next
			}
		magnitude[i,j]=max(sst[i,j,1,],na.rm=T)-min(sst[i,j,1,],na.rm=T)
		anomfil[i,j,1,]=fil.anom(anom[i,j,1,],fil)
	}
}

var.put.nc(test,'time',time)
var.put.nc(test,'zlev',var.get.nc(d1,'zlev'))
var.put.nc(test,'lat',var.get.nc(d1,'lat'))
var.put.nc(test,'lon',var.get.nc(d1,'lon'))
var.put.nc(test,'sst',sst)
var.put.nc(test,'magnitude',magnitude)
var.put.nc(test,'anom',anom)
var.put.nc(test,'anomfil',anomfil)
var.put.nc(test,'err',err)



close.nc(test)
close.nc(d1)
