library(RNetCDF)
setwd("~/Dropbox/netcdf_test/")

d1=open.nc('ersst.185401.nc')
d2=open.nc('ersst.185402.nc')

test=create.nc('test.nc')

#define dimensions
dim.def.nc(test,"time",1)
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
var.def.nc(test,"err","NC_SHORT",c(3,2,1,0))

#add fillValues attribute to variables
att.put.nc(test,'sst','_FillValue',"NC_SHORT",-999)
att.put.nc(test,'anom','_FillValue',"NC_SHORT",-999)
att.put.nc(test,'err','_FillValue',"NC_SHORT",-999)

#add values to variable
var.put.nc(test,'time',var.get.nc(d1,'time'))
var.put.nc(test,'zlev',var.get.nc(d1,'zlev'))
var.put.nc(test,'lat',var.get.nc(d1,'lat'))
var.put.nc(test,'lon',var.get.nc(d1,'lon'))
var.put.nc(test,'sst',var.get.nc(d1,'sst',collapse=F))
var.put.nc(test,'anom',var.get.nc(d1,'anom',collapse=F))
var.put.nc(test,'err',var.get.nc(d1,'err',collapse=F))


close.nc(test)
close.nc(d1)
close.nc(d2)