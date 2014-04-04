createNetCDF2d<-function(path,filename,time,height,lat,lon,var.names,variables,null_value=NULL,refs=NULL) {
	## variables should be a list with each element stores value for one variable
	## in this function, variables depend on [longitude, latitude, elevation, time], so each element of variables should be a high-dimensional array with dimension [1:Nlon, 1:Nlat, 1:Nheight, 1:Ntime], so that variables[[i]][,,1,1] is a full grid of the global.
	library(RNetCDF)
	setwd(path)
	
	file=create.nc(filename)
	Ntime=length(time)
	Nlat=length(lat)
	Nlon=length(lon)
	Nheight=length(height)
	Nvar=length(var.names)
	
	
	##define dimensions
	dim.def.nc(file,'time',Ntime)
	dim.def.nc(file,'zlev',Nheight)
	dim.def.nc(file,'lat',Nlat)
	dim.def.nc(file,'lon',Nlon)
	
	#define variables
	var.def.nc(file,"time","NC_FLOAT",0)
	var.def.nc(file,"zlev","NC_FLOAT",1)
	var.def.nc(file,"lat","NC_FLOAT",2)
	var.def.nc(file,"lon","NC_FLOAT",3)
	for (i in 1:Nvar) {
		var.def.nc(file,var.names[i],"NC_DOUBLE",c(3,2,1,0))
	}
	
	#add FillValue if applicable
	if (!is.null(null_value)) {
		if (Nvar!=length(null_value)) stop('Null value vector and Variable vector are of different length.')
		for (i in 1:length(null_value)) {
			## add fillValues attribute to variables
			att.put.nc(file,var.names[i],'_FillValue',"NC_DOUBLE",null_value[i])
		}
	}
	
	#add values to variable
	var.put.nc(file,'time',time)
	var.put.nc(file,'zlev',height)
	var.put.nc(file,'lat',lat)
	var.put.nc(file,'lon',lon)
	for (i in 1:Nvar) {
		var.put.nc(file,var.names[i],variables[[i]])
	}
	
	if(!is.null(refs)) att.put.nc(file,"NC_GLOBAL","Data source","NC_CHAR",refs)
		
	close.nc(file)
}