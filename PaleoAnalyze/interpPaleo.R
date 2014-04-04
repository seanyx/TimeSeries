interpPaleo<-function(sig,xout='min',method='linear') {
	xdmin=min(diff(sig[,1]))
	xdmean=mean(diff(sig[,1]))
	xmin=min(sig[,1])
	xmax=max(sig[,1])
	if (xout=='min') {
		xomin=seq(xmin,xmax,by=xdmin)
		y=approx(sig[,1],sig[,2],xomin,method=method)
	}
	if (xout=='mean') {
		xomean=seq(xmin,xmax,by=xdmean)
		y=approx(sig[,1],sig[,2],xomean,method=method)
	}
	if (xout!='mean' & xout!='min') y=approx(sig[,1],sig[2,],xout,method=method)
	#names(y)=c('t','y')
	return(data.frame(t=y[[1]],y=y[[2]]))
}

interPaleoAll<-function(data,xout='mean',method='linear') {
	source('http://www.unc.edu/~yangxiao/PaleoR/PaleoSearch/xyrange.R')
	N=length(data)
	ran=xyrange(data)
	xlim=ran['xi']$xi
	if(is.character(xout) & xout=='min') xinter=ran['xInterval']$xInterval[1]
	if(is.character(xout) & xout=='mean') xinter=ran['xInterval']$xInterval[2]
	if(xout!='mean' & xout!='min') xinter=mean(diff(xout))
	xouti=seq(xlim[1],xlim[2],by=xinter)
	dataout=vector(mode='list',length=N)
	
	for (i in 1:N) {
		y=approx(data[[i]][,1],data[[i]][,2],xout=xouti,method=method)	
		dataout[[i]]=data.frame(t=y[[1]],y=y[[2]])
	}
	return(dataout)
}
