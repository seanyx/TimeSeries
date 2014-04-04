bwfilter<-function(data,cut,type='low',PLOT=T) {
	### butterworth filter function
	library(signal)
	source('http://www.unc.edu/~yangxiao/PaleoR/PaleoSearch/xyrange.R')
	dt=mean(diff(data[[1]][,1]))
	fNy=.5/dt
	N=length(data)
	dataout=vector('list',N)
	res=vector('list',N)
	
	if (type=='low' | type=='high') {
		
		fc=cut
		fil.but=butter(4,fc/fNy,type=type,plane='z')
		datadm=data
		for (i in 1:N) {
			m=mean(data[[i]][,2])
			datadm[[i]][,2]=data[[i]][,2]-m
			dataout[[i]]=data.frame(t=datadm[[1]][,1],y=filtfilt(fil.but,datadm[[i]][,2]))
			dataout[[i]][,2]=dataout[[i]][,2]+m
			rest=data[[i]][,2]-dataout[[i]]$y
			res[[i]]=data.frame(t=data[[1]][,1],y=rest)
			}
		if (type=='low') filter.title=paste('lowpass filter:',fc)
		if (type=='high') filter.title=paste('highpass filter:',fc)
	}	
	
		
 	if (type=='pass') {
		 fp=cut
		 fil.but=butter(4,fp[2]/fNy,type='low',plane='z')
		 datadm=data
		 for (i in 1:N) {
			 m=mean(data[[i]][,2])
			 datadm[[i]][,2]=data[[i]][,2]-m
			 dataout[[i]]=filtfilt(fil.but,datadm[[i]][,2])
			 }
			 fil.but=butter(4,fp[1]/fNy,type='high',plane='z')
		 for (i in 1:N) {
			 dataout[[i]]=data.frame(t=data[[1]][,1],y=filtfilt(fil.but,dataout[[i]]))
			 dataout[[i]][,2]=dataout[[i]][,2]+m
			 rest=data[[i]][,2]-dataout[[i]]$y
			 res[[i]]=data.frame(t=data[[1]][,1],y=rest)
			 }
		 filter.title=paste('bandpass filter:',fp[1],'to',fp[2])
	 }
	if (PLOT) {
		dev.new()
		plot(y~t,data[[1]],type='l',col='grey',main=filter.title,xlab='Year BP')
		## lines(y~t,res[[1]],col='blue',lwd=1)
		lines(y~t,dataout[[1]],col='red',lwd=1.5)
	}
	invisible(list(dataout,res))
}
