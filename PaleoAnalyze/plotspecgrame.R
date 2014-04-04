plot.specgram<-function(data,freq,wd,title,contour=F,n=15,ts=F,taper=F) {
## function to plot spectrograme using LS method
	source('http://www.unc.edu/~yangxiao/PaleoR/PaleoAnalyze/LSFTnew.R')
	N=length(data)
	dname=names(data)
	
	if (missing(title)) title=dname
	colrgb=as.vector(col2rgb('orange'))
	tcol=rgb(colrgb[1],colrgb[2],colrgb[3],alpha=100,max=255)
	for (i in 1:N) {
		if(taper) data[[i]]$y=spec.taper(data[[i]]$y-mean(data[[i]]$y))
		if(!taper) {}
		ls=LSspecgram(data[[i]],freq,wd)
		trange=range(ls[[1]])
		yrange=range(data[[i]]$y)
		if (contour) filled.contour(ls[[1]],ls[[2]],ls[[3]],plot.axes={axis(1);axis(2);contour(ls[[1]],ls[[2]],ls[[3]],nlevels=n[2],col='grey',add=T);contour(ls[[1]],ls[[2]],ls[[3]],nlevels=n[1],lwd=1.5,add=T)},nlevels=n[3],col=rainbow(n[3]),main=title[i],xlab='Years BP',ylab='Frequency')

		if (!contour) filled.contour(ls[[1]],ls[[2]],ls[[3]],nlevels=n,col=rainbow(n),main=title[i],xlab='Years BP',ylab='Frequency')

		if (ts) {
		quartz(width=12,height=6)
		plot(y~t,data[[i]],type='l',main=title[i],axes=F,ann=F)
		title(main=title[i],xlab='Years BP',ylab='')
		axis(1,pretty(trange))
		axis(2,pretty(yrange))
		rect(trange[1],yrange[1],trange[2],yrange[2],col=tcol,border=NA)
		}
	}
}
