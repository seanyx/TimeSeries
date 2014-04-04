phshift<-function(y,pd,f=1){
	### return signal of pi*pd shift; pd=-1/2 is the Im part of hilbert transform
	library(seewave)
	
	ymean=mean(y)
	y=y-ymean
	hy=hilbert(y,f)
	my=abs(hy)*cos( (pd*pi) + (Arg(hy)) )
	my=my+ymean
	return(my)
}


fitPhase<-function(y1,y2,N,ReAll=F,PLOT=F) {
	ycor=vector(length=N)
	phrange=seq(-1,1,length=N)
	for (i in 1:N) {
		y2m=phshift(y2,pd=phrange[i])
		ycor[i]=cor(y1,y2m)
	}
	ycormax=max(ycor)
	ind=which(ycor==ycormax)
	phmax=phrange[ind]
	
	if (PLOT) {
		ylim1=c(min(min(y1),min(y2)),max(max(y1),max(y2)))
		ylim2=c(min(min(y1),min(phshift(y2,pd=phmax))),max(max(y1),max(phshift(y2,pd=phmax))))
		dev.new()
		par(mfrow=c(2,1))
		plot(y1,type='l',col='black',ann=F,axes=T,ylim=ylim1)
		lines(y2,col='blue')
		plot(y1,type='l',col='black',ann=F,axes=T,ylim=ylim2)
		lines(phshift(y2,pd=phmax),col='red')
		par(mfrow=c(1,1))
	}
	
	if (!ReAll) return(list(phmax,ycormax))
	
	if (ReAll) return(list(phmax,ycormax,phrange,ycor))	
	
}
