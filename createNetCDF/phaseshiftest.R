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


fitPhase<-function(y1,y2,N,ReAll=F) {
	ycor=vector(length=N)
	phrange=seq(-1,1,length=N)
	for (i in 1:N) {
		y2m=phshift(y2,pd=phrange[i])
		ycor[i]=cor(y1,y2m)
	}
	ycormax=max(ycor)
	ind=which(ycor==ycormax)
	phmax=phrange[ind]
	if(any(phmax==1)) phmax=1
	
	if (!ReAll) return(list(phmax,ycormax))
	
	if (ReAll) return(list(phmax,ycormax,phrange,ycor))	
	
}