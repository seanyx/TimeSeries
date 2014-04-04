plotMatrix=function(m,cols){
	Ncol=length(cols)
	m.min=min(m)
	m.max=max(m)
	m.dim=dim(m)

	dx=1/(m.dim[2]-1)/2
	dy=1/(m.dim[1]-1)/2
	axis2=seq(0-dy,1+dy,length=m.dim[1]+1)
	axis3=seq(0-dx,1+dx,length=m.dim[2]+1)
	label.axis2=rev(seq(0,m.dim[1],by=1))
	label.axis3=seq(0,m.dim[2],by=1)
	D=range(m)

	d4=1/(Ncol-1)/2
	legend.abline=seq(0-d4,1+d4,length=Ncol+1)
	k=(1+d4+d4)/(D[2]-D[1])
	key=seq(D[1],D[2],length=5)
	axis4=(key-D[1])*k-d4



	layout(matrix(data=c(1,2),nrow=1,ncol=2),width=c(4,1),height=c(1,1))
	par(mar = c(3,2,3,0))
	# cols=rainbow(n=10,s=1,start=0,end=0.8)
	image(prepMatrix(m),col=cols,ann=F,axes=F)
	axis(2,at=axis2,label=label.axis2)
	axis(3,at=axis3,label=label.axis3)
	abline(v=axis3,lty=3)
	abline(h=axis2,lty=3)
	box()

	par(mar = c(3,1,3,5))
	image(matrix(1:Ncol,ncol=Ncol,nrow=1),col=cols,axes=F)
	axis(4,at=axis4,label=format(key,digits=3))
	abline(h=legend.abline,lty=3)
	box()
	layout(1)
	par(mar=c(5,4,4,2)+0.01)
}

