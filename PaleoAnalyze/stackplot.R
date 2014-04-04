stackplot<-function(data,xlim,env='laptop',point=F,mdata,main,color) {
	
	source('http://www.unc.edu/~yangxiao/PaleoR/PaleoSearch/xyrange.R')
	if (!missing(mdata)) data=data[intersect(mdata[,1],names(data))]
	
	if (missing(xlim)) {
		print('Please input the plot range:')
		m=scan(what='character',sep='&',nlines=1,strip.white=T)
		if (m=='all') {	
			t=xyrange(data)
			xlim=as.numeric(unlist(t['xo']))
			}
		if (m!='all') {
			xlim=as.numeric(unlist(strsplit(m,'to')))
			l=c(paste('I:',xlim[1],sep=''),paste('J:',xlim[2],sep=''))
			mdata=query(mdata,t=l,output=F)
		}
	}
	

	if (!missing(mdata)) data=data[intersect(mdata[,1],names(data))]
	
	dname=names(data)
	N=length(data)
	
	
	## set up working environment
	opar=par(no.readonly=T)
	if (env=='awesome12') { 
		color=hcl(h=seq(0,360,length=N),c=70,l=100)
 		par(mai=c(0,1,1,1),ann=F,mfrow=c(N+2,1),xaxs='r',yaxs='r',xaxt='s',bty='n',bg='black',fg='white',font=2,col.axis='white',col.lab='white',cex.axis=1.5)	
		} 
 	if (env=='laptop' | env=='imac') {	
		if (missing(color)) color=hcl(h=seq(0,360,length=N),c=90,l=70)
		par(mai=c(0,1,1,1),ann=F,mfrow=c(N+2,1),xaxs='r',yaxs='r',xaxt='s',bty='n')	
		}

### title line
plot(y~t,data[[1]],type='n',axes=F,col=color[1],xlim=xlim)
if (!missing(main)) title(main=main)
axis(side=3,at=format(seq(xlim[1],xlim[2],length=10),digits=5))
par(mai=c(0,1,0,1))

	for (i in 1:N) {	
	plot(y~t,data[[i]],type='l',axes=F,col=color[i],xlim=xlim)
	if(point) points(y~t,data[[i]],col=color[i])
	axis(side=-2*(i%%2)+4,at=pretty(data[[i]][,2]))
	mtext(side=-2*(i%%2)+4,dname[i],line=3)
	}
	
###xlab and xaxis
	par(mai=c(1,1,0,1))
	plot(y~t,data[[1]],type='n',axes=F,col=color[1],xlim=xlim)
	axis(side=1,at=format(seq(xlim[1],xlim[2],length=10),digits=5))
	mtext(side=1,"Year BP",line=3)
	
	par(opar)
}
