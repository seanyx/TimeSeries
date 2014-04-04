plotlog <- function(x,y,base=10,xlim,ylim,grid=T,lcol,add=F,set.frame=F,...) {
	if(missing(lcol)) lcol='black'
	if (!add) {
		xrange=range(x)
		yrange=range(y)
		
		if(!missing(xlim)) xrange=xlim
		if(!missing(ylim)) yrange=ylim

		if(xrange[1]<0 | yrange[1]<0) stop('range < 0')

		xlog=log(xrange,base)
		ylog=log(yrange,base)

		xtickr=floor(xlog[1]):ceiling(xlog[2])
		ytickr=floor(ylog[1]):ceiling(ylog[2])

		nx=length(xtickr)-1
		ny=length(ytickr)-1

		xtick=vector(mode='numeric')
		ytick=xtick
		for ( i in 1:nx) {
			temp=log(seq(base^xtickr[i],base^xtickr[i+1]-base^xtickr[i],by=base^xtickr[i]),base)
			xtick=c(xtick,temp)
		}

		for (i in 1:ny) {
			temp=log(seq(base^ytickr[i],base^ytickr[i+1]-base^ytickr[i],by=base^ytickr[i]),base)
			ytick=c(ytick,temp)
		}
		Xtick=xtickr
		Ytick=ytickr
		labX=as.character(base^Xtick)
		labY=as.character(base^Ytick)
		
		xlim=xlog
		ylim=ylog

		plot(xlim,ylim,type='n',xlim=xlim,ylim=ylim,axes=F,...)
		axis(1,at=xtick,label=F,lwd.ticks=.5)
		axis(1,at=Xtick,label=labX,lwd.ticks=1.5)
		axis(2,at=ytick,label=F,lwd.ticks=.5)
		axis(2,at=Ytick,label=labY,lwd.ticks=1.5)
		if (grid) {
			abline(v=xtick,col='grey',lwd=.5)
			abline(h=ytick,col='grey',lwd=.5)
		}
	}

	if(!set.frame) {
		X=log(x,base)
		Y=log(y,base)
		lines(X,Y,col=lcol,lwd=1.5)
	}

}
