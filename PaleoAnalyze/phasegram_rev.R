phasegram_rev<-function(data,frange,nf,tw,toverlap=.5,M=50,taper=T) {
	# revised phasegram of between two time series

	if (taper) {
		data[[1]]$y=spec.taper(data[[1]]$y,p=.1)
		data[[2]]$y=spec.taper(data[[2]]$y,p=.1)
	}	
	na=names(data)
	trange=range(data[[1]]$t)
	N=length(data[[1]]$t)
	twd=length(which(data[[1]]$t-data[[1]]$t[1]<=tw))
	t1=floor(twd*toverlap)
	tind=seq(1,N-twd,by=(twd-t1))
	f=seq(frange[1],frange[2],length=nf)

	ef=matrix(nrow=(nf-1),ncol=length(tind))
	mean1=mean(data[[1]]$y)
	mean2=mean(data[[2]]$y)
	etotal=sum((data[[1]]$y-mean1)^2+(data[[2]]$y-mean2)^2)

	ef2=vector(length=nf-1)
	fp=ef

	for (i in 1:(nf-1)) {
		f2=c(f[i],f[i+1])
		datafil=bwfilter(data,cut=f2,type='pass',PLOT=F)
		data1=datafil[[1]]
		ef2[i]=sum((data1[[1]]$y-mean1)^2+(data1[[2]]$y-mean2)^2)


		for (j in 1:length(tind)) {
			data2=data1
			data2[[1]]=data2[[1]][(tind[j]:(tind[j]+twd-1)),]
			data2[[2]]=data2[[2]][(tind[j]:(tind[j]+twd-1)),]
			
			ef[i,j]=sum((data2[[1]]$y-mean(data2[[1]]$y))^2+(data2[[2]]$y-mean(data2[[2]]$y))^2)

			if(sd(data2[[1]]$y)>1e-9 & sd(data2[[2]]$y)>1e-9) {
				y=fitPhase(data2[[1]]$y,data2[[2]]$y,N=M,PLOT=F)
				fp[i,j]=y[[2]]*(sin(y[[1]][1]*pi))^4
			}

		}
	}
	
ef1=apply(ef,MARGIN=1,sum)
u=ef*1/ef1*(ef2*etotal/sum(ef2))
#test weighed by power
# for (i in 1:(length(f)-1)) {
#	fp[i,]=fp[i,]*sqrt(ef1[i])
#	}

dev.new()
filled.contour(data[[1]]$t[tind+floor(twd/2)],f[-nf],t(fp),color=rainbow,ann=T,axes=T)
title(main=paste('Phase shift between',na[1],'and',na[2],'\n with time window',twd,'pt and',nf,'frequency intervals and time window overlap',toverlap,'via M2'))

#dev.new()
#filled.contour(data[[1]]$t[tind+floor(twd/2)],f[-nf],t(u),color=heat.colors,ann=T,axes=T)



#dev.new()
#plot(f[-nf],ef1,type='h',col='black',ann=F,axes=F)

}
