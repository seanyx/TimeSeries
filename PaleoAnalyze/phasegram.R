phasegram<-function(data,frange,nf,tw,toverlap=.5,M=50,taper=T) {
	# phasegram of between two time series
	source('http://www.unc.edu/~yangxiao/PaleoR/PaleoAnalyze/matrix2array.R')
	library(ggplot2)
	
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

	fp=ef
	for (i in 1:length(tind)) {
		data1=data
		data1[[1]]=data1[[1]][(tind[i]:(tind[i]+twd-1)),]
		data1[[2]]=data1[[2]][(tind[i]:(tind[i]+twd-1)),]
	for (j in 1:(length(f)-1)) {
		f2=c(f[j],f[j+1])
		datafil=bwfilter(data1,cut=c(f2[1],f2[2]),type='pass',PLOT=F)
		data2=datafil[[1]]
		
		ef[j,i]=sum((data2[[1]]$y-mean(data2[[1]]$y))^2+(data2[[2]]$y-mean(data2[[2]]$y))^2)/etotal

		if(sd(data2[[1]]$y)>1e-9 & sd(data2[[2]]$y)>1e-9) {
			y=fitPhase(data2[[1]]$y,data2[[2]]$y,N=M,PLOT=F)
			fp[j,i]=y[[2]]*(sin(y[[1]][1]*pi))^4
			}
	}
}

ef1=apply(ef,MARGIN=1,sum)

#test weighed by power
# for (i in 1:(length(f)-1)) {
#	fp[i,]=fp[i,]*sqrt(ef1[i])
#	}

dev.new()
filled.contour(data[[1]]$t[tind+floor(twd/2)],f[-nf],t(fp),color=rainbow,ann=T,axes=T)
title(main=paste('Phase shift between',na[1],'and',na[2],'\n with time window',twd,'pt and',nf,'frequency intervals and time window overlap',toverlap,'via M1'))

##
ls=LSspecgram(data[[1]],freq=seq(frange[1],frange[2],length=nf-1),wd=twd,overlap=toverlap)
dev.new()
image(data[[1]]$t[tind+floor(twd/2)],f[-nf],t(fp))
contour(ls[[1]],ls[[2]],ls[[3]],add=T)
print(ls[[1]])

# test ggplot
dev.new()
x=data[[1]]$t[tind+floor(twd/2)]
y=f[-nf]
m=fp
re=matrix2array(m,y,x)
re=as.data.frame(re)
names(re)=c('x','y','z')
p=ggplot(re,aes(x,y,z=z))
p+geom_tile(aes(fill=z))+scale_fill_gradient2(low='darkblue',high='darkred',mid='white',midpoint=.3)

#dev.new()
#plot(f[-nf],ef1,type='h',col='black',lwd=2,ann=F,axes=F)

}
