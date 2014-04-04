### spectrogram of phase shift

load('MDgripage')

data=data4
na=names(data)
N=length(data[[1]]$t)
frange=c(1/10000,1/800)
trange=range(data[[1]]$t)

tw1=20000
twd=length(which(data[[1]]$t-data[[1]]$t[1]<=tw1))
toverlap=.5
t1=floor(twd*toverlap)
tind=seq(1,N-twd,by=(twd-t1)) # time index
fwd=1/20000
foverlap=.6
f1=floor(fwd*foverlap)
nf=100
f=seq(frange[1],frange[2],length=nf) # exact frequency
M=50

ef=matrix(nrow=(nf-1),ncol=length(tind))
etotal=sum((data[[1]]$y-mean(data[[1]]$y))^2+(data[[2]]$y-mean(data[[2]]$y))^2)

fp=matrix(nrow=(length(f)-1),ncol=length(tind))
for (i in 1:length(tind)) {
	data1=data
	data1[[1]]=data1[[1]][(tind[i]:(tind[i]+twd-1)),]
	data1[[2]]=data1[[2]][(tind[i]:(tind[i]+twd-1)),]
	for (j in 1:(length(f)-1)) {
		f2=c(f[j],f[j+1])
		datafil=bwfilter(data1,cut=c(f2[1],f2[2]),type='pass',PLOT=F)
		data2=datafil[[1]]
		
		ef[j,i]=sum((data2[[1]]$y-mean(data2[[1]]$y))^2+(data2[[2]]$y-mean(data2[[2]]$y))^2)/etotal
		
		if (sd(data2[[1]]$y)>1e-9 & sd(data2[[2]]$y)>1e-9) {
		y=fitPhase(data2[[1]]$y,data2[[2]]$y,N=M,PLOT=F)
		fp[j,i]=y[[2]]*(sin(y[[1]][1]*pi))^4
		}
	}
}

ef1=apply(ef,MARGIN=1,sum)

dev.new()
#pdf(file=paste('/plots/',na[1],'_phase_spectrogram.pdf',sep=''),width=15,height=15)
filled.contour(data[[1]]$t[tind+floor(twd/2)],f[-nf],t(fp),color=heat.colors,ann=T,axes=T)
title(main=paste('Phase shift between',na[1],'and',na[2],'\n with time window',twd,'pt and',nf,'frequency intervals and time window overlap',toverlap))
#dev.off()
