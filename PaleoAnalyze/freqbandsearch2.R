#### search between gisp2 and byrd

#load('/Volumes/YANGXIAO/alldata/gisp2byrdsync')
setwd('/Users/yangxiao/Dropbox/Research/PaleoSearch backup/alldata')
load('gisp2byrdsync')
data=data2

source('phaseshiftest.R')
source('filterPaleo.R')
source('xyrange.R')
source('interpPaleo.R')

fp=c(1/10000,1/800)
N=50
df=(fp[2]-fp[1])/(N-1)
f1=seq(fp[1],fp[2],length=N)
M=50
fpi2=matrix(nrow=N,ncol=N)


for (i in 1:N) {
	f2=seq(f1[i],fp[2],by=df)
	n2=length(f2)
	for (j in 1:(N+1-i)) {
		datafil=bwfilter(data,cut=c(f1[i],f2[j]),type='pass',PLOT=F)
		data2=datafil[[1]]
		y=fitPhase(data2[[1]]$y,data2[[2]]$y,N=M,PLOT=F)
		fpi2[i,(j+i-1)]=y[[2]]*(sin(y[[1]]*pi))^4
	}
	
}

f2=f1
quartz(width=12*1.2,height=12)
filled.contour(f1,f2,fpi2,color=rainbow,ann=T,axes=T)
title(xlab='lower frequency',ylab='higher frequency',main='pi/2 phase shift between GISP2 and BYRD (methane sync) data\n Search range: 1/10ky to 1/800/n (50 intervals) with sin^4',cex.lab=1.5)
