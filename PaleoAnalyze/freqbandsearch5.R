setwd('/Users/yangxiao/Dropbox/Research/PaleoSearch backup')
ngrip=read.table('alldata/NGRIP_sync.txt',skip=0)
ngrip =data.frame(t= ngrip[,1],y= ngrip[,2])
domec=read.table('alldata/DomeC_sync.txt',skip=0)
domec=data.frame(t= domec[,1],y= domec[,2])
data4=list(ngrip=ngrip,domec=domec)
data4[[2]]=data4[[2]][1:length(data4[[1]]$t),]
save(data4,file='alldata/ngripdomecsync')


#### search between benthic and planktonic in synchronized NGRIP and DOME C

#load('/Volumes/YANGXIAO/alldata/gisp2byrdsync')
setwd('/Users/yangxiao/Dropbox/Research/PaleoSearch backup/alldata')
load('ngripdomecsync')
data=data4

stackplot(data=data4)

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
		if (sd(data2[[1]]$y)>1e-9 & sd(data2[[2]]$y)>1e-9) {
		y=fitPhase(data2[[1]]$y,data2[[2]]$y,N=M,PLOT=F)
		fpi2[i,(j+i-1)]=y[[2]]*(sin(y[[1]]*pi))^4
		}
	}
	
}

f2=f1
quartz(width=12*1.2,height=12)
filled.contour(f1,f2,fpi2,color=rainbow,ann=T,axes=T)
title(xlab='lower frequency',ylab='higher frequency',main='pi/2 phase shift between NGRIP and DOME C (methane sync) data\n Search range: 1/10ky to 1/800 (50 intervals) with sin^4',cex.lab=1.5)




f=c(0.00036,0.0006)
datafil=bwfilter(data,cut=c(f[1],f[2]),type='pass',PLOT=F)
data2=datafil[[1]]
r=fitPhase(data2[[1]]$y,data2[[2]]$y,N=200,PLOT=F)
r1=phshift(data2[[2]]$y,r[[1]])
quartz(width=12,height=6)
plot(y~t,data2[[1]],type='l',col='black',xlab='Year BP',ylab='')
lines(y~t,data2[[2]],col='grey')
lines(data2[[1]]$t,r1,col='red')
title(main=paste('Red is', format(r[[1]]/.5,digits=3) ,'(*pi/2) phase shift of grey, correlation between red and black:',format(r[[2]],digits=3)))
legend('bottomright',bty='n',legend=c('MDbenthic','MDplanktonic'),col=c('black','grey'),lty=1)
