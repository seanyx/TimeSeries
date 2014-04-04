#### search between grip and byrd

#load('/Volumes/YANGXIAO/alldata/gripbyrdsync')
setwd('/Users/yangxiao/Dropbox/Research/PaleoSearch backup/alldata')
load('gripbyrdsync')

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
		fpi2[i,(j+i-1)]=y[[2]]*(sin(y[[1]]*pi))^3
	}
	
}

f2=f1
quartz(width=12*1.2,height=12)
filled.contour(f1,f2,fpi2,color=rainbow,ann=T,axes=T)
title(xlab='lower frequency',ylab='higher frequency',main='pi/2 phase shift between GRIP and BYRD (methane sync) data\n Search range: 1/10ky to 1/800/n (50 intervals) with sin^3',cex.lab=1.5)

### output table
p=locator(type='p',pch=20,col='white',cex=1.3)
text(p$x,p$y,labels=as.character(1:length(p$x)),pos=3,col='white',cex=1.3)

n3=length(p$x)
corph=vector('list',length=n3)
for (i in 1:n3) {
	f=c(p$x[i],p$y[i])
	datafil=bwfilter(data,cut=c(f[1],f[2]),type='pass',PLOT=F)
	data2=datafil[[1]]
	corph[[i]]=fitPhase(data2[[1]]$y,data2[[2]]$y,N=200,PLOT=F)
}
k=unlist(corph)
k1=matrix(k,ncol=2,byrow=T)
fpi21=k1[,2]*(sin(k1[,1]*pi))^3
k2=cbind(p$x,p$y)
k=data.frame(k2,k1,fpi21)
names(k)=c('bp_low','bp_high','phase shift','maxcorr','fpi2')


#### plot examples
f=c(0.00038,0.00051)
datafil=bwfilter(data,cut=c(f[1],f[2]),type='pass',PLOT=F)
data2=datafil[[1]]
r=fitPhase(data2[[1]]$y,data2[[2]]$y,N=200,PLOT=F)
r1=phshift(data2[[2]]$y,r[[1]])
quartz(width=12,height=6)
plot(y~t,data2[[1]],type='l',col='black',xlab='Year BP',ylab='')
lines(y~t,data2[[2]],col='grey')
lines(data2[[1]]$t,r1,col='red')
title(main=paste('Red is', format(r[[1]]/.5,digits=3) ,'(*pi/2) phase shift of grey, correlation between red and black:',format(r[[2]],digits=3)))
legend('bottomright',bty='n',legend=c('GRIP','BYRD'),col=c('black','grey'),lty=1)