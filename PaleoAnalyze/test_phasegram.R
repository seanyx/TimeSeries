source('http://www.unc.edu/~yangxiao/PaleoR/PaleoAnalyze/phasegram_rev.R')
source('http://www.unc.edu/~yangxiao/PaleoR/PaleoAnalyze/phasegram.R')
source('http://www.unc.edu/~yangxiao/PaleoR/PaleoAnalyze/filterPaleo.R')
source('http://www.unc.edu/~yangxiao/PaleoR/PaleoAnalyze/phaseshiftest.R')

load(url('http://www.unc.edu/~yangxiao/PaleoR/yxdata/gripbyrdsync'))
load(url('http://www.unc.edu/~yangxiao/PaleoR/yxdata/gisp2byrdsync'))
load(url('http://www.unc.edu/~yangxiao/PaleoR/yxdata/ngripdomecsync'))
load(url('http://www.unc.edu/~yangxiao/PaleoR/yxdata/MDgripage'))
load(url('http://www.unc.edu/~yangxiao/PaleoR/yxdata/TR163_22')) #TR163_22


t=seq(10000,80000,by=50)
f=c(1/1500,1/6000)
sg1=sin(2*pi*f[1]*t)+sin(2*pi*f[2]*t)
sg2=cos(2*pi*f[1]*t)+cos(2*pi*f[2]*t)


#
sg1=sin(2*pi*f[1]*t)
sg2=cos(2*pi*f[1]*t)

# sharp change
sg1a=sin(2*pi*f[2]*t[1:500])
sg2a=cos(2*pi*f[2]*t[1:500])
sg1b=sin(2*pi*f[1]*t[501:1401])
sg2b=cos(2*pi*f[1]*t[501:1401])

sg1=c(sg1a,sg1b)
sg2=c(sg2a,sg2b)

# linear change in frequency
sg1=sin(2*pi*((1e-7/14)*(t-t[1])+f[2])*t)
sg2=cos(2*pi*((1e-7/14)*(t-t[1])+f[2])*t)




data=list(data.frame(t=t,y=sg1),data.frame(t=t,y=sg2))

phasegram(data=data,frange=c(1/10000,1/800),nf=100,tw=10000,toverlap=.3)
phasegram_rev(data=data,frange=c(1/10000,1/800),nf=100,tw=10000)



source('/Users/yangxiao/Dropbox/Research/UltimateR/LSFTnew.R')
datafil=bwfilter(data4,cut=c(1/10000,.0009),type='pass',PLOT=T)
data=datafil[[1]]
freq=seq(0,1e-3,length=100)
ls=LSspecgram(data[[1]],freq,wd=150)
dev.new()
filled.contour(ls[[1]],ls[[2]],ls[[3]],ylim=c(0,8e-4),plot.axes={axis(1);axis(2);contour(ls[[1]],ls[[2]],ls[[3]],nlevels=5,col='grey',add=T);contour(ls[[1]],ls[[2]],ls[[3]],nlevel=1,level=20,lwd=1.5,add=T)},nlevels=15,col=rainbow(15))






## compare two spectrograms

load(url('http://www.unc.edu/~yangxiao/PaleoR/yxdata/ngripdomecsync'))
datafil=bwfilter(data4,cut=.0009,type='low',PLOT=T)
data=datafil[[1]]
freq=seq(0,1e-3,length=100)
ls1=LSspecgram(data[[1]],freq,wd=150)
ls2=LSspecgram(data[[2]],freq,wd=150)
dev.new()
indt=which(ls1[[1]]<=75000)
filled.contour(ls1[[1]][indt],ls1[[2]],(ls2[[3]]-ls1[[3]])[indt,],ylim=c(0,8e-4),plot.axes={axis(1);axis(2);contour(ls1[[1]][indt],ls1[[2]],(ls2[[3]]-ls1[[3]])[indt,],nlevels=5,col='grey',add=T);contour(ls1[[1]][indt],ls1[[2]],(ls2[[3]]-ls1[[3]])[indt,],nlevel=1,level=20,lwd=1.5,add=T)},levels=c(-50,seq(-20,40,length=14),65),col=c('grey',rainbow(13),'white'))

dev.new()
indt=which(ls1[[1]]<=75000)
n=10
filled.contour(ls1[[1]][indt],ls1[[2]],(ls2[[3]]/(ls1[[3]]+n))[indt,],ylim=c(0,8e-4),plot.axes={axis(1);axis(2);contour(ls1[[1]][indt],ls1[[2]],(ls2[[3]]/(ls1[[3]]+n))[indt,],nlevels=5,col='grey',add=T);contour(ls1[[1]][indt],ls1[[2]],(ls2[[3]]/(ls1[[3]]+n))[indt,],nlevel=1,level=20,lwd=1.5,add=T)},levels=c(-50,seq(-20,40,length=14),65),col=c('black',rainbow(13),'white'))


dev.new()
indt=which(ls1[[1]]<=75000)
n=10
filled.contour(ls1[[1]][indt],ls1[[2]],(ls2[[3]]*(ls1[[3]]))[indt,],ylim=c(0,8e-4),plot.axes={axis(1);axis(2);contour(ls1[[1]][indt],ls1[[2]],(ls2[[3]]*(ls1[[3]]))[indt,],nlevels=5,col='grey',add=T);contour(ls1[[1]][indt],ls1[[2]],(ls2[[3]]*(ls1[[3]]))[indt,],nlevel=1,level=20,lwd=1.5,add=T)},levels=c(-50,seq(-20,40,length=14),65),col=c('black',rainbow(13),'white'))

