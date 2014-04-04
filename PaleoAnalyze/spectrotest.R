## check spectrogram from Dr. Rial
library(seewave)
library(signal)
source('http://www.unc.edu/~yangxiao/PaleoR/PaleoAnalyze/filterPaleo.R')
load(url('http://www.unc.edu/~yangxiao/PaleoR/yxdata/ngripdomecsync'))
data=data4

r1=bwfilter(data,cut=c(1/6000,1/1000),type='pass',PLOT=T)
datafil=r1[[1]]
dt=mean(diff(datafil[[1]]$t))

y=spectro(datafil[[2]]$y,1/dt,wl=100,overlap=.3,zp=16)

indf=which(y$freq<=1e-6)
filled.contour(y$time+data[[1]]$t[1],y$freq[indf]*1000,t(y$amp[indf,]),color=rainbow,ylim=c(0,1e-3),nlevels=20)
