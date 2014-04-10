## this code is for filtering the anomaly SST data using Butterworth
## bandpass filter with corner periodicity 1 - 10 years
load('SSTXiao.RData')
ndim=dim(sst)
nt=ndim[4]
nlat=ndim[2]
nlon=ndim[1]

## set up the reference location
inilon=320 
inilat=60 
ini_ind=c(which(lon==inilon),which(lat==inilat))

## take one time series for filtering
tt=seq(from=1854,to=2014,length=nt)
yy=anom[ini_ind[1],ini_ind[2],1,]
dev.new()
plot(tt,yy,type='l',axes=F)
axis(1,at=pretty(range(tt),n=10))
axis(2,at=pretty(range(yy)))

## filter the sample time series according to the range of periodicity interested in
source('~/Documents/TimeSeries/PaleoAnalyze/filterPaleo.R')
temp=as.data.frame(cbind(tt,yy))
names(temp)=c('t','y')
templ=list(temp)
fil=c(0.1,0.5)
tempfil=bwfilter(templ,cut=fil,type='pass',PLOT=T)

## prep for the following loop
anomfil2to10=sst
t=tt
NAs=rep(NA,length(t))

## loop thru all time series for filtering
# print(paste('loop starts at',Sys.time()))
# for (i in 1:nlon) {
#         for (j in 1:nlat) {
#                 if(any(is.na(anom[i,j,1,])) | any(sst[i,j,1,]==-180)) {
#                         anomfil2to10[i,j,1,]=NAs
#                         next
#                 }
#                 y=anom[i,j,1,]
#                 temp=as.data.frame(cbind(t,y))
#                 templ=list(temp)
# 
#                 tempfil=bwfilter(templ,cut=fil,type='pass',PLOT=F)
#                 anomfil2to10[i,j,1,]=tempfil[[1]][[1]]$y
#         }
# }
# print(paste('loop ends at',Sys.time())) # used ~6min
# save(anomfil2to10,file='anomfil2to10.RData')
load('anomfil2to10.RData')

## sample output
inilon=320
inilat=40
ini_ind=c(which(lon==inilon),which(lat==inilat))

dev.new()
plot(t,anom[ini_ind[1],ini_ind[2],1,],type='l',col='grey')
lines(t,anomfil2to10[ini_ind[1],ini_ind[2],1,],col='red')

dev.new()
image(lon,lat,anomfil2to10[,,1,1])
dev.new()
image(lon,lat,anom[,,1,1])
