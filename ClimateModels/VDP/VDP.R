## RK4 solution for Van der Pol oscillator (Rial & Saha, 2011), no noise situation
## 

## import insolation
## link: ftp://ftp.ncdc.noaa.gov/pub/data/paleo/insolation/
inso=read.table("~/Documents/TimeSeries/ClimateModels/VDP/orbit91",skip=3)[1:101,6]
inso=rev(inso)
inso=(inso-mean(inso))/sd(inso)

f1=function(a,b,c,x,y,t,D,h) -b*x^3+a*x-y+c+D*sqrt(h)*rnorm(1)
f2=function(omega,A,x,t,inso,DC) omega^2*x-A*(approx(1:101,inso,t/1000,rule=2)[[2]]+DC)

## input parameters
## A and a together control the amount of frequency modulation
## b mostly controls the amplitude
## increase a will also increase the natural periodicity
h=1
fc=1/1500
omega=2*pi*fc
A=2.1e-6
x0=0.001
y0=0.001
a=0.1
b=1
c=0
D=1e-3
DC=0.6

N=100000
xn=vector(length=N+1)
yn=vector(length=N+1)

xn[1]=x0
yn[1]=y0

for (i in 1:N) {
	
	m=h*c(f1(a,b,c,xn[i],yn[i],i*h,D,h), f2(omega,A,xn[i],i*h,inso,DC))
	n=h*c(f1(a,b,c,xn[i]+1/2*m[1],yn[i]+1/2*m[2],i*h+1/2*h,D,h), f2(omega,A,xn[i]+1/2*m[1],i*h+1/2*h,inso,DC))
	j=h*c(f1(a,b,c,xn[i]+1/2*n[1],yn[i]+1/2*n[2],i*h+1/2*h,D,h), f2(omega,A,xn[i]+1/2*n[1],i*h+1/2*h,inso,DC))
	k=h*c(f1(a,b,c,xn[i]+j[1],yn[i]+j[2],i*h+h,D,h), f2(omega,A,xn[i]+j[1],i*h+h,inso,DC))
	
	xn[i+1]=xn[i]+1/6*(m[1]+2*n[1]+2*j[1]+k[1])
	yn[i+1]=yn[i]+1/6*(m[2]+2*n[2]+2*j[2]+k[2])
	
}

para=paste("A=",A,"; a=",a,"; Tc=",1/fc,"; b=",b,"; c=",c,"; h=",h,"; N=",N,"; D=",D,'; DC=',DC,sep='')
t=seq(0,h*N,by=h)

dev.new(width=12,height=6)
par(mfrow=c(2,1))
plot(t,xn,type='l')
lines(t,approx(1:101,inso+DC,t/1000,rule=2)[[2]]/20,col='red')
abline(h=0,col='green')
title(main=para)

ngrip=read.table("~/Documents/TimeSeries/ClimateModels/VDP/ngrip-d18o-50yr.txt",header=F,skip=80,col.names=c('t','y'),sep='')
ngrip=ngrip[,2]
tt=seq(0,100000,by=50)
ngrip=ngrip[(1:length(tt))*2]
plot(rev(tt),ngrip,type='l',main='NGRIP')


par(mfrow=c(2,1))
# dev.new()
# plot(t,yn,type='l')
# dev.new()
# plot(xn,yn,type='l')
