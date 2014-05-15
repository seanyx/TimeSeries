## RK4 solution for Van der Pol oscillator (Rial & Saha, 2011), no noise situation
## 
f1=function(a,b,c,x,y,t) -b*x^3+a*x-y+c
f2=function(omega,A,x,t,fm=10000) omega^2*x-A*sin(2*pi*fm*t)

## input parameters
fm=1/10000
h=10
fc=1/1500
omega=2*pi*fc
A=1e-6
x0=0.001
y0=0.001
b=1
c=0

N=4000
xn=vector(length=N+1)
yn=vector(length=N+1)







a=seq(0.008,0.12,by=0.01)
Na=length(a)
t=seq(0,h*N,by=h)
dev.new()
plot(t,xn,type='l',ylim=c(-0.20,(Na+1)*0.5),ann=F,axes=F)
title(main=para)
axis(2,at=seq(from=0,by=0.3,length=Na),labels=as.character(a))
axis(1,pretty(t))

for (d in 1:Na) {
	xn[1]=x0
	yn[1]=y0
	for (i in 1:N) {
	
		m=h*c(f1(a[d],b,c,xn[i],yn[i],i*h), f2(omega,A,xn[i],i*h))
		n=h*c(f1(a[d],b,c,xn[i]+1/2*m[1],yn[i]+1/2*m[2],i*h+1/2*h), f2(omega,A,xn[i]+1/2*m[1],i*h+1/2*h))
		j=h*c(f1(a[d],b,c,xn[i]+1/2*n[1],yn[i]+1/2*n[2],i*h+1/2*h), f2(omega,A,xn[i]+1/2*n[1],i*h+1/2*h))
		k=h*c(f1(a[d],b,c,xn[i]+j[1],yn[i]+j[2],i*h+h), f2(omega,A,xn[i]+j[1],i*h+h))
	
		xn[i+1]=xn[i]+1/6*(m[1]+2*n[1]+2*j[1]+k[1])
		yn[i+1]=yn[i]+1/6*(m[2]+2*n[2]+2*j[2]+k[2])
	
	}
	lines(t,xn+rep(0.5*d,N+1))
	lines(t,0.1*sin(2*pi*fm*t),col='red')
}
para=paste("A=",A,"; Tm=",1/fm,"; Tc=",1/fc,"; b=",b,"; c=",c,"; h=",h,"; N=",N,sep='')


# dev.new()
# plot(t,yn,type='l')
# dev.new()
# plot(xn,yn,type='l')
