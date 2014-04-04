data=load('~/Dropbox/class2014spring/Inverse Theory/hw2/HW2/profile.mat.RDATA')
x=G$x
t=G$t
N=length(t)
## model t = t0 + s2 * x
plot(x,t,type='p',pch=20)
sigma.t=0.1



## model fitting weighted by W
## t=Gb
## b=(t0,s2)'
W=diag(rep(1/sigma.t,N))
G=cbind(rep(1,N),x)
Gw=W %*% G
tw=W %*% t
b= solve( t(Gw) %*% Gw ) %*% t(Gw) %*% tw

## plot the data, the fitted model and the residuals
t.hat=G %*% b
res2=t-t.hat
library(pracma)
par(mfrow=c(2,1))
errorbar(x[,1],t[,1],xerr=NULL,yerr=rep(sigma.t,N),ann=F)
title(xlab='Distance x/km',ylab='Arrival time t/s',main='Seismic data and fitted linear model')
points(x,t.hat, pch=20, col='blue')
lines(x,t.hat, col='blue',lty=2)
legend('topleft',legend=c('Observed data','Model predicted data'),pch=c(1,20),col=c('black','blue'))

plot(x,t-t.hat,type='p',pch=15,xlab='Distance x/km',ylab='Residual arrival time t/s')
title(main='Residuals from fitted linear model')
par(mfrow=c(1,1))


## b. covariance matrix
cov.b=sigma.t^2 * solve( t(G) %*% G )
cor_s2t0=cov.b[1,2]/sqrt(cov.b[1,1]*cov.b[2,2]) # 2.43
cor.b=cov.b
diag(cor.b)=c(1,1)
cor.b[2,1]=cor_s2t0
cor.b[1,2]=cor_s2t0

print('The model parameter correlation matrix:')
cor.b

## Since the correlation between the two model parameters is close to -1, the projection of the ellipsoid into (t0, s2) plane shoule be needle-like with its long principal axis having a negative slope

## c. error ellipsoid

## to diagonalize covariance matrix C-1 gives the directions of the semiaxes for the error ellipsoid are
eigen.b=eigen(qr.solve(cov.b))
u=eigen.b$vectors
print('The directions of the semiaxes for the error ellipsoid are:')
u
print('with corresponding eigenvalues:')
lamda=eigen.b$values
lamda


delta=sqrt(qchisq(0.95,2))
semiaxes=delta/sqrt( lamda )
theta=seq(0,2*pi,by=0.01)
alpha=pi-atan(-u[2,2]/u[1,2])
Rotation=matrix(c(cos(alpha),-sin(alpha),sin(alpha),cos(alpha)),byrow=T,ncol=2)

r.temp=matrix(rep(0,times=2*length(theta)),ncol=2)

r.temp[,1]=semiaxes[2]*cos(theta)
r.temp[,2]=semiaxes[1]*sin(theta)

r= t(Rotation %*% t(r.temp))

plot(b[1]+r[,1],b[2]+r[,2],type='l' , ann=FALSE,col='red',lwd=2, asp = 1)
# alternative: library(ellipse); dev.new(); plot(ellipse(cor.b,scale=sqrt(diag(cov.b)),centre=b,level=.95),type='l')
points(b[1],b[2],pch=20,cex=1.5,col='red')
title(xlab='t0',ylab='s2',main='95% error ellipsoid on (t0,s2) plane')
# segments(b[1],b[2],u[1,1]+b[1],u[2,1]+b[2])
# segments(b[1],b[2],u[1,2]+b[1],u[2,2]+b[2])
range.t0=range(r[,1])
range.s2=range(r[,2])
abline(v=range.t0+b[1],col='blue',lty=2)
abline(h=range.s2+b[2],col='blue',lty=2)
axis(3,at=range.t0+b[1],labels=as.character(format(range.t0+b[1],digits=3)),col.axis='blue')
axis(4,at=range.s2+b[2],labels=as.character(format(range.s2+b[2],digits=3)),col.axis='blue')
## the conservative 95% interval for each parameter is shown in the plot as the bounding box for the error ellipsoid
print('The estimated value for t0:')
b[1]
print('with the 95% interval:')
range.t0+b[1]
print('The estimated value for s2:')
b[2]
print('with the 95% interval:')
range.s2+b[2]

## Evaluate the P-value for this model
dof=N-ncol(G)
chi2=sum( (tw-Gw %*% b)^2 )
p=1-pchisq(chi2,dof)
print('The p-value for the chi-square test is:')
p
## too small???

Nloop=1000
set.seed(1987)
t.random=rnorm(N*Nloop,mean=0,sd=sigma.t)
t.ran=matrix(t.random,ncol=N)
b.ran=matrix(NA,ncol=ncol(G),nrow=Nloop)
chi2.ran=array(dim=Nloop)

for (i in 1:Nloop) {
	t.temp=t.hat+t.ran[i,]
	tw.temp=W %*% t.temp
	b.ran[i,]= solve( t(Gw) %*% Gw ) %*% t(Gw) %*% tw.temp
	chi2.ran[i]=sum( (t.temp-G %*% b.ran[i,])^2 )/sigma.t^2
	}

range.expand<-function(lim,p=0.15) {
	len=max(lim)-min(lim)
	c(min(lim)-p*len,max(lim)+p*len)
}


library(ellipse)
plot(ellipse(cor.b,scale=sqrt(diag(cov.b)),centre=b,level=.95), type='l', xlim=range.expand(range(b[1]+r[,1])), ylim=range.expand(range(b[2]+r[,2])),lwd=2,col='red',ann=F,asp=1)
title(xlab='t0',ylab='s2',main='95% error ellipsoid on (t0,s2) plane')
points(b.ran[,1],b.ran[,2],pch=3,cex=.7,col='blue')
points(b[1],b[2],pch=20,cex=1.5,col='red')
legend('topright',legend=c('Monte Carlo simulations','Least-square estimation'),pch=c(3,20),col=c('blue','red'),cex=0.7)

hist(chi2.ran,freq=F)
temp.x=seq(0,40,by=0.1)
temp.y=dchisq(temp.x,dof)
lines(temp.x,temp.y,col='red')
legend('topright',legend='Theoretical chi^2 denstiy with DOF=4',lty=1,col='red')

## 1-norm estimation
b1=irls(Gw,tw,1e-10,1e-10,1,250)
print('1-norm solution to Gw %*% b1 = tw')
b1
## plot comparing the 2-norm and 1-norm estimation
t.hat1=G %*% b1
res1=t-t.hat1
plot(x,t,type='p',pch=15,ann=F)
title(xlab='Distance x/km',ylab='Arrival time t/s',main='Seismic data and fitted linear models')
lines(x,t.hat,col='red',lty=2,lwd=2)
lines(x,t.hat1,col='blue',lty=3,lwd=2)
legend('topleft',legend=c('Prediction from 1-norm estimation','Prediction from 2-norm estimation'),lty=c(3,2),col=c('blue','red'),lwd=2)
## residual comparison
res=t(cbind(res2,res1))
barplot(res,beside=T,col=c('red','darkblue'),ylab='Residual',xlab='# of observation',main='Residual from two models')
legend('topleft',legend=c('Residual for the 2-norm estimation','Residual for the 1-norm estimation'),cex=.8,pch=15,col=c('red','darkblue'))

##
Nloop1=1000
set.seed(1987)
t.random=rnorm(N*Nloop1,mean=0,sd=sigma.t)
t.ran=matrix(t.random,ncol=N)
b1.ran=matrix(NA,ncol=ncol(G),nrow=Nloop1)
b1.chi2.ran=array(dim=Nloop1)

for (i in 1:Nloop1) {
	t.temp=t.hat1+t.ran[i,]
	tw.temp=W %*% t.temp
	b1.ran[i,]=irls(Gw,tw.temp,1e-10,1e-10,1,250)
	b1.chi2.ran[i]=sum( (t.temp-G %*% b1.ran[i,])^2 )/sigma.t^2
}

b1.ave=matrix( rep(apply(b1.ran,2,mean),Nloop1), ncol=ncol(G),byrow=T)
A=b1.ran-b1.ave
cov.b1=t(A) %*% A /Nloop1

## 95% confidence interval estimated from covariance matrix
print('95% parameter confidence intervals (b1-, b1.est, b1+) on 1-norm solution')
del = 1.96*sqrt(diag(cov.b1));
cbind(b1-del , b1 , b1+del)

# error ellipsoid
plot(ellipse(cov.b1,centre=b1,level=.95), type='l', xlim=range(b1.ran[,1]), ylim=range(b1.ran[,2]),lwd=2,col='red',ann=F)
title(xlab='t0',ylab='s2',main='95% error ellipsoid on (t0,s2) plane')
points(b1.ran[,1],b1.ran[,2],pch=3,cex=.7,col='blue')
points(b1[1],b1[2],pch=20,cex=1.5,col='red')
legend('topright',legend=c('Monte Carlo simulations (1st-norm)','Least-square estimation'),pch=c(3,20),col=c('blue','red'),cex=0.7)


## histogram for b1.chi2.ran
hist(b1.chi2.ran,freq=F)
dof1=N-2
hist(b1.ran[,1])
hist(b1.ran[,2])

misfit1=G%*%b1-t
misfit1.tol=sum(abs(misfit1))
barplot(as.vector(abs(misfit1)/misfit1.tol),col='darkblue',ylab='1st-norm misfit (%)')
print('Contribution to the 1st-norm misfit for each data point is: (%)')
abs(misfit1)/misfit1.tol*100