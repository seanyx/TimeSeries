## RK4 solution for Van der Pol oscillator (Rial & Saha, 2011), no noise situation
## 

## import insolation
## link: ftp://ftp.ncdc.noaa.gov/pub/data/paleo/insolation/
inso=read.table("~/Documents/TimeSeries/ClimateModels/VDP/orbit91",skip=3)[1:101,6]
inso=rev(inso)
inso=(inso-mean(inso))/sd(inso)

f1=function(omega1,M1,u2,t,inso,DC,D) -omega1^2*u2+M1*(approx(1:101,inso,t/1000,rule=2)[[2]]+DC)+D*rnorm(1)
f2=function(miu1,u1,u2,u3,u4,q1,q2,u1prime,u3prime) miu1*(u2-1/3*u2^3)+u1+q1*(u1prime-u3prime)+q2*(u1-u3)
f3=function(omega2,M2,u4,t,inso,DC,D) -omega2^2*u4+M2*(approx(1:101,inso,t/1000,rule=2)[[2]]+DC)+D*rnorm(1)
f4=function(miu2,u1,u2,u3,u4,q1,q2,u1prime,u3prime) miu2*(u4-1/3*u4^3)+u3+q1*(u3prime-u1prime)+q2*(u3-u1)

## input parameters
## A and a together control the amount of frequency modulation
## b mostly controls the amplitude
## increase a will also increase the natural periodicity
fc1=1/5 # natural freqeuncy for northern oscillator
fc2=1/8 # natural freqeuncy for southern oscillator
omega1=2*pi*fc1
omega2=2*pi*fc2
M1=0 # forcing amplitude for northern oscillator
M2=0 # forcing amplitude for southern socillator
miu1=0 # positive coeff for ice albedo and greenhouse gas positive feedbacks, and the nonlinear term to control the positive feedback
miu2=0 # the same as a1
q1=0.1 # amplitude for the dissipative coupling term
q2=5 # amplitude for the reactive coupling term
D=0 # noise level
DC=0.6


h=0.1 # integration step
N=1000
u1=vector(length=N+1)
u2=u1
u3=u1
u4=u1

## initial conditions
u1[1]=0.05
u2[1]=0
u3[1]=0.05
u4[1]=0
for (i in 1:N) {
	
	u1temp=f1(omega1,M1,u2[i],h*(i-1),inso,DC,D)
	u3temp=f3(omega2,M2,u4[i],h*(i-1),inso,DC,D)
	#         u2temp=f2(miu1,u1[i],u2[i],u3[i],u4[i],q1,q2,u1temp,u3temp)
	#         u4temp=f4(miu2,u1[i],u2[i],u3[i],u4[i],q1,q2,u1temp,u3temp)

	m=h*c(f1(omega1,M1,u2[i],h*(i-1),inso,DC,D), f3(omega2,M2,u4[i],h*(i-1),inso,DC,D), 
	      f2(miu1,u1[i],u2[i],u3[i],u4[i],q1,q2,u1temp,u3temp),
	      f4(miu2,u1[i],u2[i],u3[i],u4[i],q1,q2,u1temp,u3temp))
	n=h*c(f1(omega1,M1,u2[i]+1/2*m[3],h*(i-1)+1/2*h,inso,DC,D), f3(omega2,M2,u4[i]+1/2*m[4],h*(i-1)+1/2*h,inso,DC,D), 
	      f2(miu1,u1[i]+1/2*m[1],u2[i]+1/2*m[3],u3[i]+1/2*m[2],u4[i]+1/2*m[4],q1,q2,u1temp,u3temp),
	      f4(miu2,u1[i]+1/2*m[1],u2[i]+1/2*m[3],u3[i]+1/2*m[2],u4[i]+1/2*m[4],q1,q2,u1temp,u3temp))
	j=h*c(f1(omega1,M1,u2[i]+1/2*n[3],h*(i-1)+1/2*h,inso,DC,D), f3(omega2,M2,u4[i]+1/2*n[4],h*(i-1)+1/2*h,inso,DC,D), 
	      f2(miu1,u1[i]+1/2*n[1],u2[i]+1/2*n[3],u3[i]+1/2*n[2],u4[i]+1/2*n[4],q1,q2,u1temp,u3temp),
	      f4(miu2,u1[i]+1/2*n[1],u2[i]+1/2*n[3],u3[i]+1/2*n[2],u4[i]+1/2*n[4],q1,q2,u1temp,u3temp))
	k=h*c(f1(omega1,M1,u2[i]+j[3],h*(i-1)+h,inso,DC,D), f3(omega2,M2,u4[i]+j[4],h*(i-1)+h,inso,DC,D), 
	      f2(miu1,u1[i]+j[1],u2[i]+j[3],u3[i]+j[2],u4[i]+j[4],q1,q2,u1temp,u3temp),
	      f4(miu2,u1[i]+j[1],u2[i]+j[3],u3[i]+j[2],u4[i]+j[4],q1,q2,u1temp,u3temp))

	u1[i+1]=u1[i]+1/6*(m[1]+2*n[1]+2*j[1]+k[1])
	u2[i+1]=u2[i]+1/6*(m[3]+2*n[3]+2*j[3]+k[3])
	u3[i+1]=u3[i]+1/6*(m[2]+2*n[2]+2*j[2]+k[2])
	u4[i+1]=u4[i]+1/6*(m[4]+2*n[4]+2*j[4]+k[4])
	
}

dev.new(); plot(u1,type='l',main='u13'); lines(u3,col='red')
dev.new(); plot(u3,type='l',main='u23'); lines(u2,col='red')

# para=paste("A=",A,"; a=",a,"; Tc=",1/fc,"; b=",b,"; c=",c,"; h=",h,"; N=",N,"; D=",D,'; DC=',DC,sep='')
# t=seq(0,h*N,by=h)

# dev.new(width=12,height=6)
# par(mfrow=c(2,1))
# plot(t,xn,type='l')
# lines(t,approx(1:101,inso+DC,t/1000,rule=2)[[2]]/20,col='red')
# abline(h=0,col='green')
# title(main=para)

# ngrip=read.table("~/Documents/TimeSeries/ClimateModels/VDP/ngrip-d18o-50yr.txt",header=F,skip=80,col.names=c('t','y'),sep='')
# ngrip=ngrip[,2]
# tt=seq(0,100000,by=50)
# ngrip=ngrip[(1:length(tt))*2]
# plot(rev(tt),ngrip,type='l',main='NGRIP')
# par(mfrow=c(1,1))
