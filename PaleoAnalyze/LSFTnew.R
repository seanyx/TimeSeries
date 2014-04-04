










#### from "spectral and cross-spectral analysis of uneven time series with the smoothed LS periodogram and Monte Carlo evaluation of statistical significance"

LSFT<-function(signal,freq) {
# given:
# signal #dataframe
# freq #frequency vector
	if (is.list(signal)) signal=data.frame(cbind(signal[[1]],signal[[2]]))
I=vector(mode='numeric')
z=signal[,2]
t=signal[,1]

mz=mean(z)
s2=var(z)
w=2*pi*freq

Nf=length(freq)
for (i in 1:Nf) {
	SA=sum(sin(2*w[i]*t))
	SC=sum(cos(2*w[i]*t))
	tao=atan(SA/SC)/2/w[i]
	
	SS=sum((sin(w[i]*(t-tao)))^2)
	CC=sum((cos(w[i]*(t-tao)))^2)
	
	ZS=sum((z-mz)*sin(w[i]*(t-tao)))
	ZC=sum((z-mz)*cos(w[i]*(t-tao)))
	
	I[i]=1/(2*s2)*(ZC^2/CC+ZS^2/SS)
}
I[1]=0
data.frame(f=freq,I=I)
}


SmoothLSFT<-function(LSFT,m,method) {
	# method: arithmetic mean OR triangular (or Bartlett) mean
	I=LSFT$I
	NI=length(I)
	Nw=2*m+1
	
	Is=vector(mode='numeric')
	
	if (method=='arithmetic') lamda=rep(1/Nw,length=Nw)
	if (method=='bartlett') lamda=c(seq(0,1/m,length=(m+1)),rev(seq(0,1/m,length=(m+1)))[-1])
	
	for (j in 1:NI) {
		k=seq(max(1,j-m),min(NI,j+m))
		l=seq(1,Nw)
		if ((j-m)<1) l=l[(j-m-1):(-1)]
		if ((j+m)>NI) l=l[1:(Nw-(j+m-NI))]
		Is[j]=sum(lamda[l]*I[k])
	}
	Is=data.frame(f=LSFT$f,Is=Is)
}


gSmooth<-function(y,m,method) {
	I=y
	NI=length(I)
	Nw=2*m+1
	
	Is=vector(mode='numeric')
	
	if (method=='arithmetic') lamda=rep(1/Nw,length=Nw)
	if (method=='bartlett') lamda=c(seq(0,1/m,length=(m+1)),rev(seq(0,1/m,length=(m+1)))[-1])
	
	for (j in 1:NI) {
		k=seq(max(1,j-m),min(NI,j+m))
		l=seq(1,Nw)
		if ((j-m)<1) l=l[(j-m-1):(-1)]
		if ((j+m)>NI) l=l[1:(Nw-(j+m-NI))]
		Is[j]=sum(lamda[l]*I[k])
	}
	Is
}

Permutation<-function(signal,freq,smooth=T,m,method,N) {
	I=LSFT(signal,freq)
	Is=SmoothLSFT(I,m,method)
	y=signal$y
	t=signal$t
	Ny=length(y)
	Nf=length(freq)
	J=matrix(0,nrow=Nf,ncol=N)
	
	for (i in 1:N) {
	yi=sample(y,size=Ny)
	yi=data.frame(t=t,y=yi)
	Ii=LSFT(yi,freq)
	Isi=SmoothLSFT(Ii,m,method)
	j=rep(0,Nf)
	j[which(Isi$Is>=Is$Is)]=1
	J[,i]=j
	}
	
	ASL=apply(J,1,mean)
	ACL=(1-ASL)*100
	
	data.frame(f=freq,ACL=ACL)
	
}



CrossLS<-function(signal1,signal2,freq,smooth=T,m,method) {
	z1=signal1$y
	t1=signal1$t
	z2=signal2$y
	t2=signal2$t
	
	mz1=mean(z1)
	mz2=mean(z2)
	s21=var(z1)
	s22=var(z2)
	w=2*pi*freq
	c=vector(mode='numeric')
	q=c
	I1=c
	I2=c
	
	
	Nf=length(freq)
	
	for (i in 1:Nf) {
	
	###signal 1
	SA1=sum(sin(2*w[i]*t1))
	SC1=sum(cos(2*w[i]*t1))
	tao1=atan(SA1/SC1)/2/w[i]
	
	SS1=sum((sin(w[i]*(t1-tao1)))^2)
	CC1=sum((cos(w[i]*(t1-tao1)))^2)
	
	ZS1=sum((z1-mz1)*sin(w[i]*(t1-tao1)))
	ZC1=sum((z1-mz1)*cos(w[i]*(t1-tao1)))
	
	a1=ZC1/sqrt(CC1)
	b1=ZS1/sqrt(SS1)
	
	###signal 2
	SA2=sum(sin(2*w[i]*t2))
	SC2=sum(cos(2*w[i]*t2))
	tao2=atan(SA2/SC2)/2/w[i]
	
	SS2=sum((sin(w[i]*(t2-tao2)))^2)
	CC2=sum((cos(w[i]*(t2-tao2)))^2)
	
	ZS2=sum((z2-mz2)*sin(w[i]*(t2-tao2)))
	ZC2=sum((z2-mz2)*cos(w[i]*(t2-tao2)))
	
	a2=ZC2/sqrt(CC2)
	b2=ZS2/sqrt(SS2)
	
	### cross statistics
	c[i]=a1*a2+b1*b2
	q[i]=a1*b2-a2*b1
	
	I1[i]=1/2*(a1^2+b1^2)
	I2[i]=1/2*(a2^2+b2^2)
	
	}
	
	if(smooth) {
		c=gSmooth(c,m,method)
		q=gSmooth(q,m,method)
	}
	
	alpha2=c^2+q^2
	C=alpha2/(I1*I2)/4
	phi=(atan(-q/c)+pi/2)/2/pi*360
	
	
	
	data.frame(f=freq,crossA=alpha2,Coh=C,Pha=phi)
}



# Gap filling and noise reduction of unevenly sampled data by means of the LS periodogram

complexLSFT<-function(signal,ofac=2,hifac=1) {
	## ofac: oversampling factor
	## hifac: integer, 1 for frequencies up the Nyquist frequency
	t=signal$t
	y=signal$y
	
	n=length(t)
	nout=.5*ofac*hifac*n
	
	tstart=t[1]
	t=t-tstart
	
	tave=.5*(max(t)+min(t))
	ave=mean(y)
	vari=var(y)
	
	tdif=max(t)-min(t)
	deltaf=1/(tdif*ofac)
	pnow=seq(1:nout)*deltaf
	
	arg=2*pi*((t-tave)*deltaf)
	wpr=-2*(sin(.5*arg))^2
	wpi=sin(arg)
	wr=cos(arg)
	wi=wpi
	yy=y-ave
	
	py=vector(mode='numeric')
	ph=vector(mode='numeric')
	ph1=vector(mode='numeric')
	teta=vector(mode='numeric')
		
	for (i in 1:nout) {
		sumsh=sum(wr*wi)
		sumc=sum((wr-wi)*(wr+wi))
		wtau=.5*atan2(2*sumsh,sumc)
		swtau=sin(wtau)
		cwtau=cos(wtau)
		ss=wi*cwtau-wr*swtau
		cc=wr*cwtau+wi*swtau
		
		sums=sum(ss^2)
		sumc=sum(cc^2)
		sumsy=sum(yy*ss)
		sumcy=sum(yy*cc)
		
		wtemp=wr
		wr=wr*wpr-wi*wpi+wr
		wi=wi*wpr+wtemp*wpi+wi
		
		iy=sumsy/sqrt(sums)
		ry=sumcy/sqrt(sumc)
		py[i]=.5*(ry^2+iy^2)/vari
		
		
		# here, the FFT phase is computed from the Lomb-Scargle Phase at each new frequency 'pnow' by adding the phase shift 'arg0'
		
		phLS=atan2(iy,ry)
		arg0=2*pi*(tave+tstart)*pnow[i]+wtau
		arg1=2*pi*tave*pnow[i]+wtau
		ph[i]=(phLS+arg0) %% (2*pi)
		ph1[i]=(phLS+arg1) %% (2*pi)
		teta[i]=arg1 %% (2*pi)
	}
	
	px=pnow
	dim=2*nout+1
	fac=sqrt(vari*dim/2)
	a=fac*sqrt(py)
	Fx=a*cos(ph1)
	Fy=a*sin(ph1)
	ph=(ph+5*2*pi) %% (2*pi)
	wk1=px
	wk2=py
	
# Fourier spectrum F: arrangement of the Lomb-Scargle periodogram   
# as a Fourier spectrum in the manner of Matlab:
# (it is not fully clear yet if and how the complex Fourier spectrum 
# can be exactly constructed from the Lomb-Scargle periodogram.  
# The present heuristic approach works well for the FFT back transformation   
# of F and reconstruction of an evenly spaced series 'yfit' in the time domain (after 
# multiplication by a constant, 'yfit' fits to 'y') 
	
	Fxr=rev(Fx)
	Fxr=Fxr[-1]
	Fyr=rev(Fy)
	Fyr=Fyr[-1]
	
# complex Fourier spectrum which corresponds to the Lomb-Scargle periodogram:
	FFT=c(complex(real=ave,imaginary=0),complex(real=Fx,imaginary=Fy),complex(real=Fxr,imaginary=-Fyr))
	
	FFT1=rev(c(complex(real=Fx,imaginary=Fy),complex(real=Fxr,imaginary=-Fyr),complex(real=ave,imaginary=0)))
	
	list(f=pnow,ifft=FFT,fft=FFT1,I=py,PH=ph)
}
	
	
	# I=vector(mode='numeric')
	# a=vector(mode='numeric')
	# b=vector(mode='numeric')
	# Aft=vector(mode='numeric')
	# Rft=vector(mode='numeric')
	# Ift=vector(mode='numeric')
	
	
	# z=signal$y
	# t=signal$t
	# Nz=length(z)
	
	# mz=mean(z)
	# s2=var(z)
	
	# freq=1/(t[Nz]-t[1])/ofac*seq(1,Nz*ofac/2)
	# w=2*pi*freq
	# t.ave=(t[1]+t[Nz])/2

	# Nf=length(freq)
	# for (i in 1:Nf) {
		# SA=sum(sin(2*w[i]*(t-t.ave)))
		# SC=sum(cos(2*w[i]*(t-t.ave)))
		# tao=atan(SA/SC)/2/w[i]
	
		# SS=sum((sin(w[i]*(t-tao-t.ave)))^2)
		# CC=sum((cos(w[i]*(t-tao-t.ave)))^2)
	
		# ZS=sum((z-mz)*sin(w[i]*(t-tao-t.ave)))
		# ZC=sum((z-mz)*cos(w[i]*(t-tao-t.ave)))
	
		# I[i]=1/(2*s2)*(ZC^2/CC+ZS^2/SS)
		
		# a[i]=sqrt(2/Nz)*ZC/sqrt(CC)
		# b[i]=sqrt(2/Nz)*ZS/sqrt(SS)		
		
		# ### phase part
		
		# Aft[i]=sqrt(a[i]^2+b[i]^2)
		
		
		# PHIft=-atan(b[i]/a[i])-w[i]*t.ave-tao
		
		# Rft[i]=Aft[i]*cos(PHIft)
		# Ift[i]=Aft[i]*sin(PHIft)
	# }
	
	# plxfft=c(complex(real=0,imaginary=0),complex(real=Rft,imaginary=Ift))
	
	# complexLS=list(f=freq,fft=c(plxfft,rev(Conj(plxfft))),Aft)
	
	# I=data.frame(f=freq,I=I)
	
	# return(list(I=I,complexLS=complexLS))
	
# }

plot.complexLS<-function(complexLS,doPower=T,doPhase=F,doSmooth=-1,...) {
	f=complexLS$f
	fft=complexLS$fft
	Nf=length(f)
	
	Power=(Mod(fft[1:Nf]))^2
	Phase=Arg(fft[1:Nf])
	
	if(doSmooth!=-1) {
		Power=gSmooth(Power,doSmooth,'arithmetic')
		Phase=gSmooth(Phase,doSmooth,'arithmetic')
	}
	
	if (doPower) {
		plot(f,Power,type='l',xlab='Frequency',ylab='Power',...)
		points(f,Power,cex=.2)
		}
		
	if (doPhase) {
		dev.new()
		plot(f,Phase/pi*180,type='l',xlab='Frequency',ylab='Phase (degree)',...)
		points(f,Phase,cex=.2)
	}
	
	invisible(data.frame(f=f,Power=Power,Phase=Phase))
}



### spectrogram for irregular sampled time series
LSspecgram<-function(signal,freq,wd,overlap=.5) {
	t=signal$t
	Nt=length(t)
	ns=floor(wd*overlap)+1
	no.start=seq(1,Nt,by=ns)
	Nf=length(freq)
	
	Ncol=(Nt-(Nt%%ns))/ns-1
	
	sgA=matrix(nrow=Nf,ncol=Ncol)
	sgP=sgA
	t.out=vector(length=Ncol)
	
	for (i in 1:Ncol) {
		sig=signal[(no.start[i]:(no.start[i]+wd-1)),]
		t.out[i]=(t[no.start[i]]+t[no.start[i]+wd-1])/2
		SIG=LSFT(sig,freq)
		sgA[,i]=SIG$I
		#sgP[,i]=SIG$
	}
  
	invisible(list(x=t.out,y=freq,t(sgA)))
}