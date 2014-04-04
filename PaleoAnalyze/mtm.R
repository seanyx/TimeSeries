mtm <- function(signal,dt,kind,nwin,npi,inorm,...) {
	library(RSEIS)
	y=signal-mean(signal)
	nn=next2(length(y))
	Mspec=mtapspec(y,dt,klen=nn,MTP=list(kind=kind,nwin=nwin,npi=npi,inorm=inorm))
	f=Mspec$freq
	amp=Mspec$spec[1:length(f)]

	#squig=list(y=y,dt=dt)
	#conf=autoreg(squig,numf=length(Mspec$freq),pord=1,PLOT=F)
	#AutoR=log10(conf$amp)

	f=f[-1]
	amp=amp[-1]

	return(list(f,amp))
}
