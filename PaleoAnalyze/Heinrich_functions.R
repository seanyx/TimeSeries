addHE<-function(range,color,ird=T,axis3=T,age_model='default') {
		# heinrich timing #1
	heinrich1=c(11600,14500,23300,28900,38200,45200,52100,58200,70000,74000,82500) # lower bound for heinrich events (0,1,2,3,4,5,5a,6)  adapted from Rashid 2003 on gisp time scale. The 7,8,9 from Dr. Rial's plot
	heinrich2=c(12500,15900,24100,30600,39100,46000,52900,59200,72000,75500,84000)
	
	if ( age_model == 'default' ) {
	# heinrich timing #2
	heinrich1=c(11500,14500,22000,27000,36000,45000,52500,57500,70000,74000,82500) # all read from Dr. Rial's ArcticWarmDomec.ai file
	heinrich2=c(12500,15500,22500,29000,37500,46500,54000,59000,72000,75500,84000)
	}
	
	if ( age_model == 'aicc2012') {
		
		heinrich1=c(11000,14000,23000,29000,37500,46000,53000,59000,72000,75500,84000)
		heinrich2=c(13000,18000,24500,30500,39500,47500,55000,62000,73000,77000,87000)
		
	}
	
	heinrich_ran=cbind(heinrich1,heinrich2)
	heinrich=apply(heinrich_ran,1,mean)
	hname=c('H0','H1','H2','H3','H4','H5','H5a','H6','H7a','H7b','H8')
	# lcol=rainbow(n1*n2,v=.5,alpha=.5)

	# IRDs
	irds=c(23500,26000,30500,31800,33500,39200,41500)
	irdn=c('d','e','f','g','h','j','k')

	for ( i in 1:length(heinrich)) polygon(c(rep(heinrich_ran[i,1],2),rep(heinrich_ran[i,2],2)),c(range[1],rep(range[2],2),range[1]),border=F,col=color,bg='white') # plot the heinrich events with range
	if (axis3 & ird) axis(3,at=c(heinrich,irds),label=c(hname,irdn))
	if (axis3 & !ird) axis(3,at=heinrich,label=hname)
	if (ird) abline(v=irds,lty=2,col='black')
}

standardize<-function(x) (x-mean(x))/sd(x)
Tdiff<-function(s,n) (s-n)+abs(s-n)
Eest<-function(s,n) s^2+n^2

reconEMD <- function(signal,tt=NULL,recon_imfs,plot=F) {
		
	## synthesize a signal through sum of its imfs after empirical mode decomposition
		
	# signal: vector of vector
	# tt: time index
	# recon <- imfs: index of imfs that will be used for reconstruction
		
	library(EMD)
	
	decomp=emd(signal,tt=tt)
	imfs=decomp$imf
	residue=decomp$residue
	nimf=decomp$nimf
			
	if (max(recon_imfs) > nimf) stop('Reconstructing imfs exceed number of imfs generated!')
				
	recon=apply(imfs[,recon_imfs],1,sum)
			
	print(paste('EMD results in',nimf,'IMFs. The imfs used for reconstruction is',min(recon_imfs),'-',max(recon_imfs)))
					
	if(plot) {
		yrange=range(c(recon,signal))
		xrange=c(1,length(signal))
		plot(xrange,yrange,type='n',xlab='Points',ylab='Input signal (grey); Reconstrution (black)')
		lines(signal,col='grey')
		lines(recon,col='black',lwd=2)
		}
					
	return(recon)
}


innerproduct <- function(signal1,signal2,cos=F) {
	## calculate the inner product of signal1 and signal2
	
	# signal1, signal2: input signals, must be of the same length
		
	n1=length(signal1)
	n2=length(signal2)
	if (n1!=n2) stop('input signals are of different length!')
		
	scalarproduct=sum(signal1*signal2)
			
	if (!cos) return(scalarproduct)
	if (cos) {
		norm1=sqrt(sum(signal1^2))
		norm2=sqrt(sum(signal2^2))
		return(scalarproduct/(norm1*norm2))	
	}
				
}


mvsproduct <- function(signal1,signal2,wd,incre,cos=F,plot=T) {
		
	## calculate inner product of two signal (high dimensional vector) through a moving window with window length 'wd' and moving increment 'incre'
		
	# signal1, signal2: input signals, must be of the same length
	# wd: window width (in points)
	# incre: increments (in points) of the moving window
	# plot: optional plot for fast examining the results
		
	n1=length(signal1)
	n2=length(signal2)
	if (n1!=n2) stop('input signals are of different length!')
	if (wd > n1) stop('window length is larger than the length of the signal!')
			
	n=n1
	ind1=seq(1,n-wd+1,by=incre) # index of the beginning point for each window
	ind2=ind1+wd-1 # index of the end point for each window
	ind=ind1+floor(wd/2) # index of the mid point for each window
					
	mvsproduct=ind
					
	for (i in 1:length(ind)) {
		indtemp=ind1[i]:ind2[i]
		mvsproduct[i]=innerproduct(signal1[indtemp],signal2[indtemp],cos=cos)
	}
						
	if(plot) plot(ind,mvsproduct,type='l',xlab='midpoint index',ylab='moving window scalar product')
							
	return(list(ind,mvsproduct))
							
}
