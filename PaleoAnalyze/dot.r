dot<-function(sg1,sg2) {
  ## R function for calculating dot product between two signals
  # 1. standardize signals; 2. dot product; 3. normalize the results
  mean1=mean(sg1)
  mean2=mean(sg2)
  sd1=sd(sg1)
  sd2=sd(sg2)
  sg1n=(sg1-mean1)/sd1
  sg2n=(sg2-mean2)/sd2
  
  dot=sum(sg1n*sg2n)
  dotn=2*dot/(sum(sg1n^2)+sum(sg2n^2))

  ## method 2
  m=sqrt(sum(sg1^2))*sqrt(sum(sg2^2))
  n=sum(sg1*sg2)
  c=n/m
  theta=acos(c)*180/pi
  
  return(list(dotn,c,theta))
}

angle<-function(t,sg1,sg2,wd,olap=.5,n=10,PLOT=T,...) {
## R function calculating dot product with a moving window
t1=floor(wd/2)
dt=floor(wd*(1-olap))
ind=seq(1,length(sg1)-wd,by=dt)
tind=ind+t1

angle=vector(length=length(ind))
for (i in 1:length(ind)) {
	angle[i]=dot(sg1[ind[i]:(ind[i]+wd)],sg2[ind[i]:(ind[i]+wd)])[[3]]
	}

tt=t[tind]

if (PLOT) {
	quartz(width=16,height=9)
	yrange=c(min(c(scale(sg1),scale(sg2))),max(c(scale(sg1),scale(sg2))))
	plot(t,scale(sg1),type='l',ylim=c(yrange[1],yrange[2]*1.05),ann=F,axes=F)
	lines(t,scale(sg2),lty=2)
	axis(1,pretty(tt))
	axis(2,pretty(yrange))
	title(xlab='Years BP',ylab='',...)
	rect(min(tt),yrange[1],max(tt),yrange[2],border='grey')
	y=yrange[2]*1.05
	twd=wd/length(sg1)*(max(t)-min(t))
	x1=tind-floor(wd/4)
	x2=tind+floor(wd/4)
	y1=rep(y,length(x1))
	y2=y1
	cr=c(90-n,90+n)
	scol=vector(mode='character',length=length(x1))
	scol[which(angle>cr[2] | angle<cr[1])]='red'
	scol[which(angle<=cr[2] & angle>=cr[1])]='green'
	segments(t[x1],y1,t[x2],y2,col=scol,lwd=3)
}

return(data.frame(t=tt,y=angle))
}
