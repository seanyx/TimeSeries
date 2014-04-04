matrix2array<-function(m,row,col) {

	nr=nrow(m)
	nc=ncol(m)
	t1=matrix(nrow=nr*nc,ncol=3)
	for (i in 1:nr) {
		t1[((i-1)*nc+1):(i*nc),]=cbind(col,rep(row[i],nc),m[i,])
	}
	return(t1)

}
