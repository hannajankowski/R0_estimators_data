
WP_unknown		<-	function(NT, B=100, shape.max=10, scale.max=10, tol=0.999){
	
	shape		<-	seq(0, shape.max, length.out=B+1)
	scale		<-	seq(0, scale.max, length.out=B+1)
	shape		<-	shape[-1]
	scale		<-	scale[-1]
		
	resLL		<-	matrix(0,B,B)
	resR0		<-	matrix(0,B,B)
		
	for(i in 1:B){
	for(j in 1:B){	
		range.max	<-	ceiling(qgamma(tol, shape=shape[i], scale=scale[j]))
		p			<-	diff(pgamma(0:range.max, shape=shape[i], scale=scale[j]))
		p			<-	p/sum(p)
		mle			<-	WP_known(NT, p)
		resLL[i,j]	<-	computeLL(p, NT, mle$R)
		resR0[i,j]	<-	mle$R
	}		
#		print(i)
	}
	
	J0			<-	which.max(resLL)
	R0hat		<-	resR0[J0]	
	JJ			<-	which(resLL==resLL[J0], arr.ind=TRUE)
#	JJ			<-	which(resLL==max(resLL), arr.ind=TRUE)
	range.max	<-	ceiling(qgamma(tol, shape=shape[JJ[1]], scale=scale[JJ[2]]))
	p			<-	diff(pgamma(0:range.max, shape=shape[JJ[1]], scale=scale[JJ[2]]))
	p			<-	p/sum(p)
	
	return(list(Rhat=R0hat, J0=J0, ll=resLL, Rs=resR0, scale=scale, shape=shape, JJ=JJ, p=p, range.max=range.max))
}


