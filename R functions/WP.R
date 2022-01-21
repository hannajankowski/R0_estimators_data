
WP	<-	function(NT, mu="NA", method="unknown", search=list(B=100, shape.max=10, scale.max=10), tol=0.999){
	
	if(method=="unknown"){
		
		print("You have assumed that the serial distribution is unknown.")	
		res				<-	WP_unknown(NT=NT, B=search$B, shape.max=search$shape.max, scale.max=search$scale.max, tol=tol)
		Rhat			<-	res$Rhat
		p				<-	res$p
		range.max		<-	res$range.max
		JJ				<-	res$JJ
		
		
	}
	
	if(method=="known"){
		
		if(mu=="NA"){

			res		<-	"NA"
			print("For method=known, the mean of the serial distribution must be specified.")
						
		} else {
			
		print("You have assumed that the serial distribution is known.")	
		
		range.max	<-	ceiling(qexp(tol, rate=1/mu))
		p			<-	diff(pexp(0:range.max, 1/mu))
		p			<-	p/sum(p)
		res			<-	WP_known(NT=NT, p=p)
		Rhat		<-	res$Rhat
		JJ			<-	NA
		}
		
	}
	
return(list(Rhat=Rhat, check=length(JJ), SD=list(supp=1:range.max, pmf=p), inputs=list(NT=NT, mu=mu, method=method, search=search, tol=tol)))
	
}