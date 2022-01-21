seqB		<-	function(NT, mu, kappa=20){
	
	gamma	<-	1/mu
	
	if(length(NT)<2) {
		print("Warning: length of NT should be at least two.")
	}
	else{
	if(min(NT)>0){
		times		<-	1:length(NT)		
		tau			<-	diff(times)  
	}
	
	group			<-	FALSE
	
	if(min(NT)==0){
		times	<-	which(NT>0)
		NT		<-	NT[times]
		tau		<-	diff(times)
		group	<-	TRUE
	}
	
	R			<-	seq(0, kappa, 0.01)
	prior0		<-	rep(1,kappa/0.01+1)
	prior0		<-	prior0/sum(prior0)
	
	k			<-	length(NT)-1
	R0.post		<-	matrix(0, nrow=k, ncol=length(R))
	
	prior		<-	prior0
	posterior	<-	seq(0, length(prior0))
	
	for(i in 1:k){
		
	mm1			<-	NT[i]			
	mm2			<-	NT[i+1]		
	lambda		<-	tau[i]*gamma*(R-1)
	lambda		<-	log(mm1)+lambda
	loglik		<-	mm2*lambda-exp(lambda)
			
# numerical issues solve???

#	if((mm1==0)*(mm2==0)==0){
	maxll			<-	max(loglik)
	const			<-	0
	if(maxll>700){
	const			<-	maxll-700		
		}	
	loglik			<-	loglik-const
#	}
	
# end numerical issues solve???
	
	posterior	<-	exp(loglik)*prior
	posterior	<-	posterior/sum(posterior)

	
	prior		<-	posterior
		
	}
	
	Rhat		<-	sum(R*posterior)
	
return(list(Rhat=Rhat, posterior=list(supp=R, pmf=posterior), group=group, inputs=list(NT=NT, mu=mu, kappa=kappa)))
}	
}


