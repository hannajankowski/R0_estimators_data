computeLL	<-	function(p, NT, R0){

	k		<-	length(p)
	TT		<-	length(NT)-1
	mu_t	<-	rep(0, TT)
	for(i in 1:TT){
		Nt		<-	NT[i:max(1,i-k+1)]
#		print(Nt)
#		print(p[1:min(k,i)])
		mu_t[i]	<-	sum(p[1:min(k,i)]*Nt)	
	}
	mu_t	<-	R0*mu_t
	LL		<-	-sum(mu_t)+sum(NT[-1]*log(mu_t))

	return(LL)
	
}


