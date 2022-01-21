WP_known		<-	function(NT, p){
	
	k		<-	length(p)
	TT		<-	length(NT)-1
	mu_t	<-	rep(0, TT)
	for(i in 1:TT){
		Nt		<-	NT[i:max(1,i-k+1)]
#		print(Nt)
#		print(p[1:min(k,i)])
		mu_t[i]	<-	sum(p[1:min(k,i)]*Nt)	
	}
	Rhat	<-	sum(NT[-1])/sum(mu_t)
	return(list(Rhat=Rhat))
		
}


